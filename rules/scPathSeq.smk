from os.path import join, basename

CR_SAMPLE_ODIR = "{patient}-{sample}"
CR_BAM_FILE = join(CR_SAMPLE_ODIR, "outs", "possorted_genome_bam.bam")
CR_UNMAPPED_BAM_FILE = join(CR_SAMPLE_ODIR, "outs", "unmapped.bam")
CR_UNMAPPED_BAI_FILE = join(CR_SAMPLE_ODIR, "outs", "unmapped.bam.bai")
CR_UNMAPPED_TRIMMED_BAM = join(CR_SAMPLE_ODIR, "outs", "unmapped.trimmed.bam")
CR_UNMAPPED_TRIMMED_QNAME_SORTED_BAM = join(CR_SAMPLE_ODIR, "outs", "unmapped.trimmed.qname.sorted.bam")
UNMAPPED_FQ1 = join("FASTQ", "unmapped", "polyA_TSO_removed" "{patient}-{sample}_1.fastq.gz")
TRIMMED_FQ1 = join("FASTQ", "unmapped", "trimmed" "{patient}-{sample}_1.fastq.gz")
FAILED_READS_FILE = join("FASTQ", "unmapped", "trimmed" "{patient}-{sample}_failed.fastq.gz")
FASTP_JSON_REPORT = join("FASTQ", "unmapped", "trimmed" "{patient}-{sample}-report.json")
FASTP_HTML_REPORT = join("FASTQ", "unmapped", "trimmed" "{patient}-{sample}-report.html")
UNALIGNED_BAM = join("output", "BAM", "{patient}-{sample}-unaligned.bam")

PATHSEQ_BAM = join("output", "PathSeq", "{patient}-{sample}", "pathseq.bam")
PATHSEQ_TAG_BAM = join("output", "PathSeq", "{patient}-{sample}", "pathseq_with_tags.bam")
PATHSEQ_TAG_BAI = join("output", "PathSeq", "{patient}-{sample}", "pathseq_with_tags.bam.bai")
PATHSEQ_CELL_BAM = join("output", "PathSeq", "{patient}-{sample}-{cell}", "pathseq_with_tags.bam")
PATHSEQ_CELL_SCORE = join("output", "PathSeq", "{patient}-{sample}-{cell}", "pathseq.txt")

localrules: PathSeqScoreSpark, split_PathSeq_BAM_by_CB_UB, filter_aligned_reads, trim_reads, sort_by_query_name, convert_to_fastq


rule filter_aligned_reads:
    input:
        CR_BAM_FILE
    output:
        CR_UNMAPPED_BAM_FILE,
        CR_UNMAPPED_BAI_FILE
    shell:
        "module load samtools && "
        "samtools view -h -b -f 4 {input} > {output[0]} && "
        "samtools index {output[0]}"

rule trim_reads:
    conda:
        "../envs/pysam.yaml"
    input:
        CR_UNMAPPED_BAM_FILE,
        CR_UNMAPPED_BAI_FILE
    output:
        CR_UNMAPPED_TRIMMED_BAM
    script:
        "../src/trim_TSO_polyA_sequences.py"

rule sort_by_query_name:
    input:
        CR_UNMAPPED_TRIMMED_BAM
    output:
        CR_UNMAPPED_TRIMMED_QNAME_SORTED_BAM
    shell:
        "module load samtools && "
        "samtools sort -n -o {output} {input}"

rule convert_to_fastq:
    input:
        CR_UNMAPPED_TRIMMED_QNAME_SORTED_BAM
    output:
        UNMAPPED_FQ1
    shell:
        "module load bedtools && "
        "bamToFastq -i {input} -fq {output[0]}"

rule run_fastp:
    conda:
        join(ENV_DIR, "fastp.yml")
    input:
        UNMAPPED_FQ1
    output:
        TRIMMED_FQ1,
        FAILED_READS_FILE,
        FASTP_JSON_REPORT,
        FASTP_HTML_REPORT
    threads:
        6
    shell:
        "fastp -w {threads} "
        "--unqualified_percent_limit 40 " # filter reads where 40% of bases have phred quality < 15
        "--cut_tail " # use defaults --cut_window_size 4 --cut_mean_quality 20
        "--low_complexity_filter " # filter reads with less than 30% complexity (30% of the bases are different from the preceeding base)
        "--length_required 25 "
        "--disable_adapter_trimming " # we already trimmed polyA tail and TSO
        "-i {input[0]} -o {output[0]} --failed_out {output[1]} "
        "-j {output[2]} -h {output[3]}"

rule FastqToBam:
    input:
        TRIMMED_FQ1
    output:
        UNALIGNED_BAM
    shell:
        "module load picard && "
        "java -jar $PICARDJARPATH/picard.jar FastqToSam "
        "F1={input} O={output} "
        "SM={wildcards.patient}.{wilcards.sample} "
        "RG={wildcards.patient}.{wilcards.sample} "


rule PathSeqScoreSpark:
    input:
        bam_file = PATHSEQ_CELL_BAM,
        taxonomy_db = config["PathSeq"]["taxonomy_db"]
    output:
        pathseq_output = PATHSEQ_CELL_SCORE
    run:
        shell(
            "module load GATK/4.1.6.0 && "
            "gatk PathSeqScoreSpark "
            "--unpaired-input '{input.bam_file}' "
            "--taxonomy-file {input.taxonomy_db} "
            "--scores-output '{output.pathseq_output}' "
            '--java-options "-Xmx30g -Xms30G -XX:+UseG1GC -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2" '
            '--spark-master local[2] ' + config["params"]["PathSeqScore"]
        )



rule split_PathSeq_BAM_by_CB_UB:
    conda:
        "../envs/pysam.yaml"
    input:
        PATHSEQ_TAG_BAM,
        PATHSEQ_TAG_BAI
    output:
        PATHSEQ_CELL_BAM
    script:
        "../src/split_PathSeq_BAM_by_CB_UB.py"

rule add_CB_UB_tags_to_PathSeq_BAM:
    conda:
        "../envs/pysam.yaml"
    input:
        CR_BAM_FILE,
        PATHSEQ_BAM
    output:
        PATHSEQ_TAG_BAM,
    script:
        "../src/add_tags_to_PathSeq_bam.py"

rule index_PathSeq_BAM_with_tags:
    input:
        PATHSEQ_TAG_BAM
    output:
        PATHSEQ_TAG_BAI
    shell:
        "module load samtools && "
        "samtools index {input}"
