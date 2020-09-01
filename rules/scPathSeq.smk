from os.path import join, basename

CR_SAMPLE_ODIR = "{patient}-{sample}"
CR_BAM_FILE = join(CR_SAMPLE_ODIR, "outs", "possorted_genome_bam.bam")
CR_UNMAPPED_BAM_FILE = join(CR_SAMPLE_ODIR, "outs", "unmapped.bam")
CR_UNMAPPED_BAM_FILE = join(CR_SAMPLE_ODIR, "outs", "unmapped.bam.bai")
CR_UNMAPPED_TRIMMED_BAM = join(CR_SAMPLE_ODIR, "outs", "unmapped.trimmed.bam")
CR_UNMAPPED_TRIMMED_BAM = join(CR_SAMPLE_ODIR, "outs", "unmapped.trimmed.bam")
UNMAPPED_FQ1 = join(CR_SAMPLE_ODIR, "outs", "unmapped.trimmed.fq1")
UNMAPPED_FQ2 = join(CR_SAMPLE_ODIR, "outs", "unmapped.trimmed.fq2")

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
    script:
        "module load samtools && "
        "samtools sort -n {input} {output}"

rule convert_to_fastq:
    input:
        CR_UNMAPPED_TRIMMED_QNAME_SORTED_BAM
    output:
        UNMAPPED_FQ1,
        UNMAPPED_FQ2
    shell:
        "module load bedtools && "
        "bamToFastq -i {input} -fq {output[0]} -fq2 {output[1]}"


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
