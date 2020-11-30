from os.path import join, basename


# params
GATK_VERSION = "4.1.8.1"

rule PathSeqPipelineSpark:
    input:
        bam_file = config["PathSeq"]["bam_file"],
        host_bwa_image = config["PathSeq"]["host_img"],
        microbe_bwa_image = config["PathSeq"]["microbe_bwa_image"],
        microbe_dict_file = config["PathSeq"]["microbe_dict"],
        host_hss_file = config["PathSeq"]["host_bfi"],
        taxonomy_db = config["PathSeq"]["taxonomy_db"]
    params:
        host_bwa_image = basename(config["PathSeq"]["host_img"]),
        microbe_bwa_image = basename(config["PathSeq"]["microbe_bwa_image"]),
        microbe_dict_file = basename(config["PathSeq"]["microbe_dict"]),
        host_hss_file = basename(config["PathSeq"]["host_bfi"]),
        taxonomy_db = basename(config["PathSeq"]["taxonomy_db"])
    output:
        pathseq_bam = join("output", "PathSeq", "{patient}-{sample}-{plate}", "pathseq.bam"),
        pathseq_output = join("output", "PathSeq", "{patient}-{sample}-{plate}", "pathseq.txt"),
        filter_metrics = join("output", "PathSeq", "{patient}-{sample}-{plate}", "filter-metrics.txt"),
        score_metrics = join("output", "PathSeq", "{patient}-{sample}-{plate}", "score-metrics.txt"),
    benchmark:
        "benchmarks/{patient}-{sample}-{plate}.PathSeqPipelineSpark_host_filter_single.txt"
    run:
        shell("mkdir /lscratch/$SLURM_JOBID/tmp")
        shell("cp {input.host_bwa_image} /lscratch/$SLURM_JOBID/")
        shell("cp {input.microbe_bwa_image} /lscratch/$SLURM_JOBID/")
        shell("cp {input.microbe_dict_file} /lscratch/$SLURM_JOBID/")
        shell("cp {input.host_hss_file} /lscratch/$SLURM_JOBID/")
        shell("cp {input.taxonomy_db} /lscratch/$SLURM_JOBID/")
        shell(
            "module load GATK/4.1.8.1 && "
            "gatk PathSeqPipelineSpark "
            "--filter-duplicates false "
            "--min-score-identity .7 "
            "--input '{input.bam_file}' "
            "--filter-bwa-image /lscratch/$SLURM_JOBID/{params.host_bwa_image} "
            "--kmer-file /lscratch/$SLURM_JOBID/{params.host_hss_file} "
            "--microbe-bwa-image /lscratch/$SLURM_JOBID/{params.microbe_bwa_image} "
            "--microbe-dict /lscratch/$SLURM_JOBID/{params.microbe_dict_file} "
            "--taxonomy-file /lscratch/$SLURM_JOBID/{params.taxonomy_db} "
            "--output '{output.pathseq_bam}' "
            "--scores-output '{output.pathseq_output}' "
            "--filter-metrics '{output.filter_metrics}' "
            "--score-metrics '{output.score_metrics}' "
            '--java-options "-Xmx96g -Xms96G -Djava.io.tmpdir=/lscratch/$SLURM_JOBID/tmp -XX:+UseG1GC -XX:ParallelGCThreads=8 -XX:ConcGCThreads=2" '
            '--spark-master local[8] '
            + config["params"]["PathSeq"]
        )

# -v means to only report those entries in A that have no overlap in B
rule identify_reads_with_vector_contamination:
    group:
        "scPathSeq"
    input:
        join("output", "PathSeq", "{patient}-{sample}-{plate}", "pathseq.bam"),
        "/data/Robinson-SB/run-VecScreen/output/microbev1-vecscreen-combined-matches.bed"
    output:
        temp(join("output", "PathSeq", "{patient}-{sample}-{plate}", "pathseq.contaminants.bam")),
    shell:
        "module load bedtools && "
        "bedtools intersect -abam {input[0]} -b {input[1]} > {output[0]}"

rule get_query_names_for_vector_contaminants:
    input:
        join("output", "PathSeq", "{patient}-{sample}-{plate}", "pathseq.contaminants.bam")
    output:
        join("output", "PathSeq", "{patient}-{sample}-{plate}", "contaminants.qname.txt")
    shell:
        "module load samtools && "
        "samtools view {input} | cut -f 1 > {output}"

rule filter_vector_contaminant_reads:
    input:
        join("output", "PathSeq", "{patient}-{sample}-{plate}", "pathseq.bam"),
        join("output", "PathSeq", "{patient}-{sample}-{plate}", "contaminants.qname.txt")
    output:
        join("output", "PathSeq", "{patient}-{sample}-{plate}", "pathseq.filtered.bam")
    shell:
        "module load picard && "
        "java -jar $PICARDJARPATH/picard.jar FilterSamReads "
        "I={input[0]} O={output} READ_LIST_FILE={input[1]}"
        "FILTER=excludeReadList"

rule split_PathSeq_BAM_by_RG:
    group:
        "scPathSeq"
    input:
        pathseq_bam = join("output", "PathSeq", "{patient}-{sample}-{plate}", "pathseq.filtered.bam"),
    output:
        temp(join("output", "PathSeq", "{patient}-{sample}-{plate}-{cell}", "pathseq.bam")),
    shell:
        "module load samtools && "
        "samtools view -h -b -r {wildcards.cell} {input} > {output}"

rule extract_paired_reads:
    group:
        "scPathSeq"
    input:
        join("output", "PathSeq", "{patient}-{sample}-{plate}-{cell}", "pathseq.bam"),
    output:
        temp(join("output", "PathSeq", "{patient}-{sample}-{plate}-{cell}", "pathseq.paired.bam")),
    shell:
        "module load samtools && "
        "samtools view -h -b -f 1 {input} > {output}"

rule sort_paired_reads:
    group:
        "scPathSeq"
    input:
        join("output", "PathSeq", "{patient}-{sample}-{plate}-{cell}", "pathseq.paired.bam"),
    output:
        temp(join("output", "PathSeq", "{patient}-{sample}-{plate}-{cell}", "pathseq.paired.sorted.bam")),
    shell:
        "module load samtools && "
        "samtools sort -n -o {output} {input} "

rule extract_unpaired_reads:
    group:
        "scPathSeq"
    input:
        join("output", "PathSeq", "{patient}-{sample}-{plate}-{cell}", "pathseq.bam"),
    output:
        temp(join("output", "PathSeq", "{patient}-{sample}-{plate}-{cell}", "pathseq.unpaired.bam"))
    shell:
        "module load samtools && "
        "samtools view -h -b -F 1 {input} > {output}"

rule score_PathSeq_cell_BAM:
    group:
        "scPathSeq"
    input:
        paired_bam = join("output", "PathSeq", "{patient}-{sample}-{plate}-{cell}", "pathseq.paired.sorted.bam"),
        unpaired_bam = join("output", "PathSeq", "{patient}-{sample}-{plate}-{cell}", "pathseq.unpaired.bam"),
        taxonomy_db = config["PathSeq"]["taxonomy_db"]
    output:
        pathseq_output = join("output", "PathSeq", "{patient}-{sample}-{plate}-{cell}", "pathseq.txt"),
    shell:
        "module load GATK/4.1.8.1 && "
        "gatk PathSeqScoreSpark "
        "--min-score-identity .7 "
        "--unpaired-input '{input.unpaired_bam}' "
        "--paired-input '{input.paired_bam}' "
        "--taxonomy-file {input.taxonomy_db} "
        "--scores-output '{output.pathseq_output}' "
        '--java-options "-Xmx5g -Xms5G -XX:+UseG1GC -XX:ParallelGCThreads=2 -XX:ConcGCThreads=2" '
        "--conf spark.port.maxRetries=64 "
        '--spark-master local[2] ' + config["params"]["PathSeqScore"]
