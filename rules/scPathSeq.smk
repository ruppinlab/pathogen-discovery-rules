from os.path import join, basename


PATHSEQ_BAM = join("output", "PathSeq", "{patient}-{sample}", "pathseq.bam")
PATHSEQ_TAG_BAM = join("output", "PathSeq", "{patient}-{sample}", "pathseq_with_tags.bam")
PATHSEQ_TAG_BAI = join("output", "PathSeq", "{patient}-{sample}", "pathseq_with_tags.bam.bai")
PATHSEQ_CELL_BAM = join("output", "PathSeq", "{patient}-{sample}-{cell}", "pathseq_with_tags.bam")
PATHSEQ_CELL_SCORE = join("output", "PathSeq", "{patient}-{sample}-{cell}", "pathseq.txt")

localrules: PathSeqScoreSpark, split_PathSeq_BAM_by_CB_UB

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
