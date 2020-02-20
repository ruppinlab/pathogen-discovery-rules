include: "PathSeq.smk"

SC_PATHSEQ_OUTPUT = join("output", "scPathSeq", "{patient}.{sample}")

FILTERED_FASTA_FILE = join(SC_PATHSEQ_OUTPUT, "CB.UB.STAR.filtered.fasta")
FILTERED_PATHSEQ_FASTA_FILE = join(SC_PATHSEQ_OUTPUT, "PathSeq.filtered.fasta")
PATHSEQ_QNAME_FILE = join(SC_PATHSEQ_OUTPUT, "QNAME.PathSeq.tsv")
PATHSEQ_YP_FILE = join(SC_PATHSEQ_OUTPUT, "YP.PathSeq.tsv")
STAR_QNAME_FILE = join(SC_PATHSEQ_OUTPUT, "QNAME.STAR.tsv")
STAR_CB_FILE = join(SC_PATHSEQ_OUTPUT, "CB.STAR.tsv")
STAR_UB_FILE = join(SC_PATHSEQ_OUTPUT, "UB.STAR.tsv")
QNAME_CB_UB_FILE = join(SC_PATHSEQ_OUTPUT, "QNAME.CB.UB.STAR.tsv")
PATHSEQ_QNAME_YP_FILE = join(SC_PATHSEQ_OUTPUT, "QNAME.YP.PathSeq.tsv")
FILTERED_PATHSEQ_READS_FILE = join(SC_PATHSEQ_OUTPUT, "Filtered_PathSeq_reads.tsv")

rule combine_PathSeq_STAR:
    conda:
        "../envs/pandas.yml"
    input:
        qname = STAR_QNAME_FILE,
        cb = STAR_CB_FILE,
        ub = STAR_UB_FILE,
        path_qname = PATHSEQ_QNAME_FILE,
        yp = PATHSEQ_YP_FILE
    output:
        QNAME_CB_UB_FILE,
        PATHSEQ_QNAME_YP_FILE,
        FILTERED_PATHSEQ_READS_FILE
    script:
        "../src/combine_PathSeq_STAR.py"

rule extract_QNAME_YP_from_PathSeq_FASTA:
    input:
        FILTERED_PATHSEQ_FASTA_FILE
    output:
        PATHSEQ_QNAME_FILE,
        PATHSEQ_YP_FILE
    shell:
        "cat {input} | cut -f 1 > {output[0]} && "
        "cat {input} | grep -E -o 'YP:Z:([0-9]+,?)+' > {output[1]}"

rule extract_QNAME_from_BAM_FILE:
    input:
        FILTERED_FASTA_FILE
    output:
        STAR_QNAME_FILE,
        STAR_CB_FILE,
        STAR_UB_FILE
    shell:
        "cat {input} | cut -f 1 > {output[0]} && "
        "cat {input} | grep -E -o 'CB:Z:[G,A,T,C]+' > {output[1]} && "
        "cat {input} | grep -E -o 'UB:Z:[G,A,T,C]+' > {output[2]}"

# filter only the reads that PathSeq mapped to microbial genomes
rule filter_PathSeq_bam:
    input:
        join("output", "PathSeq", "{patient}.{sample}.pathseq.bam")
    output:
        FILTERED_PATHSEQ_FASTA_FILE
    shell:
        "module load samtools && "
        "samtools view {input} | grep -E 'YP:Z:*' > {output}"

# we only want reads with UB and CB tags
rule filter_BAM_FILE:
    input:
        join("output", "STARsolo", "{patient}.{sample}", "Aligned.sortedByCoord.out.bam")
    output:
        FILTERED_FASTA_FILE
    shell:
        "module load samtools && "
        "samtools view {input} | grep -E 'UB:Z:*' | grep -E 'CB:Z:*' > {output}"
