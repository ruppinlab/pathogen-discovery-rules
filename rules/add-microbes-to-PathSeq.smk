from os.path import join

NCBI_REFSEQ_FTP = "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/"
BACTERIA_REFSEQ_URL = join(NCBI_REFSEQ_FTP, "bacteria")

BACTERIA_FNA_FILE = join("raw", "bacteria", "bacteria.{fn}.1.genomic.fna")
BACTERIA_TO_ADD_FA_FILE = join("data", "{microbe}_to_add.fa")
PATHSEQ_FA_FILE = "/data/Robinson-SB/PathSeq-data/pathseq_microbe.fa"
PATHSEQ_BACTERIA_FA_FILE = join("data", "{microbe}.fa")
# as of right now, there are 1954 bacteria files


localrules: download_bacteria_fna_refseq

include: "make-PathSeq-files.smk"

rule combine_fa_files:
    input:
        PATHSEQ_FA_FILE,
        BACTERIA_TO_ADD_FA_FILE
    output:
        PATHSEQ_BACTERIA_FA_FILE
    shell:
        "cat {input[0]} {input[1]} > {output}"



rule filter_bacteria_fna_refseq:
    conda:
        "../envs/biopython.yml"
    params:
        microbes_of_interest = config["PathSeq"]["microbes_of_interest"]
    input:
        expand(BACTERIA_FNA_FILE, fn=range(1,1955))
    output:
        BACTERIA_TO_ADD_FA_FILE
    script:
        "../src/filter_fna_files.py"


rule download_bacteria_fna_refseq:
    params:
        url = BACTERIA_REFSEQ_URL
    output:
        temp(BACTERIA_FNA_FILE)
    shell:
        "wget -O - {params.url}/bacteria.{wildcards.fn}.1.genomic.fna.gz | "
        "gunzip -c > {output}"
