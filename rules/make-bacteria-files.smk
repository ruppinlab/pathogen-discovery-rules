from os.path import join

NCBI_REFSEQ_FTP = "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/"
BACTERIA_REFSEQ_URL = join(NCBI_REFSEQ_FTP, "bacteria")

BACTERIA_FNA_FILE = join("raw", "bacteria", "bacteria.{fn}.1.genomic.fna")
BACTERIA_OF_INTERST_FA_FILE = join("output", "bacteria", "{microbe}.fa")
# as of right now, there are 1954 bacteria files
# rule all:
#     input:
#         expand(BACTERIA_FNA_FILE, fn=range(1,5)),
#         COMBINED_BACTERIA_FNA_FILE

# TODO - look for contigs belonging to specific species

include: "make-PathSeq-files.smk"


rule filter_bacteria_fna_refseq:
    conda:
        "../envs/biopython.yml"
    params:
        microbes_of_interest = config["PathSeq"]["microbes_of_interest"]
    input:
        expand(BACTERIA_FNA_FILE, fn=range(1,1955))
    output:
        BACTERIA_OF_INTERST_FA_FILE
    script:
        "../src/filter_fna_files.py"


rule download_bacteria_fna_refseq:
    params:
        url = BACTERIA_REFSEQ_URL
    output:
        BACTERIA_FNA_FILE
    shell:
        "wget -O - {params.url}/bacteria.{wildcards.fn}.1.genomic.fna.gz | "
        "gunzip -c > {output}"
