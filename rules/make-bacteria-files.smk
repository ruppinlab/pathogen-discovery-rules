from os.path import join

NCBI_REFSEQ_FTP = "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/"
BACTERIA_REFSEQ_URL = join(NCBI_REFSEQ_FTP, "bacteria")

BACTERIA_FNA_FILE = join("raw", "bacteria", "bacteria.{fn}.1.genomic.fna")
BACTERIA_OF_INTERST_FNA_FILE = join("output", "bacteria", "{microbe}.genomic.fna")
# as of right now, there are 1954 bacteria files
# rule all:
#     input:
#         expand(BACTERIA_FNA_FILE, fn=range(1,5)),
#         COMBINED_BACTERIA_FNA_FILE

# TODO - look for contigs belonging to specific species

include: "make-PathSeq-files.smk"


rule all:
    input:
        expand(BACTERIA_OF_INTERST_FNA_FILE, microbe="salmonella-Aulicino2018")


rule filter_bacteria_fna_refseq:
    conda:
        "../envs/biopython.yml"
    params:
        microbes_of_interest = ["Salmonella enterica subsp. enterica serovar Typhimurium str. SL1344", "Salmonella enterica subsp. enterica serovar Typhimurium str. LT2"]
    input:
        expand(BACTERIA_FNA_FILE, fn=range(1,1955))
    output:
        BACTERIA_OF_INTERST_FNA_FILE
    script:
        "../src/filter_fna_files.py"


rule download_bacteria_fna_refseq:
    wildcard_constraints:
        microbe = "bacteria"
    params:
        url = BACTERIA_REFSEQ_URL
    output:
        BACTERIA_FNA_FILE
    shell:
        "wget -O - {params.url}/bacteria.{wildcards.fn}.1.genomic.fna.gz | "
        "gunzip -c > {output}"
