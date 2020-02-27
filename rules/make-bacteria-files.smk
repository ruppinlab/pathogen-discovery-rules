from os.path import join

NCBI_REFSEQ_FTP = "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/"
BACTERIA_REFSEQ_URL = join(NCBI_REFSEQ_FTP, "bacteria")

BACTERIA_FNA_FILE = join("raw", "bacteria", "bacteria.{fn}.1.genomic.fna")
COMBINED_BACTERIA_FNA_FILE = join("raw", "bacteria", "combined.bacteria.genomic.fna")
# as of right now, there are 1954 bacteria files
# rule all:
#     input:
#         expand(BACTERIA_FNA_FILE, fn=range(1,5)),
#         COMBINED_BACTERIA_FNA_FILE

# TODO - look for contigs belonging to specific species

rule filter_bacteria_fna_refseq:
    conda:
        "../envs/biopython.yml"
    input:
        expand(BACTERIA_FNA_FILE, fn=range(1,1955))
    script:
        "../src/filter_fna_files.py"

# rule combine_bacteria_fna_files:
#     input:
#         expand(BACTERIA_FNA_FILE, fn=range(1,5))
#     output:
#         COMBINED_BACTERIA_FNA_FILE
#     shell:
#         """
#         for f in {input}
#         do
#             cat "$f" >> {output}
#         done
#         """

rule download_bacteria_fna_refseq:
    params:
        url = BACTERIA_REFSEQ_URL
    output:
        BACTERIA_FNA_FILE
    shell:
        "wget -O - {params.url}/bacteria.{wildcards.fn}.1.genomic.fna.gz | "
        "gunzip -c > {output}"
