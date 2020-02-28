from os.path import join

# URLs
# currently using release 97 - usure how to include this as part of the URL
NCBI_REFSEQ_FTP = "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/"
VIRUS_REFSEQ_URL = join(NCBI_REFSEQ_FTP, "viral")
REFSEQ_CATALOG_URL = join(NCBI_REFSEQ_FTP, "release-catalog/RefSeq-release98.catalog.gz")
NCBI_TAX_DUMP_URL = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

# Directories
RAW_DIR = "raw"
DATA_DIR = "data"

# Files to download from NCBI FTP
# Virus data
VIRUS_FNA_FILE = join(RAW_DIR, "viral.{fn}.1.genomic.fna")

# output files
VIRAL_FASTA_FILE = join(DATA_DIR, "viral.fa")



rule combine_virus_fna_files:
    input:
        expand(VIRUS_FNA_FILE, fn=range(1,4))
    output:
        VIRAL_FASTA_FILE
    shell:
        """
        for f in {input}
        do
            cat "$f" >> {output}
        done
        """

rule download_fna_refseq:
    params:
        url = NCBI_REFSEQ_FTP
    output:
        VIRAL_FNA_FILE,
    shell:
        "wget -O - {params.url}/viral.{wildcards.fn}.1.genomic.fna.gz | "
        "gunzip -c > {output}"
