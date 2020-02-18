from os.path import join

# URLs
# currently using release 97 - usure how to include this as part of the URL
NCBI_REFSEQ_FTP = "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/"
VIRUS_REFSEQ_URL = join(NCBI_REFSEQ_FTP, "viral")
REFSEQ_CATALOG_URL = join(NCBI_REFSEQ_FTP, "release-catalog/RefSeq-release98.catalog.gz")
NCBI_TAX_DUMP_URL = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

# Directories
RAW_DIR = "raw"

# Files to download from NCBI FTP
# Virus data
VIRUS_FNA1 = join(RAW_DIR, "viral.1.1.genomic.fna")
VIRUS_FNA2 = join(RAW_DIR, "viral.2.1.genomic.fna")
VIRUS_FNA3 = join(RAW_DIR, "viral.3.1.genomic.fna")

# other downloaded files
REFSEQ_CATALOG = join(RAW_DIR, "RefSeq-release98.catalog.gz")
NCBI_TAX_DUMP = join(RAW_DIR, "ncbi_taxdump.tar.gz")

# output files
VIRAL_TAXONOMY = join(DATA_DIR, "viral_taxonomy.db")
VIRAL_FASTA_FILE = join(DATA_DIR, "viral.fa")
VIRAL_FASTA_IDX_FILE = join(DATA_DIR, "viral.fa.fai")
VIRAL_FASTA_DICT_FILE = join(DATA_DIR, "viral.dict")
VIRAL_BWA_IMAGE_INDEX = join(DATA_DIR, "viral.fa.img")

# rules for generating PathSeq data
rule download_NCBI_taxonomy_dump:
    params:
        NCBI_TAX_DUMP_URL
    output:
        NCBI_TAX_DUMP
    shell:
        "wget {params} -O {output}"

rule download_RefSeq_accession_catalog:
    params:
        REFSEQ_CATALOG_URL
    output:
        REFSEQ_CATALOG
    shell:
        "wget {params} -O {output}"

# rules for building viral files
rule build_viral_taxonomy_file:
    input:
        fa = VIRAL_FASTA_FILE,
        catalog = REFSEQ_CATALOG,
        taxdump = NCBI_TAX_DUMP,
        fai = VIRAL_FASTA_IDX_FILE,
        dict = VIRAL_FASTA_DICT_FILE
    output:
        VIRAL_TAXONOMY
    shell:
        "module load GATK/4.1.3.0 && "
        "gatk PathSeqBuildReferenceTaxonomy "
        "-R '{input.fa}' "
        "--refseq-catalog '{input.catalog}' "
        "--tax-dump '{input.taxdump}' "
        "-O '{output}'"

rule create_fasta_dict:
    input:
        VIRAL_FASTA_FILE
    output:
        VIRAL_FASTA_DICT_FILE
    shell:
        "module load picard && java -jar $PICARDJARPATH/picard.jar "
        "CreateSequenceDictionary R= '{input}' O= '{output}'"

rule create_viral_fasta_index_file:
    input:
        VIRAL_FASTA_FILE
    output:
        VIRAL_FASTA_IDX_FILE
    shell:
        "module load samtools && samtools faidx {input}"

rule combine_virus_fna:
    input:
        fna1 = VIRUS_FNA1,
        fna2 = VIRUS_FNA2,
        fna3 = VIRUS_FNA3
    output:
        VIRAL_FASTA_FILE
    shell:
        "cat {input.fna1} {input.fna2} {input.fna3} > {output}"

rule download_virus_fna_refseq:
    params:
        url = VIRUS_REFSEQ_URL
    output:
        fna1 = VIRUS_FNA1,
        fna2 = VIRUS_FNA2,
        fna3 = VIRUS_FNA3
    shell:
        "wget -O - {params.url}/viral.1.1.genomic.fna.gz | "
        "gunzip -c > {output.fna1} && "
        "wget -O - {params.url}/viral.2.1.genomic.fna.gz | "
        "gunzip -c > {output.fna2} && "
        "wget -O - {params.url}/viral.3.1.genomic.fna.gz | "
        "gunzip -c > {output.fna3}"

rule build_viral_BWA_image:
    input:
        VIRAL_FASTA_FILE
    output:
        VIRAL_BWA_IMAGE_INDEX
    shell:
        "module load GATK/4.1.3.0 && "
        "gatk BwaMemIndexImageCreator -I {input} -O {output}"
