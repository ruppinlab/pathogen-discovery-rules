from os.path import join

# URLs
# currently using release 97 - usure how to include this as part of the URL
REFSEQ_CATALOG_URL = "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/archive/RefSeq-{release}.catalog.gz"
NCBI_TAX_DUMP_URL = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_{date}.zip"

# Directories
RAW_DIR = "raw"
DATA_DIR = "data"

# other downloaded files
REFSEQ_CATALOG = join(RAW_DIR, "RefSeq-{release}.catalog.gz")
NCBI_TAX_DUMP = join(RAW_DIR, "ncbi_taxdump_{date}.zip")

# output files
TAXONOMY = join(DATA_DIR, "{microbe}_{release}_{date}_taxonomy.db")
FASTA_FILE = join(DATA_DIR, "{microbe}.fa")
FASTA_IDX_FILE = join(DATA_DIR, "{microbe}.fa.fai")
FASTA_DICT_FILE = join(DATA_DIR, "{microbe}.dict")
BWA_IMAGE_INDEX = join(DATA_DIR, "{microbe}.fa.img")

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
rule build_taxonomy_file:
    input:
        fa = FASTA_FILE,
        catalog = REFSEQ_CATALOG,
        taxdump = NCBI_TAX_DUMP,
        fai = FASTA_IDX_FILE,
        dict = FASTA_DICT_FILE
    output:
        TAXONOMY
    shell:
        "module load GATK/4.1.8.1 && "
        "gatk PathSeqBuildReferenceTaxonomy "
        "-R '{input.fa}' "
        "--refseq-catalog '{input.catalog}' "
        "--tax-dump '{input.taxdump}' "
        "-O '{output}'"

rule create_fasta_dict:
    input:
        FASTA_FILE
    output:
        FASTA_DICT_FILE
    shell:
        "module load picard && java -jar $PICARDJARPATH/picard.jar "
        "CreateSequenceDictionary R= '{input}' O= '{output}'"

rule create_fasta_index_file:
    input:
        FASTA_FILE
    output:
        FASTA_IDX_FILE
    shell:
        "module load samtools && samtools faidx {input}"

rule build_BWA_image:
    input:
        FASTA_FILE
    output:
        BWA_IMAGE_INDEX
    shell:
        "module load GATK/4.1.8.1 && "
        "gatk BwaMemIndexImageCreator -I {input} -O {output}"
