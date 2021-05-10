

CAMMIQ_INDICES_DIR = "/data/CDSL_Gurobi_users/CAMMiQ/new_indices"
GENOME_MAP = join(CAMMIQ_INDICES_DIR, "genome_map_{tax_level}.out")
CAMMIQ_INDEX_BIN1 = join(CAMMIQ_INDICES_DIR, "index_{tax_level}_u_38.bin1")
CAMMIQ_INDEX_BIN2 = join(CAMMIQ_INDICES_DIR, "index_{tax_level}_d_38.bin2")

CAMMIQ_COUNT_FILE = join("output", "CAMMiQ", "read_cnts_{tax_level}.txt")

PATHSEQ_CELL_BAM = join("output", "PathSeq", "{patient}-{sample}-{cell}", "pathseq_with_tags.bam")

CAMMIQ_FASTQ_DIR = join("output", "CAMMiQ", "fastq_dir")
PathSeq_FQ1 = join(CAMMIQ_FASTQ_DIR, "{patient}-{sample}-{cell}_1.fq")

rule convert_BAMs_to_fastq_dir:
    group:
        "convert_BAMs_to_fastq_dir"
    input:
        bam_file=PATHSEQ_CELL_BAM
    output:
        fq1=PathSeq_FQ1,
    shell:
        "module load bedtools && "
        "bamToFastq -i {input.bam_file} -fq {output.fq1}"

rule run_CAMMiQ_species_long_reads:
    input:
        fq1=expand(PathSeq_FQ1, zip, patient=cells["patient"], sample=cells["sample"], plate=cells["plate"], cell=cells["cell"]),
    params:
        GENOME_MAP,
        CAMMIQ_INDEX_BIN1,
        CAMMIQ_INDEX_BIN2,
        CAMMIQ_FASTQ_DIR + "/"
    output:
        CAMMIQ_COUNT_FILE
    shell:
        "set -e && "
        "export GUROBI_HOME=/data/CDSL_Gurobi_users/gurobi910/linux64 && "
        "export GRB_LICENSE_FILE=/data/CDSL_Gurobi_users/gurobi910/gurobi.lic && "
        "module use --prepend /data/CDSL_Gurobi_users/modules && "
        "module load gurobi/9.1.0 && "
        "module load gcc/7.4.0  && "
        "/data/CDSL_Gurobi_users/CAMMiQ/src_scrna/run_gurobi910 --query -i doubly_unique "
        "-h 26 26 -f {params[0]} {params[1]} {params[2]} "
        "-o {output} -d {params[3]} -a 10 5 -minL 51"
