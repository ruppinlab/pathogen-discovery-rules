from os.path import join, exists


# Directories
STAR_OUTPUT_DIR = join("output", "star", "{patient}-{sample}-{plate}")
STAR_PE_OUTPUT_DIR = join(STAR_OUTPUT_DIR, "_STARpe")
STAR_SE_OUTPUT_DIR = join(STAR_OUTPUT_DIR, "_STARse")

# Input files
STAR_PE_BAM_FILE = join(STAR_PE_OUTPUT_DIR, "Aligned.out.bam")
STAR_SE_BAM_FILE = join(STAR_SE_OUTPUT_DIR, "Aligned.out.bam")

# Intermediate Files
STAR_PE_UNALIGNED_BAM_FILE = join(STAR_PE_OUTPUT_DIR, "unaligned.bam")
STAR_SE_UNALIGNED_BAM_FILE = join(STAR_SE_OUTPUT_DIR, "unaligned.bam")

# Output Files
STAR_UNALIGNED_BAM_FILE = join(STAR_OUTPUT_DIR, "unaligned.bam")

# set localrules
# localrules: compute_max_readlength, calculate_max_read_length, run_star_filter_sj_pass1, filter_aligned_reads, run_star_filter_sj_se_pass1


rule combine_se_pe_unaligned_reads:
    group:
        "STAR_2pass"
    input:
        STAR_PE_UNALIGNED_BAM_FILE,
        STAR_SE_UNALIGNED_BAM_FILE,
    output:
        STAR_UNALIGNED_BAM_FILE
    benchmark:
        "benchmarks/{patient}-{sample}-{plate}.combine_se_pe_aligned_reads.benchmark.txt"
    shell:
        "module load bamtools && "
        "bamtools merge -in {input[0]} -in {input[1]} -out {output[0]}"

rule filter_aligned_pe_reads:
    group:
        "STAR_2pass"
    input:
        STAR_PE_BAM_FILE,
    output:
        STAR_PE_UNALIGNED_BAM_FILE
    benchmark:
        "benchmarks/{patient}-{sample}-{plate}.filter_pe_aligned_reads.benchmark.txt"
    shell:
        "module load bamtools && "
        "bamtools filter -tag 'uT:<=2' -in {input[0]} -out {output[0]}"

rule filter_aligned_se_reads:
    group:
        "STAR_2pass"
    input:
        STAR_SE_BAM_FILE,
    output:
        STAR_SE_UNALIGNED_BAM_FILE
    benchmark:
        "benchmarks/{patient}-{sample}-{plate}.filter_se_aligned_reads.benchmark.txt"
    shell:
        "module load bamtools && "
        "bamtools filter -tag 'uT:<=2' -in {input[0]} -out {output[0]}"
