from os.path import join, exists
import json


# Directories
ENV_DIR = join("..", "envs")
STAR_OUTPUT_DIR = join("output", "star", "{patient}-{sample}")
STAR_PASS1_OUTPUT_DIR = join(STAR_OUTPUT_DIR, "_STARpass1")
STAR_SE_PASS1_OUTPUT_DIR = join(STAR_OUTPUT_DIR, "_STARsepass1")
STAR_PASS2_OUTPUT_DIR = join(STAR_OUTPUT_DIR, "_STARpass2")
STAR_SE_PASS2_OUTPUT_DIR = join(STAR_OUTPUT_DIR, "_STARsepass2")
# Files
STAR_ENV_FILE = join(ENV_DIR, "star.yml")

# STAR output directories
STAR_GENOME_INDEX = join("output", "star-index")

# Intemediate Files
STAR_PASS1_SJ_FILE = join(STAR_PASS1_OUTPUT_DIR, "SJ.out.tab")
STAR_PASS1_SJ_FILTERED_FILE = join(STAR_PASS1_OUTPUT_DIR, "SJ.filtered.out.tab")
STAR_SE_PASS1_SJ_FILE = join(STAR_SE_PASS1_OUTPUT_DIR, "SJ.out.tab")
STAR_SE_PASS1_SJ_FILTERED_FILE = join(STAR_SE_PASS1_OUTPUT_DIR, "SJ.filtered.out.tab")
READLENGTH_HISTOGRAM = join(STAR_OUTPUT_DIR, "read-length-histogram.tsv")
SAMPLE_METADATA = join(STAR_OUTPUT_DIR, "sample-metadata.json")
# Output Files
STAR_PASS2_BAM_FILE = join(STAR_PASS2_OUTPUT_DIR, "Aligned.out.bam")
STAR_PASS2_READCOUNT_FILE = join(STAR_PASS2_OUTPUT_DIR, "ReadsPerGene.out.tab")
STAR_SE_PASS2_BAM_FILE = join(STAR_SE_PASS2_OUTPUT_DIR, "Aligned.out.bam")
STAR_SE_PASS2_READCOUNT_FILE = join(STAR_SE_PASS2_OUTPUT_DIR, "ReadsPerGene.out.tab")
HOST_READ_FILTERED_PE_BAM_FILE = join(STAR_SE_PASS2_OUTPUT_DIR, "host.reads.filtered.bam")
HOST_READ_FILTERED_SE_BAM_FILE = join(STAR_PASS2_OUTPUT_DIR, "host.reads.filtered.bam")
HOST_READ_FILTERED_BAM_FILE = join(STAR_OUTPUT_DIR, "host.reads.filtered.bam")

# set localrules
# localrules: compute_max_readlength, calculate_max_read_length, run_star_filter_sj_pass1, filter_aligned_reads, run_star_filter_sj_se_pass1

# functions
def get_sjdbOverhang(file):
    if not exists(file):
        return ""
    with open(file) as json_file:
        data = json.load(json_file)
        return data["sjdbOverhang"]

def get_fq(wildcards):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        dir = join("FASTQ", "raw")
    else:
        dir = join("FASTQ", "trimmed")
    return {
        'fq1': join(dir, "{wildcards.patient}-{wildcards.sample}_1.fastq.gz".format(wildcards=wildcards)),
        'fq2': join(dir, "{wildcards.patient}-{wildcards.sample}_2.fastq.gz".format(wildcards=wildcards))
        }


rule filter_aligned_reads:
    group:
        "STAR_2pass"
    input:
        STAR_PASS2_BAM_FILE,
        STAR_SE_PASS2_BAM_FILE,
    output:
        HOST_READ_FILTERED_PE_BAM_FILE,
        HOST_READ_FILTERED_SE_BAM_FILE,
        HOST_READ_FILTERED_BAM_FILE
    shell:
        "module load bamtools && "
        "bamtools filter -tag 'uT:<=2' -in {input[0]} -out {output[0]} && "
        "bamtools filter -tag 'uT:<=2' -in {input[1]} -out {output[1]} && "
        "bamtools merge -in {output[0]} -in {output[1]} -out {output[2]}"


rule create_star_index:
    conda:
        STAR_ENV_FILE
    input:
        config["ref"]["genome"]
    output:
        directory(STAR_GENOME_INDEX)
    threads:
        16
    shell:
        "mkdir '{output}' && STAR "
        "--runThreadN {threads} "
        "--runMode genomeGenerate "
        "--genomeDir '{output}' "
        "--genomeFastaFiles '{input}'"


# rules
rule compute_max_readlength:
    group:
        "STAR_2pass"
    conda:
        join(ENV_DIR, "pandas.yml")
    input:
        READLENGTH_HISTOGRAM,
    output:
        SAMPLE_METADATA
    script:
        "../src/extract_max_readlength.py"

rule calculate_max_read_length:
    group:
        "STAR_2pass"
    conda:
        join(ENV_DIR, "bbmap.yml")
    input:
        unpack(get_fq),
    output:
        READLENGTH_HISTOGRAM
    shell:
        "readlength.sh in='{input.fq1}' in2='{input.fq2}' bin=1 out='{output}'"


# reads for paired-end reads

rule run_star_pe_pass1:
    group:
        "STAR_2pass"
    conda:
        STAR_ENV_FILE
    input:
        unpack(get_fq),
        index = STAR_GENOME_INDEX,
        gtf = config["ref"]["annotation"],
        metadata = SAMPLE_METADATA
    output:
        STAR_PASS1_SJ_FILE
    params:
        odir = join(STAR_PASS1_OUTPUT_DIR, ""),
        sjdbOverhang = lambda wildcards, input: get_sjdbOverhang(input.metadata)
    threads:
        48
    shell:
        "STAR "
        "--runThreadN {threads} "
        "--readFilesIn '{input.fq1}' '{input.fq2}' "
        "--alignIntronMax 1000000 "
        "--alignIntronMin 20 "
        "--alignMatesGapMax 1000000 "
        "--alignSJDBoverhangMin 1 "
        "--alignSJoverhangMin 8 "
        "--alignSoftClipAtReferenceEnds Yes "
        "--genomeDir '{input.index}' "
        "--genomeLoad NoSharedMemory "
        "--limitSjdbInsertNsj 1200000 "
        "--outFileNamePrefix '{params.odir}' "
        "--outFilterIntronMotifs None "
        "--outFilterMismatchNmax 999 "
        "--outFilterMismatchNoverLmax 0.1 "
        "--outFilterMultimapNmax 20 "
        "--outSAMtype None "
        "--readFilesCommand zcat "
        "--sjdbGTFfile '{input.gtf}' "
        "--sjdbOverhang {params.sjdbOverhang}"
        + config["params"]["STAR"]

rule run_star_filter_sj_pass1:
    group:
        "STAR_2pass"
    input:
        STAR_PASS1_SJ_FILE
    output:
        STAR_PASS1_SJ_FILTERED_FILE
    shell:
        # $1~chromosomal and non-mitochondrial (regexp specific to GTF style!)
        # $5>0 canonical
        # $6==0 novel (since annotated get added from GTF)
        # $7>0 supported by at least one unique mapper
        "awk '$1~/chr([1-9][0-9]?|X|Y)/ && $5>0 && $6==0 && $7>0' "
        "{input} > {output}"

rule run_star_pe_pass2:
    group:
        "STAR_2pass"
    conda:
        STAR_ENV_FILE
    input:
        unpack(get_fq),
        index = STAR_GENOME_INDEX,
        gtf = config["ref"]["annotation"],
        sj = STAR_PASS1_SJ_FILTERED_FILE,
        metadata = SAMPLE_METADATA
    output:
        STAR_PASS2_BAM_FILE,
        STAR_PASS2_READCOUNT_FILE
    params:
        odir = join(STAR_PASS2_OUTPUT_DIR, ""),
        sjdbOverhang = lambda wildcards, input: get_sjdbOverhang(input.metadata)
    threads:
        48
    shell:
        "STAR "
        "--runThreadN {threads} "
        "--readFilesIn '{input.fq1}' '{input.fq2}' "
        "--alignIntronMax 1000000 "
        "--alignIntronMin 20 "
        "--alignMatesGapMax 1000000 "
        "--alignSJDBoverhangMin 1 "
        "--alignSJoverhangMin 8 "
        "--alignSoftClipAtReferenceEnds Yes "
        "--chimJunctionOverhangMin 15 "
        "--chimMainSegmentMultNmax 1 "
        "--chimOutType Junctions SeparateSAMold WithinBAM SoftClip "
        "--chimSegmentMin 15 "
        "--genomeDir '{input.index}' "
        "--genomeLoad NoSharedMemory "
        "--limitSjdbInsertNsj 1200000 "
        "--outFileNamePrefix '{params.odir}' "
        "--outFilterIntronMotifs None "
        "--outFilterMismatchNmax 999 "
        "--outFilterMismatchNoverLmax 0.1 "
        "--outFilterMultimapNmax 20 "
        "--outFilterType BySJout "
        "--outSAMattrRGline ID:{wildcards.patient}.{wildcards.sample} "
        "PL:illumina SM:{wildcards.patient}.{wildcards.sample} LB:RNA "
        "--outSAMattributes NH HI AS nM NM ch "
        "--outSAMstrandField intronMotif "
        "--outSAMtype BAM Unsorted "
        "--outSAMunmapped Within "
        "--quantMode GeneCounts "
        "--readFilesCommand zcat "
        "--sjdbGTFfile '{input.gtf}' "
        "--sjdbOverhang {params.sjdbOverhang} "
        "--sjdbFileChrStartEnd '{input.sj}' "
        + config["params"]["STAR"]

# reads for SE reads

rule run_star_se_pass1:
    group:
        "STAR_2pass"
    conda:
        STAR_ENV_FILE
    input:
        fq1 = join("FASTQ", "trimmed", "{patient}-{sample}_3.fastq.gz"),
        index = STAR_GENOME_INDEX,
        gtf = config["ref"]["annotation"],
        metadata = SAMPLE_METADATA
    output:
        STAR_SE_PASS1_SJ_FILE
    params:
        odir = join(STAR_SE_PASS1_OUTPUT_DIR, ""),
        sjdbOverhang = lambda wildcards, input: get_sjdbOverhang(input.metadata)
    threads:
        16
    shell:
        "STAR "
        "--runThreadN {threads} "
        "--readFilesIn '{input.fq1}' "
        "--alignIntronMax 1000000 "
        "--alignIntronMin 20 "
        "--alignMatesGapMax 1000000 "
        "--alignSJDBoverhangMin 1 "
        "--alignSJoverhangMin 8 "
        "--alignSoftClipAtReferenceEnds Yes "
        "--genomeDir '{input.index}' "
        "--genomeLoad NoSharedMemory "
        "--limitSjdbInsertNsj 1200000 "
        "--outFileNamePrefix '{params.odir}' "
        "--outFilterIntronMotifs None "
        "--outFilterMismatchNmax 999 "
        "--outFilterMismatchNoverLmax 0.1 "
        "--outFilterMultimapNmax 20 "
        "--outSAMtype None "
        "--readFilesCommand zcat "
        "--sjdbGTFfile '{input.gtf}' "
        "--sjdbOverhang {params.sjdbOverhang}"
        + config["params"]["STAR"]

rule run_star_filter_sj_se_pass1:
    group:
        "STAR_2pass"
    input:
        STAR_SE_PASS1_SJ_FILE
    output:
        STAR_SE_PASS1_SJ_FILTERED_FILE
    shell:
        # $1~chromosomal and non-mitochondrial (regexp specific to GTF style!)
        # $5>0 canonical
        # $6==0 novel (since annotated get added from GTF)
        # $7>0 supported by at least one unique mapper
        "awk '$1~/chr([1-9][0-9]?|X|Y)/ && $5>0 && $6==0 && $7>0' "
        "{input} > {output}"

rule run_star_se_pass2:
    group:
        "STAR_2pass"
    conda:
        STAR_ENV_FILE
    input:
        fq1 = join("FASTQ", "trimmed", "{patient}-{sample}_3.fastq.gz"),
        index = STAR_GENOME_INDEX,
        gtf = config["ref"]["annotation"],
        sj = STAR_SE_PASS1_SJ_FILTERED_FILE,
        metadata = SAMPLE_METADATA
    output:
        STAR_SE_PASS2_BAM_FILE,
        STAR_SE_PASS2_READCOUNT_FILE
    params:
        odir = join(STAR_SE_PASS2_OUTPUT_DIR, ""),
        sjdbOverhang = lambda wildcards, input: get_sjdbOverhang(input.metadata)
    threads:
        48
    shell:
        "STAR "
        "--runThreadN {threads} "
        "--readFilesIn '{input.fq1}' "
        "--alignIntronMax 1000000 "
        "--alignIntronMin 20 "
        "--alignMatesGapMax 1000000 "
        "--alignSJDBoverhangMin 1 "
        "--alignSJoverhangMin 8 "
        "--alignSoftClipAtReferenceEnds Yes "
        "--chimJunctionOverhangMin 15 "
        "--chimMainSegmentMultNmax 1 "
        "--chimOutType Junctions SeparateSAMold WithinBAM SoftClip "
        "--chimSegmentMin 15 "
        "--genomeDir '{input.index}' "
        "--genomeLoad NoSharedMemory "
        "--limitSjdbInsertNsj 1200000 "
        "--outFileNamePrefix '{params.odir}' "
        "--outFilterIntronMotifs None "
        "--outFilterMismatchNmax 999 "
        "--outFilterMismatchNoverLmax 0.1 "
        "--outFilterMultimapNmax 20 "
        "--outFilterType BySJout "
        "--outSAMattrRGline ID:{wildcards.patient}.{wildcards.sample} "
        "PL:illumina SM:{wildcards.patient}.{wildcards.sample} LB:RNA "
        "--outSAMattributes NH HI AS nM NM ch "
        "--outSAMstrandField intronMotif "
        "--outSAMtype BAM Unsorted "
        "--outSAMunmapped Within "
        "--quantMode GeneCounts "
        "--readFilesCommand zcat "
        "--sjdbGTFfile '{input.gtf}' "
        "--sjdbOverhang {params.sjdbOverhang} "
        "--sjdbFileChrStartEnd '{input.sj}' "
        + config["params"]["STAR"]
