from os.path import join, basename

#ruleorder: PathSeqFilterSpark_paired_only > PathSeqFilterSpark_unpaired_only > PathSeqFilterSpark

#HOST_HSS_FILE = join("output", "PathSeq", "host.hss")
#HOST_BWA_IMAGE_INDEX = join("output", "PathSeq", "host.fasta.img")

# PAIRED_FILTERED_BAM = join("output", "PathSeq", "{patient}-{sample}", "filtered-paired.bam")
# UNPAIRED_FILTERED_BAM = join("output", "PathSeq", "{patient}-{sample}", "filtered-unpaired.bam")
# PATHSEQ_FILTER_FILE = join("output", "PathSeq", "{patient}-{sample}", "filter-metrics.txt")
# PAIRED_ALIGNED_BAM = join("output", "PathSeq", "{patient}-{sample}", "aligned-paired.bam")
# UNPAIRED_ALIGNED_BAM = join("output", "PathSeq", "{patient}-{sample}", "aligned-unpaired.bam")

# samples_with_only_paired_reads = [
#     "SRR7667576", "SRR7667864", "SRR7667841", "SRR7667837", "SRR7667795", "SRR7667801",
#     "SRR7667613", "SRR7667808", "SRR7667796", "SRR7667630", "SRR7667629"
#     ]
#
# samples_with_only_unpaired_reads = [
#     "SRR7667753", "SRR7667733", "SRR7667734", "SRR7667836", "SRR7667717", "SRR7667771",
#     "SRR7667730", "SRR7667727", "SRR7667578", "SRR7667776", "SRR7667852", "SRR7667846",
#     "SRR7667853", "SRR7667847", "SRR7667573", "SRR7667569", "SRR7667855", "SRR7667804",
#     "SRR7667854", "SRR7667844", "SRR7667833", "SRR7667802", "SRR7667821", "SRR7667786",
#     "SRR7667789", "SRR7667842", "SRR7667785", "SRR7667835", "SRR7667621", "SRR7667827",
#     "SRR7667791", "SRR7667572", "SRR7667793", "SRR7667784", "SRR7667575", "SRR7667815",
#     "SRR7667825", "SRR7667857", "SRR7667863", "SRR7667625", "SRR7667623", "SRR7667797",
#     "SRR7667829", "SRR7667775", "SRR7667794", "SRR7667822", "SRR7667777", "SRR7667571",
#     "SRR7667609", "SRR7667772", "SRR7667620", "SRR7667798", "SRR7667807", "SRR7667858",
#     "SRR7667622", "SRR7667626", "SRR7667781", "SRR7667809", "SRR7667579", "SRR7667805",
#     "SRR7667774", "SRR7667830", "SRR7667867", "SRR7667574", "SRR7667614", "SRR7667616",
#     "SRR7667861", "SRR7667610", "SRR7667818", "SRR7667619", "SRR7667627", "SRR7667783",
#     "SRR7667618", "SRR7667866", "SRR7667813", "SRR7667607", "SRR7667611", "SRR7667624",
#     "SRR7667859", "SRR7667782", "SRR7667612", "SRR7667628", "SRR7667617"
#     ]
# params
GATK_VERSION = "4.1.6.0"
# rule get_num_properly_paired_microbial_reads:
#     input:
#         PATHSEQ_BAM_FILE
#     output:
#         PROPERLY_PAIRED_MICROBE_FILE
#     shell:
#         "module load samtools && "
#         "samtools view -f 0x2 {input} | wc -l > {output}"
#
# rule get_num_properly_paired_hg38_reads:
#     input:
#         STAR_PASS2_BAM_FILE
#     output:
#         PROPERLY_PAIRED_HG38_FILE
#     shell:
#         "module load samtools && "
#         "samtools view -f 0x2 {input} | wc -l > {output}"

# functions
# def get_ref_genome(wildcards):
#     try:
#         return config["PathSeq"]["genome"]
#     except:
#         return config["ref"]["genome"]

# rules for running PathSeqSpark

# rule PathSeqFilterSpark_paired_only:
#     wildcard_constraints:
#         sample = "|".join(samples_with_only_paired_reads)
#     input:
#         bam_file = config["PathSeq"]["bam_file"],
#         host_bwa_image = HOST_BWA_IMAGE_INDEX,
#         host_hss_file = HOST_HSS_FILE
#     output:
#         paired_output = PAIRED_FILTERED_BAM,
#         filter_metrics = PATHSEQ_FILTER_FILE
#     shell:
#         "module load GATK/{GATK_VERSION} && "
#         "gatk PathSeqFilterSpark "
#         "--is-host-aligned true "
#         "--input '{input.bam_file}' "
#         "--filter-bwa-image '{input.host_bwa_image}' "
#         "--kmer-file '{input.host_hss_file}' "
#         "--paired-output '{output.paired_output}' "
#         "--filter-metrics '{output.filter_metrics}' "
#         + config["params"]["PathSeq"]["filter"]
#
# rule PathSeqFilterSpark_unpaired_only:
#     wildcard_constraints:
#         sample = "|".join(samples_with_only_unpaired_reads)
#     input:
#         bam_file = config["PathSeq"]["bam_file"],
#         host_bwa_image = HOST_BWA_IMAGE_INDEX,
#         host_hss_file = HOST_HSS_FILE
#     output:
#         unpaired_output = UNPAIRED_FILTERED_BAM,
#         filter_metrics = PATHSEQ_FILTER_FILE
#     shell:
#         "module load GATK/{GATK_VERSION} && "
#         "gatk PathSeqFilterSpark "
#         "--is-host-aligned true "
#         "--input '{input.bam_file}' "
#         "--filter-bwa-image '{input.host_bwa_image}' "
#         "--kmer-file '{input.host_hss_file}' "
#         "--unpaired-output '{output.unpaired_output}' "
#         "--filter-metrics '{output.filter_metrics}' "
#         + config["params"]["PathSeq"]["filter"]


rule copy_PathSeqFilter_files_to_lscratch:
    input:
        host_bwa_image = config["PathSeq"]["host_img"],
        host_hss_file = config["PathSeq"]["host_bfi"]
    output:
        temp(touch("PathSeqFilter-hack-{batch}.txt"))
    group:
        "PathSeqFilter"
    shell:
        "mkdir /lscratch/$SLURM_JOBID/tmp && "
        "cp {input.host_bwa_image} {input.host_hss_file} /lscratch/$SLURM_JOBID"

rule PathSeqFilterSpark:
    input:
        bam_file = config["PathSeq"]["bam_file"],
        hack = "PathSeqFilter-hack-{batch}.txt"
    params:
        host_bwa_image = basename(config["PathSeq"]["host_img"]),
        host_hss_file = basename(config["PathSeq"]["host_bfi"])
    output:
        paired_output = temp(join("output", "PathSeq", "{patient}-{sample}", "{batch}-filtered-paired.bam")),
        unpaired_output = temp(join("output", "PathSeq", "{patient}-{sample}", "{batch}-filtered-unpaired.bam")),
        filter_metrics = temp(join("output", "PathSeq", "{patient}-{sample}", "{batch}-filter-metrics.txt")),
    group:
        "PathSeqFilter"
    threads:
        8
    shell:
        "module load GATK/{GATK_VERSION} && "
        "gatk PathSeqFilterSpark "
        "--input '{input.bam_file}' "
        "--filter-bwa-image /lscratch/$SLURM_JOBID/{params.host_bwa_image} "
        "--kmer-file /lscratch/$SLURM_JOBID/{params.host_hss_file} "
        "--paired-output '{output.paired_output}' "
        "--unpaired-output '{output.unpaired_output}' "
        "--filter-metrics '{output.filter_metrics}' "
        + config["params"]["PathSeq"]["filter"]

rule copy_PathSeqBwa_files_to_lscratch:
    input:
        microbe_bwa_image = config["PathSeq"]["microbe_bwa_image"],
        microbe_fasta_file = config["PathSeq"]["microbe_fasta"]
    params:
        microbe_path = config["PathSeq"]["microbe_fasta"].replace("fa", "")
    output:
        temp(touch("PathSeqBwa-hack-{batch}.txt"))
    group:
        "PathSeqBwa"
    shell:
        "mkdir /lscratch/$SLURM_JOBID/tmp && "
        "cp {params.microbe_path}* /lscratch/$SLURM_JOBID"


rule PathSeqBwaSpark:
    input:
        paired_input = PAIRED_FILTERED_BAM,
        unpaired_input = UNPAIRED_FILTERED_BAM,
        h = "PathSeqBwa-hack-{batch}.txt"
    params:
        microbe_bwa_image = basename(config["PathSeq"]["microbe_bwa_image"]),
        microbe_fasta_file = basename(config["PathSeq"]["microbe_fasta"])
    output:
        paired_output = join("output", "PathSeq", "{patient}-{sample}", "{batch}-aligned-paired.bam"),
        unpaired_output = join("output", "PathSeq", "{patient}-{sample}", "{batch}-aligned-unpaired.bam")
    threads:
        8
    group:
        "PathSeqBwa"
    shell:
        "module load GATK/{GATK_VERSION} && "
        "gatk PathSeqBwaSpark "
        "--paired-input '{input.paired_input}' "
        "--unpaired-input '{input.unpaired_input}' "
        "--microbe-fasta /lscratch/$SLURM_JOBID/{params.microbe_fasta_file} "
        "--microbe-bwa-image /lscratch/$SLURM_JOBID/{params.microbe_bwa_image} "
        "--paired-output '{output.paired_output}' "
        "--unpaired-output '{output.unpaired_output}' "
        + config["params"]["PathSeq"]["BWA"]

# rule PathSeqScoreSpark:
#     input:
#         paired_input = PAIRED_ALIGNED_BAM,
#         unpaired_input = UNPAIRED_ALIGNED_BAM,
#         taxonomy_db = config["PathSeq"]["taxonomy_db"]
#     output:
#         score_file =,
#         output_bam = ,# annotated with the NCBI taxonomy IDs of mapped organisms
#



# rule PathSeqPipelineSpark:
#     input:
#         bam_file = config["PathSeq"]["bam_file"],
#         host_bwa_image = HOST_BWA_IMAGE_INDEX,
#         microbe_bwa_image = config["PathSeq"]["microbe_bwa_image"],
#         microbe_fasta_file = config["PathSeq"]["microbe_fasta"],
#         host_hss_file = HOST_HSS_FILE,
#         taxonomy_db = config["PathSeq"]["taxonomy_db"]
#     output:
#         pathseq_bam = join("output", "PathSeq", "{patient}-{sample}", "pathseq.bam"),
#         pathseq_output = join("output", "PathSeq", "{patient}-{sample}", "pathseq.txt")
#     shell:
#         "module load GATK/{GATK_VERSION} && "
#         "gatk PathSeqPipelineSpark "
#         "--is-host-aligned true "
#         "--input '{input.bam_file}' "
#         "--filter-bwa-image '{input.host_bwa_image}' "
#         "--kmer-file '{input.host_hss_file}' "
#         "--microbe-fasta '{input.microbe_fasta_file}' "
#         "--microbe-bwa-image '{input.microbe_bwa_image}' "
#         "--taxonomy-file '{input.taxonomy_db}' "
#         "--output '{output.pathseq_bam}' "
#         "--scores-output '{output.pathseq_output}' "
#         + config["params"]["PathSeq"]["pipeline"]


# Rules for building host files

rule build_host_kmer_file:
    input:
        config["ref"]["genome"]
    output:
        config["PathSeq"]["host_bfi"]
    shell:
        "module load GATK/{GATK_VERSION} && "
        "gatk PathSeqBuildKmers "
        "--java-options '-Xmx80g' "
        "--reference '{input}' "
        "-O '{output}'"

rule build_host_BWA_image:
    input:
        config["ref"]["genome"]
    output:
        config["PathSeq"]["host_img"]
    shell:
        "module load GATK/{GATK_VERSION} && "
        "gatk BwaMemIndexImageCreator -I {input} -O {output}"
