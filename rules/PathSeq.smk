from os.path import join

HOST_HSS_FILE = join("output", "PathSeq", "host.hss")
HOST_BWA_IMAGE_INDEX = join("output", "PathSeq", "host.fasta.img")

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

rule run_PathSeq:
    input:
        bam_file = config["PathSeq"]["bam_file"],
        host_bwa_image = HOST_BWA_IMAGE_INDEX,
        microbe_bwa_image = config["PathSeq"]["microbe_bwa_image"],
        microbe_fasta_file = config["PathSeq"]["microbe_fasta"],
        host_hss_file = HOST_HSS_FILE,
        taxonomy_db = config["PathSeq"]["taxonomy_db"]
    output:
        pathseq_bam = join("output", "PathSeq", "{patient}-{sample}", "pathseq.bam"),
        pathseq_output = join("output", "PathSeq", "{patient}-{sample}", "pathseq.txt")
    shell:
        "module load GATK/4.1.3.0 && "
        "gatk PathSeqPipelineSpark "
        '--java-options "-Xmx115g -Djava.io.tmpdir=/lscratch/$SLURM_JOBID" ' # double quotes escapes the wildcard
        "--spark-master local[*] "  # use all available cores
        "--is-host-aligned true "
        "--input '{input.bam_file}' "
        "--filter-bwa-image '{input.host_bwa_image}' "
        "--kmer-file '{input.host_hss_file}' "
        "--microbe-fasta '{input.microbe_fasta_file}' "
        "--microbe-bwa-image '{input.microbe_bwa_image}' "
        "--taxonomy-file '{input.taxonomy_db}' "
        "--output '{output.pathseq_bam}' "
        "--scores-output '{output.pathseq_output}' "
        + config["params"]["PathSeq"]


# Rules for building host files
rule build_host_kmer_file:
    input:
        config["ref"]["genome"]
    output:
        HOST_HSS_FILE
    shell:
        "module load GATK/4.1.3.0 && "
        "gatk PathSeqBuildKmers "
        "--java-options '-Xmx80g' "
        "--reference '{input}' "
        "-O '{output}'"

rule build_host_BWA_image:
    input:
        config["ref"]["genome"]
    output:
        HOST_BWA_IMAGE_INDEX
    shell:
        "module load GATK/4.1.3.0 && "
        "gatk BwaMemIndexImageCreator -I {input} -O {output}"
