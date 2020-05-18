from os.path import join, basename


# params
GATK_VERSION = "4.1.6.0"

rule PathSeqPipelineSpark_batch:
    input:
        bam_file = expand(config["PathSeq"]["bam_file"], zip, sample=samples["sample"], patient=samples["patient"]),
        host_bwa_image = config["PathSeq"]["host_img"],
        microbe_bwa_image = config["PathSeq"]["microbe_bwa_image"],
        microbe_fasta_file = config["PathSeq"]["microbe_fasta"],
        host_hss_file = config["PathSeq"]["host_bfi"],
        taxonomy_db = config["PathSeq"]["taxonomy_db"]
    params:
        host_bwa_image = basename(config["PathSeq"]["host_img"]),
        microbe_bwa_image = basename(config["PathSeq"]["microbe_bwa_image"]),
        microbe_fasta_file = basename(config["PathSeq"]["microbe_fasta"]),
        host_hss_file = basename(config["PathSeq"]["host_bfi"]),
        taxonomy_db = basename(config["PathSeq"]["taxonomy_db"])
    output:
        pathseq_bam = expand(join("output", "PathSeq", "{patient}-{sample}", "pathseq.bam"), zip, sample=samples["sample"], patient=samples["patient"]),
        pathseq_output = expand(join("output", "PathSeq", "{patient}-{sample}", "pathseq.txt"), zip, sample=samples["sample"], patient=samples["patient"]),
        filter_metrics = expand(join("output", "PathSeq", "{patient}-{sample}", "filter-metrics.txt"), zip, sample=samples["sample"], patient=samples["patient"])
    run:
        shell("mkdir /lscratch/$SLURM_JOBID/tmp")
        shell("cp /data/Robinson-SB/PathSeq-data/* /lscratch/$SLURM_JOBID")
        for bam_file, pathseq_bam, pathseq_output, filter_metrics in zip(input.bam_file, output.pathseq_bam, output.pathseq_output, output.filter_metrics):
            shell(
                "module load GATK/4.1.6.0 && "
                "gatk PathSeqPipelineSpark "
                "--input '{bam_file}' "
                "--filter-bwa-image /lscratch/$SLURM_JOBID/{params.host_bwa_image} "
                "--kmer-file /lscratch/$SLURM_JOBID/{params.host_hss_file} "
                "--microbe-fasta /lscratch/$SLURM_JOBID/{params.microbe_fasta_file} "
                "--microbe-bwa-image /lscratch/$SLURM_JOBID/{params.microbe_bwa_image} "
                "--taxonomy-file /lscratch/$SLURM_JOBID/{params.taxonomy_db} "
                "--output '{pathseq_bam}' "
                "--scores-output '{pathseq_output}' "
                "--filter-metrics '{filter_metrics}' "
                '--java-options "-Xmx64g -Xms64G -Djava.io.tmpdir=/lscratch/$SLURM_JOBID/tmp -XX:+UseG1GC -XX:ParallelGCThreads=8 -XX:ConcGCThreads=2" '
                '--spark-master local[8] '
                + config["params"]["PathSeq"]
            )

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
