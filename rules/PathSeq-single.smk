from os.path import join, basename


# params
GATK_VERSION = "4.1.8.1"

rule PathSeqPipelineSpark:
    input:
        bam_file = config["PathSeq"]["bam_file"],
        microbe_bwa_image = config["PathSeq"]["microbe_bwa_image"],
        microbe_dict_file = config["PathSeq"]["microbe_dict"],
        taxonomy_db = config["PathSeq"]["taxonomy_db"]
    params:
        microbe_bwa_image = basename(config["PathSeq"]["microbe_bwa_image"]),
        microbe_dict_file = basename(config["PathSeq"]["microbe_dict"]),
        taxonomy_db = basename(config["PathSeq"]["taxonomy_db"])
    output:
        pathseq_bam = join("output", "PathSeq", "{patient}-{sample}", "pathseq.bam"),
        pathseq_output = join("output", "PathSeq", "{patient}-{sample}", "pathseq.txt"),
        filter_metrics = join("output", "PathSeq", "{patient}-{sample}", "filter-metrics.txt"),
        score_metrics = join("output", "PathSeq", "{patient}-{sample}", "score-metrics.txt"),
    run:
        shell("mkdir /lscratch/$SLURM_JOBID/tmp")
        shell("cp {input.microbe_bwa_image} /lscratch/$SLURM_JOBID/")
        shell("cp {input.microbe_dict_file} /lscratch/$SLURM_JOBID/")
        shell("cp {input.taxonomy_db} /lscratch/$SLURM_JOBID/")
        shell(
            "module load GATK/4.1.8.1 && "
            "gatk PathSeqPipelineSpark "
            "--skip-quality-filters true "
            "--filter-duplicates false "
            "--min-score-identity .7 "
            "--input '{input.bam_file}' "
            "--microbe-bwa-image /lscratch/$SLURM_JOBID/{params.microbe_bwa_image} "
            "--microbe-dict /lscratch/$SLURM_JOBID/{params.microbe_dict_file} "
            "--taxonomy-file /lscratch/$SLURM_JOBID/{params.taxonomy_db} "
            "--output '{output.pathseq_bam}' "
            "--scores-output '{output.pathseq_output}' "
            "--filter-metrics '{output.filter_metrics}' "
            "--score-metrics '{output.score_metrics}' "
            '--java-options "-Xmx64g -Xms64G -Djava.io.tmpdir=/lscratch/$SLURM_JOBID/tmp -XX:+UseG1GC -XX:ParallelGCThreads=8 -XX:ConcGCThreads=2" '
            '--spark-master local[8] '
            + config["params"]["PathSeq"]
        )
