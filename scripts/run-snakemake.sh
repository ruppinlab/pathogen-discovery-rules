#!/bin/bash

module load snakemake/5.13.0
snakemake \
--snakefile rules/build-PathSeq-host-files.smk \
--use-conda \
--nolock \
--rerun-incomplete \
--cluster-config config/cluster.json \
--cluster "sbatch --partition={cluster.partition} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.nthreads} --gres=lscratch:{cluster.gres}" \
--jobs 10 \
--latency-wait 60 \
--keep-going \
--local-cores 4 all
