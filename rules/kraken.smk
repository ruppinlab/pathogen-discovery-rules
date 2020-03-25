

# biowulf kraken help page - https://hpc.nih.gov/apps/kraken.html
# key idea - copy the DB into memory via "cp -r $DB /dev/shm"
# 20180220_standard DB is about 200GB in size
rule run_kraken:
    params:
        dbname = "20180220_standard_8GB" # this is a directory
    input:
        fq1 = expand(join("FASTQ", "raw", "{patient}-{sample}_1.fastq.gz"), zip, patient=samples["patient"], sample=samples["sample"]),
        fq2 = expand(join("FASTQ", "raw", "{patient}-{sample}_2.fastq.gz"), zip, patient=samples["patient"], sample=samples["sample"]),
    output:
        expand(join("output", "kraken", "{patient}-{sample}.kraken"))
    run:
        shell("module load kraken/1.1")
        shell("trap 'rm -rf /dev/shm/{params.dbname}' EXIT")
        shell("cp -r /fdb/kraken/{params.dbname} /dev/shm")
        for fq1, fq2, o in zip(input.fq1, input.fq2, output):
            shell(
                "kraken --db /dev/shm/{params.dbname} --fastq-input --paired "
                "--gzip-compressed --check-names --threads $SLURM_CPUS_PER_TASK "
                "--only-classified-output --output {o} {fq1} {fq2}"
                )
