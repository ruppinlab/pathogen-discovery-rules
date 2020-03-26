from os.path import join


KRAKEN_OUTPUT_FILE = join("output", "Kraken", "{patient}-{sample}", "sequences.kraken")
KRAKEN_TRANSLATE_FILE = join("output", "Kraken", "{patient}-{sample}", "sequences_mpa_report.txt")
KRAKEN_BIOM_FILE = join("output", "Kraken", "{patient}-{sample}", "sequences.biom")

# biowulf kraken help page - https://hpc.nih.gov/apps/kraken.html
# key idea - copy the DB into memory via "cp -r $DB /dev/shm"
# 20180220_standard DB is about 200GB in size
rule run_Kraken:
    params:
        dbname = config["Kraken"]["dbname"] # this is a directory
    input:
        fq1 = expand(join("FASTQ", "raw", "{patient}-{sample}_1.fastq.gz"), zip, patient=samples["patient"], sample=samples["sample"]),
        fq2 = expand(join("FASTQ", "raw", "{patient}-{sample}_2.fastq.gz"), zip, patient=samples["patient"], sample=samples["sample"]),
    output:
        expand(KRAKEN_OUTPUT_FILE, zip, patient=samples["patient"], sample=samples["sample"])
    run:
        shell("module load kraken/1.1")
        shell("trap 'rm -rf /dev/shm/{params.dbname}' EXIT")
        shell("cp -r /fdb/kraken/{params.dbname} /dev/shm")
        for fq1, fq2, o in zip(input.fq1, input.fq2, output):
            shell(
                "kraken --db /dev/shm/{params.dbname} --fastq-input --paired "
                "--gzip-compressed --check-names --output {o} "
                + config["params"]["Kraken"] + " "
                "{fq1} {fq2}"
                )

rule run_Kraken_translate:
    params:
        dbname = config["Kraken"]["dbname"]
    input:
        expand(KRAKEN_OUTPUT_FILE, zip, patient=samples["patient"], sample=samples["sample"])
    output:
        expand(KRAKEN_TRANSLATE_FILE, zip, patient=samples["patient"], sample=samples["sample"])
    run:
        shell("module load kraken/1.1")
        shell("trap 'rm -rf /dev/shm/{params.dbname}' EXIT")
        shell("cp -r /fdb/kraken/{params.dbname} /dev/shm")
        for i, o in zip(input, output):
            shell(
                "kraken-translate --db /dev/shm/{params.dbname} "
                "--mpa-format {i} > {o}"
                )


rule parse_Kraken_to_biom:
    conda:
        "../envs/kraken_biom_parser.yaml"
    input:
        KRAKEN_TRANSLATE_FILE
    output:
        KRAKEN_BIOM_FILE
    shell:
        "../src/parse_kraken_to_biom.py --kraken-translate-report-fp {input} "
        "--taxonomic-rank " + config["Kraken"]["taxonomic_rank"] + " "
        "--biom-output-fp {output}"
