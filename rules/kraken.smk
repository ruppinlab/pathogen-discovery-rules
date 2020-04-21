from os.path import join, dirname

# input files
PAIRED_FILTERED_BAM = join("output", "PathSeq", "{patient}-{sample}", "filtered-paired.bam")
UNPAIRED_FILTERED_BAM = join("output", "PathSeq", "{patient}-{sample}", "filtered-unpaired.bam")

# intermediate files
SORTED_PAIRED_FILTERED_BAM = join("output", "PathSeq", "{patient}-{sample}", "sorted-filtered-paired.bam")
PAIRED_FILTERED_FQ1 = join("output", "Kraken", "{patient}-{sample}", "filtered-paired.fq1")
PAIRED_FILTERED_FQ2 = join("output", "Kraken", "{patient}-{sample}", "filtered-paired.fq2")
UNPAIRED_FILTERED_FQ = join("output", "Kraken", "{patient}-{sample}", "filtered-unpaired.fq")

# output files
KRAKEN_PAIRED_OUTPUT_FILE = join("output", "Kraken", "{patient}-{sample}", "paired-sequences.kraken")
KRAKEN_UNPAIRED_OUTPUT_FILE = join("output", "Kraken", "{patient}-{sample}", "unpaired-sequences.kraken")
KRAKEN_PAIRED_TRANSLATE_FILE = join("output", "Kraken", "{patient}-{sample}", "paired-sequences_mpa_report.txt")
KRAKEN_UNPAIRED_TRANSLATE_FILE = join("output", "Kraken", "{patient}-{sample}", "unpaired-sequences_mpa_report.txt")
KRAKEN_PAIRED_BIOM_FILE = join("output", "Kraken", "{patient}-{sample}", "paired-sequences.biom")
KRAKEN_UNPAIRED_BIOM_FILE = join("output", "Kraken", "{patient}-{sample}", "unpaired-sequences.biom")
# scripts
KRAKEN_TO_BIOM_SCRIPT = join(dirname(srcdir("kraken.smk")), "..", "src", "parse_kraken_to_biom.py")

localrules: paired_bam_to_fastq, unpaired_bam_to_fastq

include: "PathSeq.smk"

# rule sort_bam_by_queryname:
#     input:
#         PAIRED_FILTERED_BAM
#     output:
#         temp(SORTED_PAIRED_FILTERED_BAM)
#     shell:
#         "module load samtools && "
#         "samtools sort -n {input} {output}"

rule paired_bam_to_fastq:
    input:
        PAIRED_FILTERED_BAM
    output:
        fq1=temp(PAIRED_FILTERED_FQ1),
        fq2=temp(PAIRED_FILTERED_FQ2)
    shell:
        "module load bedtools && "
        "bedtools bamtofastq -i {input} -fq {output.fq1} -fq2 {output.fq2}"

rule unpaired_bam_to_fastq:
    input:
        UNPAIRED_FILTERED_BAM
    output:
        temp(UNPAIRED_FILTERED_FQ)
    shell:
        "module load bedtools && "
        "bedtools bamtofastq -i {input} -fq {output}"

# biowulf kraken help page - https://hpc.nih.gov/apps/kraken.html
# key idea - copy the DB into memory via "cp -r $DB /dev/shm"
# 20180220_standard DB is about 200GB in size
# run Kraken for paired-end reads filtered by PathSeqSparkFilter
rule run_Kraken_paired_reads:
    params:
        dbname = config["Kraken"]["dbname"] # this is a directory
    input:
        fq1 = expand(PAIRED_FILTERED_FQ1, zip, patient=samples["patient"], sample=samples["sample"]),
        fq2 = expand(PAIRED_FILTERED_FQ2, zip, patient=samples["patient"], sample=samples["sample"]),
        db = join(config["Kraken"]["db_path"], config["Kraken"]["dbname"], "database.kdb")
    output:
        expand(KRAKEN_PAIRED_OUTPUT_FILE, zip, patient=samples["patient"], sample=samples["sample"])
    run:
        shell("trap 'rm -rf /dev/shm/{params.dbname}' EXIT")
        shell("cp -r /fdb/kraken/{params.dbname} /dev/shm")
        for fq1, fq2, o in zip(input.fq1, input.fq2, output):
            shell(
                "module load kraken/1.1 && "
                "kraken --db /dev/shm/{params.dbname} --fastq-input --paired "
                "--check-names --output {o} "
                + config["params"]["Kraken"] + " "
                "{fq1} {fq2}"
                )

rule run_Kraken_unpaired_reads:
    params:
        dbname = config["Kraken"]["dbname"] # this is a directory
    input:
        fq = expand(UNPAIRED_FILTERED_FQ, zip, patient=samples["patient"], sample=samples["sample"]),
        db = join(config["Kraken"]["db_path"], config["Kraken"]["dbname"], "database.kdb")
    output:
        expand(KRAKEN_UNPAIRED_OUTPUT_FILE, zip, patient=samples["patient"], sample=samples["sample"])
    run:
        shell("trap 'rm -rf /dev/shm/{params.dbname}' EXIT")
        shell("cp -r /fdb/kraken/{params.dbname} /dev/shm")
        for fq, o in zip(input.fq, output):
            shell(
                "module load kraken/1.1 && "
                "kraken --db /dev/shm/{params.dbname} --fastq-input "
                "--check-names --output {o} "
                + config["params"]["Kraken"] + " "
                "{fq}"
                )

rule run_Kraken_translate_paired_reads:
    params:
        dbname = config["Kraken"]["dbname"]
    input:
        expand(KRAKEN_PAIRED_OUTPUT_FILE, zip, patient=samples["patient"], sample=samples["sample"])
    output:
        expand(KRAKEN_PAIRED_TRANSLATE_FILE, zip, patient=samples["patient"], sample=samples["sample"])
    run:
        shell("trap 'rm -rf /dev/shm/{params.dbname}' EXIT")
        shell("cp -r /fdb/kraken/{params.dbname} /dev/shm")
        for i, o in zip(input, output):
            shell(
                "module load kraken/1.1 && "
                "kraken-translate --db /dev/shm/{params.dbname} "
                "--mpa-format {i} > {o}"
                )

rule run_Kraken_translate_unpaired_reads:
    params:
        dbname = config["Kraken"]["dbname"]
    input:
        expand(KRAKEN_UNPAIRED_OUTPUT_FILE, zip, patient=samples["patient"], sample=samples["sample"])
    output:
        expand(KRAKEN_UNPAIRED_TRANSLATE_FILE, zip, patient=samples["patient"], sample=samples["sample"])
    run:
        shell("trap 'rm -rf /dev/shm/{params.dbname}' EXIT")
        shell("cp -r /fdb/kraken/{params.dbname} /dev/shm")
        for i, o in zip(input, output):
            shell(
                "module load kraken/1.1 && "
                "kraken-translate --db /dev/shm/{params.dbname} "
                "--mpa-format {i} > {o}"
                )

rule parse_paired_Kraken_to_biom:
    conda:
        "../envs/kraken_biom_parser.yaml"
    input:
        KRAKEN_PAIRED_TRANSLATE_FILE
    output:
        KRAKEN_PAIRED_BIOM_FILE
    shell:
        "python {KRAKEN_TO_BIOM_SCRIPT} --kraken-translate-report-fp {input} "
        "--taxonomic-rank " + config["Kraken"]["taxonomic_rank"] + " "
        "--biom-output-fp {output}"

rule parse_unpaired_Kraken_to_biom:
    conda:
        "../envs/kraken_biom_parser.yaml"
    input:
        KRAKEN_UNPAIRED_TRANSLATE_FILE
    output:
        KRAKEN_UNPAIRED_BIOM_FILE
    shell:
        "python {KRAKEN_TO_BIOM_SCRIPT} --kraken-translate-report-fp {input} "
        "--taxonomic-rank " + config["Kraken"]["taxonomic_rank"] + " "
        "--biom-output-fp {output}"
