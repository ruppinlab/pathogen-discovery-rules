

# A Kraken DB is a directory containing at least 4 files
# database.kdb
# database.idx
# taxonomy/nodes.dmp
# taxonomy/names.dmp

# to use the NCBI taxonomy - not sure there's another way - hopefully this works
# kraken-build --download-taxonomy --db $DBNAME
# to add your own fa file, use
# kraken-build --add-to-library *.fa --db $DBNAME
rule build_Kraken_DB:
    conda:
        "../envs/kraken.yml"
    params:
        dir = join(config["Kraken"]["db_path"], config["Kraken"]["dbname"])
    input:
        fasta=config["Kraken"]["microbe_fasta"]
    output:
        db=join(config["Kraken"]["db_path"], config["Kraken"]["dbname"], "database.kdb")
    shell:
        "kraken-build --download-taxonomy --db {output.db} && "
        "kraken-build --add-to-library {input.fasta} --db {output.db} && "
        "kraken-build --build --db {output.db} " + config["params"]["Kraken"]
