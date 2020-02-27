from Bio import SeqIO

for fna_file in snakemake.input:
    for record in SeqIO.parse(fna_file, "fasta"):
        if record.description.contains("Salmonella enterica"):
            print(record.description)
