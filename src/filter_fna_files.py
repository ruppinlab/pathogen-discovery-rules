from Bio import SeqIO
from multiprocessing import Pool
import itertools


bacteria_of_interest = snakemake.params["microbes_of_interest"]


def identify_specific_sequences(f):
    output = []
    for record in SeqIO.parse(f, "fasta"):
        # NC - complete genomic molecule, usually reference assembly
        # NZ - Unfinished WGS - https://www.biostars.org/p/209159/
        if record.name.startswith("NC"):
            if any(b in record.description for b in bacteria_of_interest):
                output.append(record)
    return output


if __name__ == '__main__':
    with Pool(31) as p:
        list_of_lists = p.map(identify_specific_sequences, snakemake.input)
        sequences = list(itertools.chain(*list_of_lists))
        SeqIO.write(sequences, snakemake.output[0], "fasta")
