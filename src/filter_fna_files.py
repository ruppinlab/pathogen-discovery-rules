from Bio import SeqIO
from multiprocessing import Pool
import itertools


bacteria_of_interest = ["Salmonella enterica subsp. enterica serovar Typhimurium str. SL1344", "Salmonella enterica subsp. enterica serovar Typhimurium str. LT2"]


def identify_specific_sequences(f):
    output = []
    for record in SeqIO.parse(f, "fasta"):
        if any(b in record.description for b in bacteria_of_interest):
            output.append(record)
    return output


if __name__ == '__main__':
    with Pool(31) as p:
        list_of_lists = p.map(identify_specific_sequences, snakemake.input)
        flat_list = list(itertools.chain(*list_of_lists)) 
        for r in flat_list:
            print(r.description)
