# input FA files
HUMAN_FA = "/data/Robinson-SB/PathSeq-data/pathseq_host.fa"
ERCC_FA = "/data/Robinson-SB/PathSeq-data/ERCC92.fa"
RNA_SPIKE_FA = "/data/Robinson-SB/PathSeq-data/RNA_spike.fa"
# combined FA files
HUMAN_ERCC_FA = "/data/Robinson-SB/PathSeq-data/pathseq_host_ERCC92.fa"
HUMAN_RNA_SPIKE_FA = "/data/Robinson-SB/PathSeq-data/pathseq_host_RNA_spike.fa"
# bloom filters
HUMAN_ERCC_BFI = "/data/Robinson-SB/PathSeq-data/pathseq_host_ERCC92.bfi"
HUMAN_RNA_SPIKE_BFI = "/data/Robinson-SB/PathSeq-data/pathseq_host_RNA_spike.bfi"
# dictionaries
HUMAN_ERCC_DICT = "/data/Robinson-SB/PathSeq-data/pathseq_host_ERCC92.dict"
HUMAN_RNA_SPIKE_DICT = "/data/Robinson-SB/PathSeq-data/pathseq_host_RNA_spike.dict"
# fasta index files
HUMAN_ERCC_FAI = "/data/Robinson-SB/PathSeq-data/pathseq_host_ERCC92.fa.fai"
HUMAN_RNA_SPIKE_FAI = "/data/Robinson-SB/PathSeq-data/pathseq_host_RNA_spike.fa.fai"

localrules: combine_human_ERCC_FA, combine_human_RNA_SPIKE_FA
localrules: build_host_ERCC_reference_dict, build_host_RNA_SPIKE_reference_dict
localrules: build_host_ERCC_fai, build_host_RNA_SPIKE_fai

rule all:
    input:
        HUMAN_ERCC_BFI,
        HUMAN_RNA_SPIKE_BFI,
        HUMAN_ERCC_DICT,
        HUMAN_RNA_SPIKE_DICT,
        HUMAN_ERCC_FAI,
        HUMAN_RNA_SPIKE_FAI

# Rules for building host files
rule combine_human_ERCC_FA:
    input:
        HUMAN_FA,
        ERCC_FA
    output:
        HUMAN_ERCC_FA
    shell:
        "cat {input[0]} {input[1]} > {output}"

rule combine_human_RNA_SPIKE_FA:
    input:
        HUMAN_FA,
        RNA_SPIKE_FA
    output:
        HUMAN_RNA_SPIKE_FA
    shell:
        "cat {input[0]} {input[1]} > {output}"

# builds a hash table of host reference k-mers of length 31 where there is a
# masked base at position 16 in each k-mer
# the other supported option is to build a Bloom filter, which is faster and
# uses less memory but has a higher probability of false positives
rule build_ERCC_host_kmer_file:
    input:
        HUMAN_ERCC_FA
    output:
        HUMAN_ERCC_BFI
    shell:
        "module load GATK/4.1.8.1 && "
        "gatk PathSeqBuildKmers "
        "--java-options '-Xmx80g' "
        "--bloom-false-positive-probability 0.001 "
        "--kmer-mask 16 --kmer-size 31 "
        "--reference '{input}' "
        "-O '{output}'"

rule build_RNA_spike_host_kmer_file:
    input:
        HUMAN_RNA_SPIKE_FA
    output:
        HUMAN_RNA_SPIKE_BFI
    shell:
        "module load GATK/4.1.8.1 && "
        "gatk PathSeqBuildKmers "
        "--java-options '-Xmx80g' "
        "--bloom-false-positive-probability 0.001 "
        "--kmer-mask 16 --kmer-size 31 "
        "--reference '{input}' "
        "-O '{output}'"

# build a BWA-MEM index image file for use with GATK BWA tools
rule build_host_RNA_SPIKE_BWA_image:
    input:
        HUMAN_RNA_SPIKE_FA
    output:
        HUMAN_RNA_SPIKE_IMG
    shell:
        "module load GATK/4.1.8.1 && "
        "gatk BwaMemIndexImageCreator -I {input} -O {output}"

rule build_host_ERCC_BWA_image:
    input:
        HUMAN_ERCC_FA
    output:
        HUMAN_ERCC_IMG
    shell:
        "module load GATK/4.1.8.1 && "
        "gatk BwaMemIndexImageCreator -I {input} -O {output}"

# build a sequence dictionary for a reference sequence
# output is a SAM/BAM file that contains a header but no SAMRecords
rule build_host_ERCC_reference_dict:
    input:
        HUMAN_ERCC_FA
    output:
        HUMAN_ERCC_DICT
    shell:
        "module load GATK/4.1.8.1 && "
        "gatk CreateSequenceDictionary -R {input}"

rule build_host_RNA_SPIKE_reference_dict:
    input:
        HUMAN_RNA_SPIKE_FA
    output:
        HUMAN_RNA_SPIKE_DICT
    shell:
        "module load GATK/4.1.8.1 && "
        "gatk CreateSequenceDictionary -R {input}"

rule build_host_ERCC_fai:
    input:
        HUMAN_ERCC_FA
    output:
        HUMAN_ERCC_FAI
    shell:
        "module load samtools && "
        "samtools faidx {input}"

rule build_host_RNA_SPIKE_fai:
    input:
        HUMAN_RNA_SPIKE_FA
    output:
        HUMAN_RNA_SPIKE_FAI
    shell:
        "module load samtools && "
        "samtools faidx {input}"
