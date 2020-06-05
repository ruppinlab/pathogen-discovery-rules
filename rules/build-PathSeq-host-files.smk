# params
GATK_VERSION = "4.1.6.0"

# Rules for building host files

# builds a hash table of host reference k-mers of length 31 where there is a
# masked base at position 16 in each k-mer
# the other supported option is to build a Bloom filter, which is faster and
# uses less memory but has a higher probability of false positives
rule build_host_kmer_file:
    input:
        config["ref"]["genome"]
    output:
        config["PathSeq"]["host_bfi"]
    shell:
        "module load GATK/{GATK_VERSION} && "
        "gatk PathSeqBuildKmers "
        "--java-options '-Xmx80g' "
        "--kmer-mask 16 --kmer-size 31 "
        "--reference '{input}' "
        "-O '{output}'"

# build a BWA-MEM index image file for use with GATK BWA tools
rule build_host_BWA_image:
    input:
        config["ref"]["genome"]
    output:
        config["PathSeq"]["host_img"]
    shell:
        "module load GATK/{GATK_VERSION} && "
        "gatk BwaMemIndexImageCreator -I {input} -O {output}"

# build a sequence dictionary for a reference sequence
# output is a SAM/BAM file that contains a header but no SAMRecords
# rule build_host_reference_dict:
#     input:
#         config["ref"]["genome"]
#     output:
#         config["PathSeq"]["host_dict"]
#     shell:
#         "module load GATK/{GATK_VERSION} && "
#         "gatk CreateSequenceDictionary -R {input}"
# 
