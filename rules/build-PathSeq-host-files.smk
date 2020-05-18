# params
GATK_VERSION = "4.1.6.0"

# Rules for building host files
rule build_host_kmer_file:
    input:
        config["ref"]["genome"]
    output:
        config["PathSeq"]["host_bfi"]
    shell:
        "module load GATK/{GATK_VERSION} && "
        "gatk PathSeqBuildKmers "
        "--java-options '-Xmx80g' "
        "--reference '{input}' "
        "-O '{output}'"

rule build_host_BWA_image:
    input:
        config["ref"]["genome"]
    output:
        config["PathSeq"]["host_img"]
    shell:
        "module load GATK/{GATK_VERSION} && "
        "gatk BwaMemIndexImageCreator -I {input} -O {output}"
