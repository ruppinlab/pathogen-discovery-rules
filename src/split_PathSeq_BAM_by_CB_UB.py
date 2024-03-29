import pysam
from collections import defaultdict

# suppress incorrect error warning - https://github.com/pysam-developers/pysam/issues/939
save = pysam.set_verbosity(0)
# load and iterate through the PathSeq BAM file
pathseq_bam = pysam.AlignmentFile(snakemake.input[0], mode="rb")
# set verbosity back to original setting
pysam.set_verbosity(save)

output = []
UMI_dict = defaultdict(list)
# seg is an AlignedSegment object
for seg in pathseq_bam.fetch(until_eof=True):
    # not all records will have the CB tag and the UB tag - they should now
    if seg.has_tag("CB") and seg.has_tag("UB"):
        if (seg.get_tag(tag="CB") == snakemake.wildcards["cell"]):
            UMI_dict[seg.get_tag(tag="UB")].append(seg)


barcode_bam = pysam.AlignmentFile(snakemake.output[0], mode="wb", template=pathseq_bam)
for UMI in UMI_dict:
    #print(UMI)
    # keep one read per UMI - the read with the highest mapping quality
    UMI_reads = UMI_dict[UMI]
    UMI_read = UMI_reads[0]
    #print(UMI_read)
    for read in UMI_reads:
        #print(read)
        if read.mapping_quality > UMI_read.mapping_quality:
            UMI_read = read
    barcode_bam.write(UMI_read)
barcode_bam.close()

pathseq_bam.close()
