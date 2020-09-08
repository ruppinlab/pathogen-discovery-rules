import pysam


# load and index the CellRanger BAM file
in_bam = pysam.AlignmentFile(snakemake.input[0], mode="rb")
out_bam = pysam.AlignmentFile(snakemake.output[0], mode="wb", template=in_bam)
# seg is an AlignedSegment object
for seg in in_bam:

    # check if seg has polyA (pa) or TSO (ts) tag?
    if seg.has_tag("pa"):
        # there is a polyA tail that we need to trim from the 3' end of the read
        query_qual = seg.query_qualities
        polya_tail_len = seg.get_tag(tag="pa")
        # hard clip the poly-A tail on the 3' end of the read
        stop_trim = seg.query_length-polya_tail_len # where does the polya tail stop?
        seg.query_sequence = seg.query_sequence[0:stop_trim]
        seg.query_qualities = query_qual[0:stop_trim]
    if seg.has_tag("ts"):
        # there is a TSO that we need to trim from the 5' end of the read
        TSO_len = seg.get_tag(tag="ts")
        query_qual = seg.query_qualities
        seg.query_sequence = seg.query_sequence[TSO_len:]
        seg.query_qualities = query_qual[TSO_len:]
        # hard clip the TSO from the 5' end of the read
    if (seg.query_length > 15) and (seg.has_tag("CB")) and (seg.has_tag("UB")):
        out_bam.write(seg)

out_bam.close()
