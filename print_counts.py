import pysam

sam = pysam.AlignmentFile("SRR3189743.join.aligned.sort.bam")

# print(sam.header)

#sam.count(contig=sam.references[100])
#sam.count_coverage(contig=sam.references[100])
#sam.pileup(contig=sam.references[100])

#for c in sorted([(sam.count(contig=c), c) for c in sam.references]):
#    print(f"{c[0]}\t{c[1]}")

target_contig = sam.references[0]

# select significant column
cc = sam.count_coverage(contig=target_contig)

[ for i in range(len(cc)) ]

# divide reads into new tmp bams

# recur for new bams

sam.close()
