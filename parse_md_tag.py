import pysam
bam = pysam.AlignmentFile("SRR3189743.join.aligned.sort.bam")
target_contig = bam.references[0]
cc = bam.count_coverage(contig=target_contig)



'''
max([1, 2, 3, 2.3])
max([1, 2, 3, 2.3]) / 4
[ max([ for b in range(3)] )]

[ max([ cc[b][i] for b in range(4) ]) for i in range(10) ]
sum([1, 2, 3])
[ max([ cc[b][i] for b in range(4) ]) / sum([ cc[b][i] for b in range(4) ]) for i in range(10) ]
[ sorted([ cc[b][i] for b in range(4) ]) for i in range(10) ]
[ sorted([ cc[b][i] for b in range(4) ]) for i in range(len(cc[0])) ]
[ sorted([ cc[b][i] for b in range(4) ])[3] for i in range(len(cc[0])) ]
[ sorted([ cc[b][i] for b in range(4) ])[0] for i in range(len(cc[0])) ]
[ sorted([ cc[b][i] for b in range(4) ])[2] for i in range(len(cc[0])) ]
[ (i, sorted([ cc[b][i] for b in range(4) ])[2]) for i in range(len(cc[0])) ]
target_contig = bam.references[1]
cc = bam.count_coverage(contig=target_contig)
[ (i, sorted([ cc[b][i] for b in range(4) ])[2]) for i in range(len(cc[0])) ]
[ (i, sorted([ cc[b][i] for b in range(4) ])[3]) for i in range(len(cc[0])) ]
target_contig = bam.references[2]
cc = bam.count_coverage(contig=target_contig)
[ (i, sorted([ cc[b][i] for b in range(4) ])[3]) for i in range(len(cc[0])) ]
[ (i, sorted([ cc[b][i] for b in range(4) ])[2]) for i in range(len(cc[0])) ]
'''
