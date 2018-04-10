import sys
import pysam
import re
from collections import Counter
from itertools import takewhile

# p = read.get_aligned_pairs returns the following tuple - p: (query_pos, ref_pos)
# with with_seq, p: (query_pos, ref_pos, ref_base)

bam = pysam.AlignmentFile("data/SRR1997411/alignments/SRR1997411.join.aligned.sort.bam")

def count_coverage(reads, rlen):
    """
        Emulate pysam.AlignmentFile.count_coverage for a set of reads given as iterator.
        This avoids to create temporary bam file.
        Reads must be aligned to the same contig. TODO; add (optional) assertion.
    """
    counter = Counter()
    for read in reads:
        qseq = read.query_sequence
        for p in read.get_aligned_pairs(matches_only=True):
            if read.query_qualities[p[0]] > 14: # needed to emulate pysam impl
                counter[(p[1], qseq[p[0]])] += 1
    return [[counter[(i, b)] for i in range(rlen)] for b in "ACGT"]

def TODO_FUNC():
    # wont need this
    qaseq = reads[0].query_alignment_sequence 
    qseq = reads[1].query_sequence
    # pos_in_ref, base_in_ref(lowercase), base_in_read
    mms = [ (p[1], p[2], qseq[p[0]]) for p in reads[1].get_aligned_pairs(matches_only=True, with_seq=True) if p[2].islower() ] 

    ## determine significant ref poss; get mismatch list; count freq of each set; perform grad-asc algorithm
    # algorithm; enumerate every possible state; their association to observations; initialize;
    # calc likelihood factor; calc gradient; update estimate; loop..
    # report final estimation
    ## Or; count_coverage from list(read); split read; recur
    # read list would have difining vars set; list may be overlapping


def normalize(l):
    s = sum(l)
    if s:
        return [li / s for li in l]
    else:
        return [0 for li in l]

if __name__ == '__main__':

    b2i = {"A":0, "C":1, "G":2, "T":3}
    i2b = {0:"A", 1:"C", 2:"G", 3:"T"}

    def openFasta(path):
        """ open fasta as simple dict """
        from Bio.SeqIO import FastaIO
        with open(path) as handle:
            # trim after the first space (as in ref in bam file)
            return { item[0].split()[0]:item[1] for item in dict(FastaIO.SimpleFastaParser(handle)).items() }

    ref_dict = openFasta("data/monomers/MigaKH.HigherOrderRptMon.fa")
    for i in range(len(bam.references))[:20]:
        ref = bam.references[i]
        rlen = bam.lengths[i]

        print(f"#####\tref = {ref}\t#####")
        # print(bam.count_coverage(contig=ref))
        reads = list(bam.fetch(contig=ref))
        cc = count_coverage(reads, rlen)

        # [i,b]
        nc = [normalize([cc[b][i] for b in range(4)]) for i in range(rlen)]

        # list every position with freq(ref) < 99% or max(freq(b)) < 99%
        # non_ref_pos = [i for i in range(rlen) if nc[i][b2i[ref_dict[ref][i]]] < .50 ]
        non_ref_pos = []
        variant_pos = [i for i in range(rlen) if max(nc[i]) < .90]

        for i in non_ref_pos:
            print((i, nc[i][b2i[ref_dict[ref][i]]]))

        for i in variant_pos:
            print(f"{i}:\t{ sorted(nc[i]) }")

        #print(nc)
        sys.stdout.flush()
