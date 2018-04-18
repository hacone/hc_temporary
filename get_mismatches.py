import sys
import pysam
import re
from collections import Counter
from itertools import takewhile


# p = read.get_aligned_pairs returns the following tuple - p: (query_pos, ref_pos)
# with with_seq, p: (query_pos, ref_pos, ref_base)

b2i = {"A":0, "C":1, "G":2, "T":3}
i2b = {0:"A", 1:"C", 2:"G", 3:"T"}

bam = pysam.AlignmentFile("data/SRR1997411/alignments/SRR1997411.join.aligned.sort.bam")

def normalize(l):
    s = sum(l)
    if s:
        return [li / s for li in l]
    else:
        return [0 for li in l]

def count_coverage(reads, rlen, normalized = False):
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

    # list of 4 lists cc[i][b] TODO: make this to be numpy array of shape (rlen, 4)
    cc = [[counter[(i, b)] for b in "ACGT"] for i in range(rlen)]

    if normalized:
        # normalize each column
        return [ normalize([cc[i][b] for b in range(4)]) for i in range(rlen) ]
    else:
        return cc

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

def openFasta(path):
    """ open fasta as simple dict (refname is trimmed after the first space)"""
    from Bio.SeqIO import FastaIO
    with open(path) as handle:
        # trim after the first space (as in ref in bam file)
        return { item[0].split()[0]:item[1] for item in dict(FastaIO.SimpleFastaParser(handle)).items() }

## two possible implementations; 2nd is more extensive concept (eg, 97, 2, .5, .5)
## 1. fix varsite as maxfreq<.9 then find varbase as freq>.01
## 2. find varsite as site with 2 bases freq>.01
def detect_variant_sites():
    # TODO: use distinct monomers set later
    ref_dict = openFasta("data/monomers/MigaKH.HigherOrderRptMon.fa")

    for i in range(len(bam.references))[:20]:
        ref, rlen = bam.references[i], bam.lengths[i]

        print(f"#####\tref = {ref}\t#####")

        reads = list(bam.fetch(contig=ref))
        nc = count_coverage(reads, rlen, normalized = True)

        variant_pos = [i for i in range(rlen) if max(nc[i]) < .90]
        
        ## here is far from pythonic
        variants = []
        for vp in variant_pos:
            varbase = []
            for j in range(4):
                if nc[vp][j] > .01:
                    varbase += ((i2b[j], nc[vp][j]))
            variants += ((vp,varbase))

        print(variants)
        sys.stdout.flush()

# TODO: rel freq is important as well to feel the error rate of short reads.
def absolute_frequency_distribution():

    counter, max_cov = Counter(), 0

    for i in range(len(bam.references)): # just test in subset
        ref, rlen = bam.references[i], bam.lengths[i]
        reads = list(bam.fetch(contig=ref)) # TODO: do i need to listify this ?
        nc = count_coverage(reads, rlen, normalized = False)

        print(f"{ref}\t{i}\t{len(bam.references)}")
        sys.stdout.flush()
        ## again it's not pythonic
        for pos in range(rlen):
            for base in range(4):
                if max_cov < nc[pos][base]:
                    max_cov = nc[pos][base]
                counter[nc[pos][base]] += 1

    # report !
    for c in range(max_cov + 1):
        print(f"{c}\t{counter[c]}")


if __name__ == '__main__':

    #detect_variant_sites()

    absolute_frequency_distribution()

