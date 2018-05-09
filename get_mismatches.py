from collections import Counter, defaultdict
from itertools import takewhile
import pysam
import re
import sys

# p = read.get_aligned_pairs returns the following tuple - p: (query_pos, ref_pos)
# with with_seq, p: (query_pos, ref_pos, ref_base)

b2i = {"A":0, "C":1, "G":2, "T":3}
i2b = {0:"A", 1:"C", 2:"G", 3:"T"}

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
## 3. rel.freq must be > error rate = 0.1? (estimate?); abs.freq must be > seq depth (30?)
def detect_variant_sites(bamfile, reffile):
    """ define variant sites from short read alignment data and write out as pickle. """

    bam = pysam.AlignmentFile(bamfile)
    ref = openFasta(reffile) # I don't need this ??
    variant_sites = {} # {"Monomer1":{ "Monomer":freq_as_mon, (pos,base):freq }}

    # SNVs with rel. freq. less than this value are supposed to be errors,
    # or at least it's difficult to distinguish them with error.
    t_err = 0.02
    # SNVs with abs. freq. less than short read sequencing depth are supposed to be errors
    t_abs = 50

    for i in range(len(bam.references)):
        name = bam.references[i]
        reads = list(bam.fetch(contig=name)) # TODO: do i need to listify this ?
        counts = count_coverage(reads, bam.lengths[i], normalized = False)
        vs = defaultdict(int)
        vs["Monomer"] = len(reads)
        
        for idx in range(bam.lengths[i]):
            for r, c in sorted(enumerate(counts[idx]), key=lambda x:x[1])[:-1]:
                if (c > max([t_abs, t_err * len(reads)])):
                    vs[ (idx, i2b[r]) ] = c
        variant_sites[name] = vs

        print(f"{i}\t{name}\t{len(reads)}\t{len(list(vs))}")
        # TODO; i need frequency of variants
        sys.stdout.flush()

    import pickle
    with open(f"variant_sites_{t_abs}_{t_err}.pickle", "wb") as f:
        pickle.dump(variant_sites, f)

# TODO: rel freq is important as well to feel the error rate of short reads.
def absolute_frequency_distribution(bamfile):
    """ absolute frequency distribution grouped by each range of overall monomer coverage"""

    bam = pysam.AlignmentFile(bamfile)
    # read depth of each monomer
    mon_depth = {}
    # counters for each monomer, each rank (1-4).
    counters = [ { j : Counter() for j in range(4) } for i in range(len(bam.references)) ]
    max_freq = [ { j : 0 for j in range(4) } for i in range(len(bam.references)) ]

    for i in range(len(bam.references)):
        name = bam.references[i]
        reads = list(bam.fetch(contig=name)) # TODO: do i need to listify this ?
        mon_depth[i] = len(reads)
        counts = count_coverage(reads, bam.lengths[i], normalized = False)

        for idx in range(bam.lengths[i]):
            f = sorted(counts[idx])
            for r in range(4):
                counters[i][r][f[r]] += 1
                max_freq[i][r] = max(max_freq[i][r], f[r]) 

        print(f"{i}\t{name}\t{len(reads)}\t{max_freq[i]}")
        sys.stdout.flush()

    
    # report! for each depth range, <1000, <2000, <3000, ...
    thresholds = list(range(0, 20000, 1000)) + list(range(20000, 40000, 2000)) + list(range(40000, 100000, 10000))

    for t in range(len(thresholds) - 1):
        f = open(f"abs_freq_dist.{thresholds[t+1]}.dat", "w")
        counts_bin = [ sum( [ counters[i][r] for i in range(len(bam.references))
            if thresholds[t] <= mon_depth[i] & mon_depth[i] < thresholds[t+1] ],
            Counter()) for r in range(4) ]
        max_bin = [ max( [0] + [ max_freq[i][r] for i in range(len(bam.references))
                     if thresholds[t] <= mon_depth[i] & mon_depth[i] < thresholds[t+1] ]) for r in range(4) ]

        # for each rank (least frequent to most frequent)
        for r in range(4):
            for c in range(max_bin[r] + 1):
                if counts_bin[r][c] > 0:
                    f.write(f"{r}\t{c}\t{counts_bin[r][c]}\n")
            f.write("\n\n")

        # all ranks aggregated
        all_c = sum(counts_bin, Counter())
        for c in range(max(max_bin) + 1):
            if all_c[c] > 0:
                f.write(f"ALL\t{c}\t{all_c[c]}\n")
        f.close()


if __name__ == '__main__':

    # TODO: write menu
    import argparse
    parser = argparse.ArgumentParser(description='Define SNVs from short-read alignments.')
    parser.add_argument('action', metavar='action', type=str, help='action to perform: define-vars, freq-vars, ...')
    parser.add_argument('--bam', dest='bamfile', help='path to BAM format file')
    parser.add_argument('--ref-fa', dest='reffile', help='path to reference .fa file')
    args = parser.parse_args()

    #bamfile = "data/SRR1997411/alignments/SRR1997411.join.aligned.sort.bam"
    #reffile = "data/monomers/MigaKH.HigherOrderRptMon.fa"

    if args.action == "define-vars":
        assert args.bamfile, "bam file is missing"
        assert args.reffile, "ref file is missing"
        detect_variant_sites(args.bamfile, args.reffile)
    elif args.action == "freq-vars":
        assert args.bamfile, "bam file is missing"
        absolute_frequency_distribution(args.bamfile)

