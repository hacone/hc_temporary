import sys
import pysam
import re
from collections import Counter
from itertools import takewhile

bam = pysam.AlignmentFile("data/SRR1997411/alignments/SRR1997411.join.aligned.sort.bam")

def count_coverage(reads, rlen):
    """
        emulate pysam.AlignmentFile.count_coverage for set of reads given as iterator
        avoid to create temporary bam file.
        reads must be aligned to the same contig. TODO; add assertion
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


### extract neighboring information
def cigar_to_alnstart(cigarstr):
    # TODO: let's check if I can remove list()
    cl = list(filter(None, re.split(r"([0-9]+[A-Z])", cigarstr))) # cl = ['51M', '301S']
    ct = [ re.split(r"([0-9]+)", c)[1:] for c in cl ] # ct = [['51', 'M'], ['301', 'S']
    return sum([ int(a[0]) for a in takewhile(lambda e: e[1] != 'M', ct) ])

def parse_sa_entry(str, cigar_only = False):
    """
    parse a part of SA tag which state one alignment.
    intended usage ex.: [ parse_sa_entry(str) for str in filter(None, read.get_tag("SA").split(";")) ]

    input: a part of SA tag (each separated by ;)
    output: [alt_contig_name, alt_alnstart_inread, ori('+'/'-')]
    """
    strsp = str.split(",")
    if cigar_only:
        return [strsp[0], strsp[3], strsp[2]]
    else:
        return [strsp[0], cigar_to_alnstart(strsp[3]), strsp[2]]

def get_neighbors(read, cigar_only = False):
    """
    input: pysam.AlignedSegment
    output: [[alt_contig_name, alt_alnstart_inread, ori]] <- TODO: This should be namedtuple
    """
    if not read.has_tag("SA"):
        return []
    else:
        #return [ parse_sa_entry(str) for str in filter(None, read.get_tag("SA").split(";")) ]
        return sorted([ parse_sa_entry(str, cigar_only) for str in filter(None, read.get_tag("SA").split(";")) ],
                      key=lambda x: x[1])
                  
    #if alt_start < start: before = alt_mon...
    # or make list of (mon, start) then split at start.

if __name__ == '__main__out':

    nref = len(bam.references)
    ref_to_i = { bam.references[i] : i for i in range(nref) }

    import numpy as np
    pred_table = np.zeros((nref, nref))
    next_table = np.zeros((nref, nref))

    for ref in bam.references:
        reads = list(bam.fetch(contig=ref))
        pred_counter = Counter()
        next_counter = Counter()

        for read in reads:
            ns = get_neighbors(read)
            alnstart = cigar_to_alnstart(read.cigarstring)

            preds = sorted([n for n in ns if n[1] < alnstart], key=lambda x: x[1])
            nexts = sorted([n for n in ns if n[1] > alnstart], key=lambda x: x[1])

            #print(f"read - {read.cigarstring} aln from {alnstart}")
            #print(f"neighbors: {ns}")
            #print(f"preds: {preds}")
            #print(f"nexts: {nexts}\n")

            if preds:
                # print(f"preds: {preds[-1]}")
                pred_counter.update([preds[-1][0]])
            if nexts:
                # print(f"nexts: {nexts[0]}\n")
                next_counter.update([nexts[0][0]])

        print(f"\tFor {ref} - ")
        print(f"\tpred_counter({sum(pred_counter.values())} elems) = {pred_counter.most_common(7)}")
        print(f"\tnext_counter({sum(next_counter.values())} elems) = {next_counter.most_common(7)}\n\n")
        sys.stdout.flush()

        for item in pred_counter.items():
            pred_table[ref_to_i[ref], ref_to_i[item[0]]] = item[1]
        for item in next_counter.items():
            next_table[ref_to_i[ref], ref_to_i[item[0]]] = item[1]

        print(pred_table)
        print(next_table)

    np.save("preds.npy", pred_table)
    np.save("nexts.npy", next_table)

def normalize(l):
    s = sum(l)
    if s:
        return [li / s for li in l]
    else:
        return [0 for li in l]

if __name__ == '__main__':

    b2i = {"A":0, "C":1, "G":2, "T":3}

    def openFasta(path):
        """ open fasta as simple dict """
        from Bio.SeqIO import FastaIO
        with open(path) as handle:
            # trim after the first space (as in ref in bam file)
            return { item[0].split()[0]:item[1] for item in dict(FastaIO.SimpleFastaParser(handle)).items() }

    ref_dict = openFasta("data/monomers/MigaKH.HigherOrderRptMon.fa")
    for i in range(len(bam.references))[:10]:
        ref = bam.references[i]
        rlen = bam.lengths[i]

        print(f"#####\tref = {ref}\t#####")
        # print(bam.count_coverage(contig=ref))
        reads = list(bam.fetch(contig=ref))
        cc = count_coverage(reads, rlen)

        # [i,b]
        nc = [normalize([cc[b][i] for b in range(4)]) for i in range(rlen)]

        # list every position with freq(ref) < 99% or max(freq(b)) < 99%
        non_ref_pos = [i for i in range(rlen) if nc[i][b2i[ref_dict[ref][i]]] < .50 ]
        variant_pos = [i for i in range(rlen) if max(nc[i]) < .90]

        for i in non_ref_pos:
            print((i, nc[i][b2i[ref_dict[ref][i]]]))

        for i in variant_pos:
            print((i, nc[i]))

        #print(nc)
        sys.stdout.flush()

