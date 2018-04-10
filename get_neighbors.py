import sys
import pysam
import re
from collections import Counter
from itertools import takewhile

# p = read.get_aligned_pairs returns the following tuple - p: (query_pos, ref_pos)
# with with_seq, p: (query_pos, ref_pos, ref_base)

bam = pysam.AlignmentFile("data/SRR1997411/alignments/SRR1997411.join.aligned.sort.bam")

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

if __name__ == '__main':

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
