## This is the class representing a set of monomer-encoded PacBio reads.
from collections import Counter
from collections import namedtuple
from itertools import chain
from itertools import groupby
from itertools import islice
from itertools import zip_longest
import numpy as np
import pickle
import pysam
import re
import sys

# These are all immutable
SNV = namedtuple("SNV", ("pos", "base")) # Int, char
Monomer = namedtuple("Monomer", ("name", "snvs")) # string, list(SNV) or id: Int, list(SNV)
# TODO: should it include blast score? orientation of assignment? aligned parts info?
AssignedMonomer = namedtuple("AssignedMonomer", ("begin", "end", "monomer", "ori")) # Int, Int, Monomer, '+'/'-'
EncodedRead = namedtuple("EncodedRead", ("name", "mons", "length")) # string, list(AssignedMonomer), Int

def pairwise_encoded(er1, er2):
    # banded (real coord) dp

    def match(m1, m2):
        """ (mis)match score for pair of monomers. """
        if m1.name != m2.name:
            return -10000
        else:
            s1, s2 = { s for s in m1.snvs }, { s for s in m2.snvs }
            return len(s1 & s2) - len(s1 ^ s2)

    def backtrace(s, bt, i, j):
        r = []
        while (i > 0) & (j > 0):
            r = r + [(i, j)]
            if bt[i,j] == 0:
                i -= 1
            elif bt[i,j] == 1:
                j -= 1
            else:
                i -= 1
                j -= 1
        return sorted(r)


    valid_alignments = []
    l1, l2 = len(er1.mons), len(er2.mons)

    if len({ m.monomer.name for m in er1.mons } & { m.monomer.name for m in er2.mons }) < 4:
        return valid_alignments # not align
        #return None # not align


    # init # NOTE: I need dovetail alignments!
    # i_begin = 0, j_begin in range(l2)
    i_begin = 1 # NOTE: always start with 1
    # for j_begin in [ j for j in range(1,l2+1) if er1.mons[i_begin-1].monomer.name == er2.mons[j-1].monomer.name ]:
    # TODO: maybe I won't need alignments starting with gap-gap
    for j_begin in [ j for j in range(1,l2) if er1.mons[i_begin-1].monomer.name == er2.mons[j-1].monomer.name ]:
        print(f"i_begin = {i_begin}")
        print(f"j_begin = {j_begin}")
        s = np.zeros((l1+1)*(l2+1)).reshape(l1+1, l2+1)
        bt = np.zeros((l1+1)*(l2+1)).reshape(l1+1, l2+1)

        for i in range(l1+1): 
            s[i,0] = -10000
        for j in range(l2+1):
            s[0,j] = -10000
        for i in range(l1+1):
            s[i, j_begin-1] = -10000

        s[i_begin-1, j_begin-1] = 0 # This ensures the first monomer match

        # recur
        for i in range(1, l1+1):
            ioffset = er1.mons[i_begin-1].end - er1.mons[i-1].end
            for j in range(j_begin, l2+1):
                joffset = er2.mons[j_begin-1].end - er2.mons[j-1].end
                if abs(ioffset - joffset) > 800: # NOTE: won't allow displacement of ~5 monomers
                    s[i,j] = -10000
                else:
                    if er1.mons[i-1].monomer.name == "GAP":
                        ii = s[i-1,j] # ins in i
                    else:
                        ii = s[i-1,j]-10 

                    if er1.mons[j-1].monomer.name == "GAP":
                        ij = s[i,j-1] # ins in j
                    else:
                        ij = s[i,j-1]-10 

                    # match
                    m = s[i-1,j-1] + match(er1.mons[i-1].monomer, er2.mons[j-1].monomer)
                    s[i,j], bt[i,j] = max([(ii, 0), (ij, 1), (m, 2)])

        #print(s)
        #print(bt)
        edges = [ (s[l1, j], l1, j) for j in range(l2) ] + [ (s[i, l2], i, l2) for i in range(l1) ] + [(s[l1, l2], l1, l2)]
        max_s, max_i, max_j = max(edges)
        if max_s > 0:
            valid_alignments += [backtrace(s, bt, max_i, max_j)]

    # bt TODO: No, I need all the alignments with positive scores.
    #max_s, max_o, max_t = max([ (s[l1-1, j], j, 0) for j in range(l2) ] ++ [ (s[i, l2-1], i, 1) for i in range(l1) ])
    print(valid_alignments)

    # next I need malign to clusters, or list of alignments to clusters.

def get_read_range(aligned_segment):
    ## get read_start and read_end(aligned coordinates in an original read with hard clipped part)

    #TODO: why sometime they don't have MD tag
    if aligned_segment.cigartuples:
        first_cigar = aligned_segment.cigartuples[0]
    else:
        first_cigar = (0, None)
    #first_cigar = aligned_segment.cigartuples[0]

    if (first_cigar[0] == 5):
        # correction for hard clipping: 5 for H
        return (aligned_segment.query_alignment_start + first_cigar[1],
                aligned_segment.query_alignment_end   + first_cigar[1])
    else:
        return (aligned_segment.query_alignment_start, aligned_segment.query_alignment_end)

# this parses sam records for a read # NOTE: deprecated, not used in encode_dp
def parse_sam_records(aln, mons):
    """
    calc encoded representation of the read.
    might be OK to ignore <11bp gap """

    epos_last_monomer = 0

    for mon in sorted([ (get_read_range(m), m) for m in mons if (not m.is_secondary) and m.has_tag("MD") ], key = lambda x: x[0]):
        rs, re = mon[0]
        m = mon[1] 
        rs_full = rs - m.reference_start
        re_full = re + aln.lengths[m.reference_id] - m.reference_end

        # Insert a GAP if there's (10) bp # NOTE: some ad-hoc workaround
        if rs_full - epos_last_monomer > 10:
            yield AssignedMonomer(begin = epos_last_monomer, end = rs_full,
                    monomer = Monomer(name = "GAP", snvs = []))

        if epos_last_monomer - rs_full > 0.8 * aln.lengths[m.reference_id]:
            # skip if there is 80% overlap with previous monomer; and epos_last_monomer need not changed # NOTE: ad hoc, too
            re_full = epos_last_monomer
        else:
            qseq, rseq = m.query_alignment_sequence, m.get_reference_sequence()
            qs, rs = m.query_alignment_start, m.reference_start
            # snvs = [ SNV(pos = j, base = qseq[i-qs]) for i,j in m.get_aligned_pairs(matches_only = True) if rseq[j-rs].islower() ] # This won't work 
            snvs = [ SNV(pos = j, base = qseq[i-qs]) for i,j in m.get_aligned_pairs(matches_only = True) if qseq[i-qs].upper() != rseq[j-rs].upper() ] 
            monomer = Monomer(name = aln.references[m.reference_id], snvs = snvs)
            yield AssignedMonomer(begin = rs_full, end = re_full, monomer = monomer)

        epos_last_monomer = re_full


# TODO: handle orientation info properly
def encodeSAM_DP(samfile, usemons = None):
    """ Encode PacBio reads in SAM format (reads as query, monomers as reference)
        NOTE: this tries to select best assignment using DP
        The result is written out in `out`.
    """

    # TODO: give up this. check if it's OK to set readlen to be 0
    def readname2len(s):
        #import re
        #a = re.sub(".*/", "", s).split("_")
        #return int(a[1]) - int(a[0])
        # NOTE: for ONT
        return 0

    aln = pysam.AlignmentFile(samfile)

    usemon_ids = [ i for i, n in enumerate(aln.references) if n in usemons ] if usemons else []

    for read, mons in groupby(aln.fetch(until_eof=True), key=lambda x: x.query_name):
        #yield EncodedRead(name = read, mons = list(parse_sam_records(aln, mons)), length = readname2len(read))
        length = readname2len(read)
        mons_f = []

        #if usemons:
        #    mons = [ m for m in mons if aln.references[m.reference_id] in usemons ]

        for mon in [ m for m in mons if (not m.is_secondary) and m.has_tag("MD") and m.get_tag("BS") > 50 ]:

            if usemon_ids and mon.reference_id not in usemon_ids:
                continue

            rs, re = get_read_range(mon)
            rs_full, re_full = rs - mon.reference_start, re + aln.lengths[mon.reference_id] - mon.reference_end
            if mon.is_reverse:
                re_full, rs_full = length - rs_full, length - re_full
            score = mon.get_tag("BS")
            mons_f += [(mon, rs_full, re_full, score, not mon.is_reverse)]

        if not mons_f:
            continue

        # NOTE: here is stringent heuristic; taking ones with maximal score
        mons_f = sorted(mons_f, key = lambda x: x[1])
        mons_f = groupby(mons_f, key = lambda x: x[1])
        mons_f = [ sorted(list(sg[1]), key = lambda x: -x[3])[0] for sg in mons_f ]

        print(read, flush=True)
        dp_s = np.zeros(len(mons_f))
        dp_t = [ () for i in range(len(mons_f)) ] # empty traceback

        if False: # NOTE: pure DP. BS threshold 50. overlap penalty 2 / bp
            for i in range(len(mons_f)):
                cand = [ (dp_s[j] + (mons_f[i][3] - 50) +  min([0, mons_f[i][1] - mons_f[j][2]]) * 2, j) for j in range(i) ] + [(mons_f[i][3] - 50, -1)]
                # TODO: I may need suboptimal ones too...
                # TODO: these 3 lines can be one?
                max_s, max_j = max(cand)
                dp_s[i] = max_s
                dp_t[i] = max_j

        # NOTE: optimized impl.
        max_past = -1 # index
        cand = [] # overlapping extension is possible only from these.
        for i in range(len(mons_f)):
            if max_past != -1:
                max_past = max( [ (dp_s[j], j) for j in cand if mons_f[j][2] <= mons_f[i][1] ] + [(dp_s[max_past], max_past)] )[1]
            else:
                max_past = max( [ (dp_s[j], j) for j in cand if mons_f[j][2] <= mons_f[i][1] ] + [(-1, -1)] )[1]
            cand = [ j for j in cand if mons_f[j][2] > mons_f[i][1] ]
            max_s, max_j = dp_s[max_past] + (mons_f[i][3] - 50), max_past
            for j in cand:
                s = dp_s[j] + (mons_f[i][3] - 50) + min([0, mons_f[i][1] - mons_f[j][2]]) * 2
                if s > max_s:
                    max_s, max_j = s, j
            dp_s[i], dp_t[i] = max_s, max_j
            cand += [i]

        # traceback
        result = []
        last_i = max([ (dp_s[i], i) for i in range(len(mons_f)) ])[1]
        while last_i != -1:
            result += [mons_f[last_i]]
            last_i = dp_t[last_i]

        # TODO: put pickled file instead of printing
        #len([ m for m in result if m[3] ]) / len(result) >= 0.5

        result = sorted(result, key=lambda x:x[1])

        def r_to_m(r_cmp):
            r = r_cmp[0]
            qseq, rseq = r.query_alignment_sequence, r.get_reference_sequence()
            qs, rs = r.query_alignment_start, r.reference_start
            snvs = [ SNV(pos = j, base = qseq[i-qs]) for i,j in r.get_aligned_pairs(matches_only = True) if qseq[i-qs].upper() != rseq[j-rs].upper() ] 
            monomer = Monomer(name = aln.references[r.reference_id], snvs = snvs)
            return AssignedMonomer(begin = r_cmp[1], end = r_cmp[2], monomer = monomer, ori = '+' if r_cmp[4] else '-')

        #yield EncodedRead(name = read, mons = [ r_to_m(r) for r in result ], length = readname2len(read))
        yield EncodedRead(name = read, mons = [ r_to_m(r) for r in result ], length = length)

def correct(reads, variants):
    def correct_read(read, variants):
        return EncodedRead(
               name = read.name,
               length = read.length,
               mons = [ AssignedMonomer(
                   begin = m.begin,
                   end = m.end,
                   monomer = Monomer(
                       name = m.monomer.name,
                       snvs = [ snv for snv in m.monomer.snvs if variants[m.monomer.name][(snv.pos, snv.base)] > 0 ]))
                   for m in read.mons ]) 

    return [ correct_read(read, variants) for read in reads ]


def encoding_stats(ers, variants = None):
    """ print out some statistics of encoded reads set to see the quality of encoding. """

    print(f"name\tcoverage\tsnv_rate(on read)\tsnv_rate(on SNV sites)\tnSNV\tnsites\tmonomer_length(after recover)")
    c_total, s_total, l_total, nread= 0, 0, 0, 0
    v_total = 0
    for ER in ers:
        c, s, v = 0, 0, 0
        for m in ER.mons:
            #print(f"{m.begin}-{m.end} : {m.monomer.name} ({len(m.monomer.snvs)} SNVs)")
            if m.monomer.name != "GAP":
                c += m.end - m.begin
                if variants:
                    # v += len(list(groupby(m.monomer.snvs, lambda x: x.pos)))
                    v += len(list(groupby([ k[0] for k, v in variants[m.monomer.name].items() if v > 0 ])))
            s += len(m.monomer.snvs)
        if c > 0:
            if variants:
                print(f"{ER.name}\t{100*c/ER.length:.3f} %\t{100*s/c:.3f} %\t{100*s/v:.3f} %\t{s}\t{v}\t{c}")
            else:
                print(f"{ER.name}\t{100*c/ER.length:.3f} %\t{100*s/c:.3f} %\tNA\t{s}\tNA\t{c}")
        else:
            print(f"{ER.name}\t0 %\t0 %")
        c_total += c
        s_total += s
        l_total += ER.length
        nread += 1
        v_total += v

    if variants:
        print(f"TOTAL\t{100*c_total/l_total:.3f} %\t{100*s_total/c_total:.3f} %\t{100*s_total/v_total:.3f}\t{s_total}\t{v_total}\t{c_total}")
    else:
        print(f"TOTAL\t{100*c_total/l_total:.3f} %\t{100*s_total/c_total:.3f} %\tNA\t{s_total}\tNA\t{c_total}")

    print(f"statistics of {nread} reads.")

def kmonomers_stats(ers):
    """ Calculates k-monomer frequencies in given encoded read set. """

    for k in [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]:
        counter = Counter()
        for er in ers:
            # mons_array = [ m.monomer.name for m in er.mons if (m.monomer.name != "GAP") | ((m.monomer.end - m.monomer.begin) > 100) ]
            mons_array = [ m.monomer.name for m in er.mons if m.monomer.name != "GAP" ]
            for i in range(len(mons_array)-k+1):
                counter[tuple(mons_array[i:i+k])] += 1

        print(f"***** K = {k} *****")
        res = [ c for c in counter.most_common(500) if c[1] > 0 ]
        print(*res, sep="\n")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Monomer-encoding of reads.')
    parser.add_argument('action', metavar='action', type=str, help='action to perform: encode_dp, encode, correct, k-monomer, ...')
    parser.add_argument('--sam', dest='samfile', help='SAM format alignments to be encoded')
    parser.add_argument('--split', dest='split', help='specify to split the encoded reads pickle (for larger sam file)')
    parser.add_argument('--read', dest='readfile', help='pickled encoded reads')
    parser.add_argument('--vars', dest='varsfile', help='pickled admissible variants')
    parser.add_argument('--mons', dest='monsfile', help='list of monomers to be used; default to use all.')
    parser.add_argument('--out', dest='outfile', help='the output to which pickled encoded reads are written out.' +
            ' For correction task, corrected reads are written out to the path')
    args = parser.parse_args()

    print(args.action)

    if args.action == "encode_dp":
        # This tries to select good assignments based on blast score and overlapping among them
        assert args.samfile, "SAM file is not specified. aborting."
        assert args.outfile, "output file is not specified. aborting."

        if args.monsfile:
            usemons = [ l.strip() for l in open(args.monsfile, "r").readlines() ]
        else:
            usemons = None

        ers = encodeSAM_DP(args.samfile, usemons)
        with open(args.outfile, "wb") as f:
            pickle.dump(list(ers), f)

    elif args.action == "correct":
        assert args.readfile, "file of encoded reads is not specified. aborting."
        with open(args.readfile, "rb") as f:
            ers = pickle.load(f)

        assert args.varsfile, "file of variants is not specified. aborting."
        with open(args.varsfile, "rb") as f:
            variants = pickle.load(f)

        encoding_stats(ers)
        corrected = correct(ers, variants)
        encoding_stats(corrected, variants)
        
        with open(args.outfile, "wb") as f:
            pickle.dump(corrected, f)

    ## NOTE: replaced by `kmer` in HOR_segregation.py ??
    elif args.action == "k-monomer":
        assert args.readfile, "file of encoded reads is not specified. aborting."
        with open(args.readfile, "rb") as f:
            ers = pickle.load(f)
        kmonomers_stats(ers)

    else:
        print(f"unknown action. {args.action}")
