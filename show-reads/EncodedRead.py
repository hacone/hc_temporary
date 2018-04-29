## This is the class representing a set of monomer-encoded PacBio reads.
from collections import namedtuple
from itertools import groupby
import numpy as np
import pysam
import re
import sys

# These are all immutable
SNV = namedtuple("SNV", ("pos", "base")) # Int, char
Monomer = namedtuple("Monomer", ("name", "snvs")) # string, list(SNV) or id: Int, list(SNV)
AssignedMonomer = namedtuple("AssignedMonomer", ("begin", "end", "monomer")) # Int, Int, Monomer
EncodedRead = namedtuple("EncodedRead", ("name", "mons", "length")) # string, list(AssignedMonomer), Int

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

# this parses sam records for a read
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

        # Insert a GAP if there's (10) bp
        if rs_full - epos_last_monomer > 10:
            yield AssignedMonomer(begin = epos_last_monomer, end = rs_full,
                    monomer = Monomer(name = "GAP", snvs = []))

        if epos_last_monomer - rs_full > 0.8 * aln.lengths[m.reference_id]:
            # skip if there is 80% overlap with previous monomer; and epos_last_monomer need not changed
            re_full = epos_last_monomer
        else:
            qseq, rseq = m.query_alignment_sequence, m.get_reference_sequence()
            qs, rs = m.query_alignment_start, m.reference_start
            # snvs = [ SNV(pos = j, base = qseq[i-qs]) for i,j in m.get_aligned_pairs(matches_only = True) if rseq[j-rs].islower() ] # This won't work 
            snvs = [ SNV(pos = j, base = qseq[i-qs]) for i,j in m.get_aligned_pairs(matches_only = True) if qseq[i-qs].upper() != rseq[j-rs].upper() ] 
            monomer = Monomer(name = aln.references[m.reference_id], snvs = snvs)
            yield AssignedMonomer(begin = rs_full, end = re_full, monomer = monomer)

        epos_last_monomer = re_full


def readname2len(s):
    a = re.sub(".*/", "", s).split("_")
    return int(a[1]) - int(a[0])

def fromSAM(sam):
    aln = pysam.AlignmentFile(sam)
    return [ EncodedRead(name = read, mons = list(parse_sam_records(aln, mons)), length = readname2len(read))
            for read, mons in groupby(aln.fetch(until_eof=True), key=lambda x: x.query_name) ]

if __name__ == "__main__":

    if len(sys.argv) > 1:
        ERs = fromSAM(sys.argv[1])
    else:
        ERs = fromSAM()

    print("readname\t%covered\t%errors")

    c_total, s_total, l_total, nread= 0, 0, 0, 0
    for ER in ERs:
        c, s = 0, 0
        for m in ER.mons:
            #print(f"{m.begin}-{m.end} : {m.monomer.name} ({len(m.monomer.snvs)} SNVs)")
            if m.monomer.name != "GAP":
                c += m.end - m.begin
            s += len(m.monomer.snvs)

        if c > 0:
            print(f"{ER.name}\t{100*c/ER.length:.3f} %\t{100*s/c:.3f} %")
        else:
            print(f"{ER.name}\t0 %\t0 %")

        c_total += c
        s_total += s
        l_total += ER.length
        nread += 1

    print(f"TOTAL\t{100*c_total/l_total:.3f} %\t{100*s_total/c_total:.3f} %\t{nread} reads")

