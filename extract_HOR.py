import pysam
import itertools

# parsed = [(r, list(ms)) for r, ms in rit]

def get_read_range(aligned_segment):
    ## get read_start and read_end
    ## (aligned coordinates in an original read including hard clipped part)
    a = aligned_segment
    read_start, read_end = a.query_alignment_start, a.query_alignment_end
    # correction for hard clipping: 5 for H
    if (a.cigartuples[0][0] == 5):
        read_start += a.cigartuples[0][1]
        read_end   += a.cigartuples[0][1]
    return (read_start, read_end)

def get_encoded_read(aln, ms):
    """
    calc encoded representation of the read.
    might be OK to ignore <5bp gap """
    epos_last_monomer = 0
    for m in sorted(list(ms), key=lambda x:get_read_range(x)[0]):
        rs, re = get_read_range(m)
        rs_full = rs - m.reference_start
        re_full = re + aln.lengths[m.reference_id] - m.reference_end
        if rs_full - epos_last_monomer > 0:
            yield ("GAP", epos_last_monomer, rs_full)
        yield (aln.references[m.reference_id], rs_full, re_full)
        epos_last_monomer = re_full

aln = pysam.AlignmentFile("./show-reads/md_k8_s25_m015_143cells_primary_hd1khd")
rit = itertools.groupby(aln.fetch(until_eof=True), key=lambda x: x.query_name)
x = [(r, list(get_encoded_read(aln, ms))) for r, ms in rit]

import numpy as np

distmat = np.load("./count_vars/d0.distmat.npy")

name_to_idx = { x:i for (i, x) in enumerate(aln.references) }
print(len(name_to_idx))

## TODO: nice measures for gap by extending distmat
name_to_idx["GAP"] = 0

for r, ms_gap in x[:10]:

    ms = list(filter(lambda x: (x[0] != "GAP") | (x[2] - x[1] > 100), ms_gap))
    print(ms_gap)
    print(ms)

    for p in range(1,30):
        # Does ms have periodicity of p ? a few option to check this.
        # 1. calc a kind of power spectrum
        # l = [ distmat[name_to_idx[ms[i][0]], name_to_idx[ms[i+p][0]]] for i in range(len(ms)) if i+p < len(ms) ]
        l = []
        for i in range(len(ms)):
            if i+p < len(ms):
                if (ms[i][0] == ms[i+p][0]) & (ms[i][0] == "GAP"):
                    #l += [-0.1]
                    pass
                elif (ms[i][0] == "GAP") | (ms[i+p][0] == "GAP"):
                    #l += [10.1]
                    pass
                else:
                    l += [ distmat[name_to_idx[ms[i][0]], name_to_idx[ms[i+p][0]]] ]

        if l:
            print(f"\t\t{p} : {sum(l) / len(l)}")
            print(f"{l}\n")

        # 2. self-wrap dp ???
            
