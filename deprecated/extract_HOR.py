## This tries to extract HOR structure for each read, which might be used in (roughly) identifying overlap candidate of the reads 
import pysam
import itertools
import numpy as np

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

        if epos_last_monomer - rs_full > 0.8 * aln.lengths[m.reference_id]:
            # skip if there is 80% overlap with previous monomer; and epos_last_monomer need not changed
            re_full = epos_last_monomer
            pass
        else:
            yield (aln.references[m.reference_id], rs_full, re_full)
        epos_last_monomer = re_full

aln = pysam.AlignmentFile("./show-reads/md_k7_s25_m015_143cells_primary_hd1khd")
rit = itertools.groupby(aln.fetch(until_eof=True), key=lambda x: x.query_name)
x = [(r, list(get_encoded_read(aln, ms))) for r, ms in rit]
distmat = np.load("./count_vars/d0.distmat.npy")
name_to_idx = { x:i for (i, x) in enumerate(aln.references) }
print(len(name_to_idx))
## TODO: nice measures for gap by extending distmat
name_to_idx["GAP"] = 0

def power_spec(ms_gap):
    # Does ms have periodicity of p ? a few option to check this.
    # 1. calc a kind of power spectrum

    ms = list(filter(lambda x: (x[0] != "GAP") | (x[2] - x[1] > 100), ms_gap))
    print(ms_gap)
    print(ms)
    for p in range(1,30):
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
            fmtstr = "["
            for ll in l:
                if ll == 0:
                    fmtstr += "\033[44m 0\033[0m, "
                elif ll < 10:
                    fmtstr += f"\033[46m {ll}\033[0m, "
                elif ll < 20:
                    fmtstr += f"\033[47m{ll}\033[0m, "
                else:
                    fmtstr += f"{ll}, "

            #print(f"{l}\n")
            print(f"{fmtstr[:-2]}]\n")

def self_off_diag(ms):


    # 2. self-wrap dp ???
    def match(m1, m2):
        I = 0.1
        # match score for each case
        if (m1[0] == "GAP") & (m2[0] == "GAP"):
            return I * abs( (m1[2]-m1[1]) - (m2[2]-m2[1]) )
        elif m1[0] == "GAP":
            return np.mean(distmat, 0)[name_to_idx[m2[0]]] + I * abs( (m1[2]-m1[1]) - (m2[2]-m2[1]) )
        elif m2[0] == "GAP":
            return np.mean(distmat, 0)[name_to_idx[m1[0]]] + I * abs( (m1[2]-m1[1]) - (m2[2]-m2[1]) )
        else:
            return distmat[name_to_idx[m1[0]], name_to_idx[m2[0]]]

    def skip(m):
        # penalty for skipping this
        I = 0.1 # TODO: correlate to I above; do we have minus sized gap?
        return I * abs(m[2] - m[1])
        

    ms = list(ms)
    s = np.zeros(len(ms) * len(ms)).reshape(len(ms), len(ms))
    trace = np.zeros(len(ms) * len(ms)).reshape(len(ms), len(ms))

    # only part i < j is required. 
    # init
    s[0,0] = float("inf")
    for j in range(1, len(ms)):
        # s[i,0] = match(ms[i], ms[0])
        s[0,j] = match(ms[0], ms[j])
        s[j,j] = float("inf")

    # recur
    for i in range(1, len(ms)):
        for j in range(i+1, len(ms)):
            # match
            s_match = s[i-1, j-1] + match(ms[i], ms[j])
            # skip in i
            s_skipi = s[i-1, j] + skip(ms[i])
            # skip in j
            s_skipj = s[i, j-1] + skip(ms[j])
            # select the best
            s[i, j], trace[i, j] =  min([(s_match, 0), (s_skipi, 1), (s_skipj, 2)])

    # trace
    # from j==max
    for i in range(2, 30):
        pass
        # from s[len(ms)-i ,len(ms)-1]
            
    # report: 
    np.set_printoptions(threshold=40000)
    print("DP table is")
    print(s)

    mtab = np.zeros(len(ms) * len(ms)).reshape(len(ms), len(ms))
    for i in range(len(ms)):
        for j in range(i+1, len(ms)):
            mtab[i,j] = match(ms[i], ms[j])
    print("mtab is")
    print(mtab)

for r, ms_gap in x:
    #ms_sorted = sorted(list(ms_gap), key=lambda x:get_read_range(x)[0])
    print(f"power_spec for {r}")
    power_spec(ms_gap)
    #print(f"self_off_diag for {r}")
    #self_off_diag(ms_gap)

## TODO: write menu
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Interact with monomer databases.')
    parser.add_argument('action', metavar='action', type=str, help='action to perform: distmat, ...')
    parser.add_argument('--mons', dest='monfile', help='path to monomers.fa')
    parser.add_argument('--out', dest='outfile', help='path to output')
    args = parser.parse_args()

    print(args.action)

    if args.action == "distmat":
        assert args.monfile, "monomers database is not specified. aborting."
        if args.outfile:
            save_pairwise_edit_distance(args.monfile, args.outfile)
        else:
            save_pairwise_edit_distance(args.monfile)
    else:
        print(f"unknown action. {args.action}")

