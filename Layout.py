# This is evolved from Alignments.py, which is then evolved from SNV_distribution.py

from collections import Counter
from collections import namedtuple
import numpy as np

from EncodedRead import *
from HOR_segregation import *

Pool = namedtuple("Pool", ("hers", "arrs")) # simple class for the context; hor encoded reads, layouts, ...
pool = None # TODO; how can I make sure this is singleton (like the pattern)? namedtuple is immutable!
# bits_dict should be in context...

# a record for an alignment.
#   eff_ovlp = effective length of ovlp, f_ext/r_ext = extention forward/reverse
Aln = namedtuple("Aln", ("i", "ai", "li", "j", "aj", "lj", "k", "score", "len_ovlp", "eff_ovlp", "f_ext", "r_ext"))

# 1 - First, Detect variants
def detect_snvs(units):
    """
    input: [(read id, mon id of the head of the HOR unit)]
    return: [(monomer index in unit, position in monomer, alt. base, relative frequency)]
    """

    n_units = len(units)
    counter = Counter()
    for ri, i in units:
        for j in range(2, 12): # I'm skipping the first 2 monomers, where assignment is always wrong
            counter.update([ (j, s.pos, s.base) for s in pool.hers[ri].mons[i+j].monomer.snvs ])
    # TODO: make cleaner
    return [ (k, p, b, c/n_units) for (k, p, b), c in counter.most_common(1000) if (0.05 < c/n_units) and (c/n_units < 0.8) ]

# 2 and 3 - filling up gaps, get continuous regions.
def valid_read_regions(hers):
    """
        get valid consecutive regions after filling up 22/23/24-mons or 33/34/35-mons gaps
        this is a version for hor encoded reads, instead of a bitvector matrix
        return { ri : [[mi, mi, mi, ...], [mi, mi, ...]] } (mi = -1 for masked)
    """

    regs = {}
    for i, her in enumerate(hers):
        regs[i], ai, lh = [], -1, -36
        for h, _, t in her.hors:
            if t != "~":
                continue
            if h == lh + 12:
                regs[i][ai] += [h]
            elif h in [lh + 22, lh + 23, lh + 24]:
                regs[i][ai] += [-1, h]
            elif h in [lh + 33, lh + 34, lh + 35]:
                regs[i][ai] += [-1, -1, h]
            else:
                regs[i] += [[h]]
                ai += 1
            lh = h
        regs[i] = [ l for l in regs[i] if l ] # remove if empty
    return { k : v for k, v in regs.items() if v } # remove if empty

def calc_align(i, ai, j, aj, k, bits_dict):

    li, lj = bits_dict[(i, ai)].shape[0], bits_dict[(j, aj)].shape[0]

    if k < 0:
        len_ovlp = min(li, k+lj)
        xi = bits_dict[(i, ai)][0:len_ovlp,:]
        xj = bits_dict[(j, aj)][-k:-k+len_ovlp,:]
    else:
        len_ovlp = min(lj, -k+li) 
        xi = bits_dict[(i, ai)][k:k+len_ovlp,:]
        xj = bits_dict[(j, aj)][0:len_ovlp,]

    m = np.multiply(xi, xj) # 1: match, 0: masked, -1: mismatch
    match = np.multiply((m + 1), m) / 2 # 1: match, 0: masked/mismatch
    mismatch = np.sum(np.multiply(m, m) - match) # 1: mismatch, 0: masked/match # n01 or n10
    n11 = np.sum(np.multiply(match, (xi + 1) / 2)) # n11
    score = int(1000 * n11 / (1 + n11 + mismatch)) # per-mille !
    eff_ovlp = len_ovlp - sum([ 1 if x == 0 else 0 for x in m[:,0] ])

    return Aln(i = i, ai = ai, li = li, j = j, aj = aj, lj = lj, k = k,
        score = score, len_ovlp = len_ovlp, eff_ovlp = eff_ovlp, f_ext = max(0, k+lj-li), r_ext = max(0, -k))

def print_align(aln, bits_dict):

    pair_to_char = {
            ( 1, 1): "#", ( 1, 0): "0", ( 1, -1): ".",
            ( 0, 1): "0", ( 0, 0): "@", ( 0, -1): "0",
            (-1, 1): ",", (-1, 0): "0", (-1, -1): "_"}

    bit_to_char = { 1: "x", 0: "o", -1: " "}

    r, q = bits_dict[(aln.i, aln.ai)], bits_dict[(aln.j, aln.aj)]
    rs, qs = max(0, aln.k), max(-aln.k, 0)
    re, qe = rs + aln.len_ovlp, qs + aln.len_ovlp

    #if is_top:
    print(f"\n*\tri\tqi\trs\tre\trl\tqs\tqe\tql\tidt.\tovlp\teff.ovlp")
    print(f"ALIGN\t{aln.i}\t{aln.j}\t{rs}\t{re}\t{aln.li}\t{qs}\t{qe}\t{aln.lj}\t{aln.score/10:.1f}\t{aln.len_ovlp}\t{aln.eff_ovlp}\n")

    _l = f"== {r.shape[1]} SNVs =="
    lines = f"  i\t{_l:^40}\tj  "

    # dangling part if any
    for i in range(rs):
        lines += f"\n{i:>3}\t" + "".join([ bit_to_char[e] for e in r[i,:] ]) + f"\t***"
    for i in range(qs):
        lines += f"\n***\t" + "".join([ bit_to_char[e] for e in q[i,:] ]) + f"\t{i:<3}"

    # the aligned part
    for i in range(re-rs):
        lines += f"\n{rs+i:>3}\t" + "".join([ pair_to_char[(e, d)] for e, d in zip(r[rs+i,:], q[qs+i,:]) ]) + f"\t{qs+i:<3}"

    # the remaining dangling part if any
    for i in range(re, aln.li):
        lines += f"\n{i:>3}\t" + "".join([ bit_to_char[e] for e in r[i,:] ]) + f"\t***"

    for i in range(qe, aln.lj):
        lines += f"\n***\t" + "".join([ bit_to_char[e] for e in q[i,:] ]) + f"\t{i:<3}"

    print(lines)


def layout(alns_pers):
    """ implementing layout idea; alns_pers = initial pairwise alignment, for now. """

    # setup contexts
    bag_of_units = [ (ri, h) for ri, er in enumerate(pool.hers) for h, _, t in er.hors if t == "~" ]
    n_units = len(bag_of_units)
    print(f"{n_units} units found in {len(pool.hers)} reads.")
    snv_sites = detect_snvs(bag_of_units)
    print(f"{len(snv_sites)} SNV sites defined.")
    arrs = valid_read_regions(pool.hers)
    print(f"{len(arrs)} read has an array to be aligned.")
    print(f"{ sum([ len(v) for k, v in arrs.items() ]) } arrays in total.")

    # construct bit vectors in dict
    bits_dict = {}
    for i, ai, l in [ (i, ai, l) for i, a in arrs.items() for ai, l in enumerate(a) if len(l) > 4 ]:
        v = []
        for h in l:
            if h == -1:
                v += [0] * len(snv_sites)
            else:
                # NOTE: Was snv a dict ? # TODO: cleanup
                v += [ 1 if any([ (sp.pos, sp.base) == (p, b) for sp in pool.hers[i].mons[h+k].monomer.snvs ]) else -1 for k, p, b, f in snv_sites ]
        bits_dict[(i, ai)] = np.array(v).reshape(len(l), len(snv_sites))

    # calculate alignments globally
    long_arrays = [ (i, ai) for i, a in arrs.items() for ai, l in enumerate(a) if len(l) > 9 ]
    alns_dict = dict()
    n = 0
    for i, ai in long_arrays:
        n += 1
        print(f"aligning {i} - {ai}. {n} / {len(long_arrays)}")
        alns_dict[(i, ai)] = dict()
        for j, aj in long_arrays:
            if (i, ai) == (j, aj):
                continue
            li, lj = bits_dict[(i, ai)].shape[0], bits_dict[(j, aj)].shape[0]
            # with at least 2 units in overlap
            alns_dict[(i, ai)][(j, aj)] = sorted([ calc_align(i, ai, j, aj, k, bits_dict = bits_dict) for k in range(-lj+2, li-2) ], key = lambda x: -1 * x.score)

        print(f"{ len([ (j, aj, aln) for (j, aj), alns in alns_dict[(i, ai)].items() for aln in alns if aln.score > 500 and aln.eff_ovlp > 4 ]) } targets found.")

    import pickle
    pickle.dump(alns_dict, open(f"alns_dict.10u.pickle", "wb"))

    import sys
    sys.exit()

    goods = [ (j, aj, aln) for (j, aj), alns in alns_dict[long_arrays[0]].items() for aln in alns if aln.score > 600 ]

    print(len(goods))
    print(goods)
    print(goods[0])

    for i in range(5):
        print_align(goods[i][2], bits_dict)

    import sys
    sys.exit()
    
    pw_aln = pickle.load(open(alns_pers, "rb"))
    print(f"{ sum([ len(pw_aln[i]) for i in pw_aln.keys() ]) } alignments for {len(pw_aln)} reads loaded.")

    def Fx(x):
        """ X to S
            scoring scheme: n11 / (n10+n01+n11) """
        #assert bv1.shape == bv2.shape, "wrong shapes"
        #m = np.multiply(bv1, bv2) # 1: match, 0: masked, -1: mismatch
        #mt = np.multiply((m + 1), m) / 2 # 1: match, 0: masked/mismatch
        #n_m = np.sum(np.multiply(m, m) - mt) # n01 or n10
        #n11 = np.sum(np.multiply(mt, (bv1 + 1) / 2)) # n11
        #return int(1000 * n11 / (1 + n11 + n_m)) # per-mille !
        return 0

    def aln_x(i, j, k):
        # TODO: memoize ?
        li, lj = [ bms[x].shape[0] for x in [i, j] ]
    
        if k < 0:
            x1, x2 = bms[i][0:min([li, k+lj]),3:], bms[j][-k:min([lj, -k+li]),3:]
        else:
            x1, x2 = bms[i][k:min([li, k+lj]),3:], bms[j][ 0:min([lj, -k+li]),3:]

        m = np.multiply(x1, x2) # 1: match, 0: masked, -1: mismatch
        #mt = np.multiply((m + 1), m) / 2 # 1: match, 0: masked/mismatch
        #n_m = np.sum(np.multiply(m, m) - mt) # n01 or n10
        #n11 = np.sum(np.multiply(mt, (bv1 + 1) / 2)) # n11
        #return int(1000 * n11 / (1 + n11 + n_m)) # per-mille !

        #return Fx(x), aln_len
        return 0

    layout = []

    def get_slippy():
        """ temporarily I'm using context """

        def flatten(l):
            from itertools import chain
            return list(chain.from_iterable(l))

        def tr_ext(r): # one-step transitive extension # OBVIOUSLY TOO SLOW
            # reflexive and symmetric
            x = [ [ (i, k) for j2, k in r if j == j2 ] for i, j in r ]
            r = { (i, j) for i, j in flatten(x) }
            r |= { (j, i) for i, j in r } | { (i, i) for i, j in r } | { (j, j) for i, j in r }
            return r

        r0 = { (i, j) for i, v in pw_aln.items() for j, alns in v.items() if len([0 for a in alns if a[4] > 500 and a[5] > 4]) > 1 }



        n_r0 = -1
        while n_r0 < len(r0):
            n_r0 = len(r0)
            print(len(r0))
            r0 = tr_ext(r0)

    get_slippy()

    #r1 = { (i, j) for i in Set(r0.flatten(

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='perform pair-wise alignment among HOR encoded reads on SNV data, again') # TODO: explain
    parser.add_argument('action', metavar='action', type=str, help='action to perform: align, ...')
    parser.add_argument('--hor-reads', dest='hors', help='pickled hor-encoded long reads')
    parser.add_argument('--alns', dest='alns', help='pickled pickled alignments')
    args = parser.parse_args()

    if args.action == "align":
        assert args.hors, "need HOR-encoded reads"
        hers = pickle.load(open(args.hors, "rb"))
        if args.alns:
            pass
            #align(args.alns)
        #else:
            #align()

    if args.action == "layout":
        assert args.hors, "need HOR-encoded reads"
        hers = pickle.load(open(args.hors, "rb"))
        pool = Pool(hers = hers, arrs = None)
        if args.alns:
            layout(args.alns)
        #else:
        #    layout()
    else:
        assert None, "invalid action."
