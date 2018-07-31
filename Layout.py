# This is evolved from Alignments.py, which is then evolved from SNV_distribution.py
from collections import Counter
from collections import namedtuple
import numpy as np
import pickle

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

def get_bits_dict(snvs, regions):
    """ construct bit vectors in dict """
    bits_dict = {}
    for i, ai, l in [ (i, ai, pool.arrs[i][ai]) for i, ai in regions ]:
        v = []
        for h in l:
            if h == -1:
                v += [0] * len(snvs)
            else:
                # NOTE: Was snv a dict ? # TODO: cleanup
                v += [ 1 if any([ (sp.pos, sp.base) == (p, b) for sp in pool.hers[i].mons[h+k].monomer.snvs ]) else -1 for k, p, b, f in snvs ]
        bits_dict[(i, ai)] = np.array(v).reshape(len(l), len(snvs))
    return bits_dict

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

def double_edge_component(alns_dict):
    """
    obtain so called slippy component, whose elements are connected by multi-edges.
    edge filter is hard-coded to be at least 5 units with 50% identity. TODO: see the distribution to support this!
    """

    def is_double(alns):
        return 1 < len([ aln for aln in alns if aln.score > 500 and aln.len_ovlp > 4 ])

    nodes = set()
    for (i, ai), d in alns_dict.items():
        nodes |= { (i, ai) } | { (j, aj) for (j, aj), alns in d.items() }
    print(f"{len(nodes)} nodes in total.")

    components = []
    while len(nodes) > 0:
        k, ak = nodes.pop()
        print(f"\nchecking component for {k}-{ak} / {len(nodes)} remain to be checked.")

        curr_component = { (k, ak) } | { (j, aj) for (j, aj), alns in alns_dict[(k, ak)].items() if is_double(alns) }
        curr_checked = set()

        while len(curr_checked) < len(curr_component):
            l, al = (curr_component - curr_checked).pop()
            curr_checked |= { (l, al) }
            curr_component |= { (j, aj) for (j, aj), alns in alns_dict[(l, al)].items() if is_double(alns) }

        components += [[ e for e in curr_component ]]
        nodes -= curr_component
        print(f"{len(curr_component)} nodes in the component")

    return components

def layout(alns_dict = None):
    """ implementing layout idea """

    # setup contexts
    bag_of_units = [ (ri, h) for ri, er in enumerate(pool.hers) for h, _, t in er.hors if t == "~" ]
    n_units = len(bag_of_units)
    print(f"{n_units} units found in {len(pool.hers)} reads.")

    snv_sites = detect_snvs(bag_of_units)
    print(f"{len(snv_sites)} SNV sites defined globally.")

    arrs = valid_read_regions(pool.hers)
    print(f"{len(arrs)} read has an array to be aligned.")
    print(f"{ sum([ len(v) for k, v in arrs.items() ]) } in total.")

    pool = Pool(hers = hers, arrs = arrs) # NOTE: I wanted to express global variables ...

    # bits vectors are constructed for all regions, on SNVs defined globally.
    regs = [ (i, ai) for i, a in pool.arrs.items() for ai, l in enumerate(a) if len(l) > 4 ]
    bits_dict = get_bits_dict(snv_sites, regs)

    # calculate alignments globally
    long_arrays = [ (i, ai) for i, a in arrs.items() for ai, l in enumerate(a) if len(l) > 9 ]

    # TODO: once again abstract the logic to calc all-vs-all alns_dict for subset of reads (i mean, slippy parts).
    if not alns_dict:
        alns_dict = dict()
        n = 0
        for i, ai in long_arrays:
            n += 1
            print(f"aligning {i} - {ai}. {n} / {len(long_arrays)}")
            alns_dict[(i, ai)] = dict()
            for j, aj in [ (j, aj) for j, aj in long_arrays if not (j, aj) == (i, ai)]:
                li, lj = bits_dict[(i, ai)].shape[0], bits_dict[(j, aj)].shape[0]
                # with at least 2 units in overlap
                alns = [ calc_align(i, ai, j, aj, k, bits_dict = bits_dict) for k in range(-lj+2, li-2) ]
                alns = [ aln for aln in alns if aln.score > 400 and aln.eff_ovlp > 2 ]
                if alns:
                    alns_dict[(i, ai)][(j, aj)] = sorted(alns, key = lambda x: -1 * x.score)
            print(f"{ len([ (j, aj, aln) for (j, aj), alns in alns_dict[(i, ai)].items() for aln in alns if aln.score > 400 and aln.eff_ovlp > 2 ]) } targets found.")

        pickle.dump(alns_dict, open(f"alns_dict.10u.400-2.pickle", "wb"))

    slippies = double_edge_component(alns_dict)
    print(f"{ len(slippies) } components found.")
    print(slippies)

    for ic, comp in enumerate(slippies):

        if len(comp) < 2:
            continue
        units_in_comp = [ (i, mi) for i, ai in comp for mi in arrs[i][ai] if mi > -1 ]

        snvs = detect_snvs(units_in_comp)
        print(f"\n{len(snvs)} SNVs for the Component {ic}: {len(comp)} reads; {len(units_in_comp)} units.")

        for k, p, b, f in sorted(snvs):
            print(f"{k}\t{p}\t{b}\t{f:.2f}")

    import sys
    sys.exit()

    # TODO: I need to define layout, alignment between layout/read vs layout/read, check for internal consistency, ... then final iterative assembly procedure.
    #layout = []

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
        pool = Pool(hers = hers, arrs = None) # NOTE: I wanted to express global variables ...
        if args.alns:
            alns_dict = pickle.load(open(args.alns, "rb"))
            layout(alns_dict)
        else:
            layout()
    else:
        assert None, "invalid action."
