# Utilities to analyze HOR encoded reads before alignment.
# As a standalone script, this outputs relevant data and figures available at this stage.

from scipy.stats import binom
from collections import Counter
import hashlib
from collections import namedtuple
import numpy as np
import pickle
# from underscore import _ as us

# TODO: I just need datatype definition. That can be separated from other codes.
from EncodedRead import *
from HOR_segregation import *

# NOTE: This new definition is not compatible with one in Alignment.py
Aln = namedtuple("Aln", ("i", "ai", "li", "j", "aj", "lj", "k", "lov", "eov", "f_ext", "r_ext", "n11", "n00", "nmis", "score"))

# TODO: I'll work on this.
# TODO: can this be abstracted to work with other chromosomes?

def fillx(read, verbose = False):

    """
    this is more like `fill`ing appropriate number of `*` HOR units.
    get an array from a read. returns (h, s, t) : midx, size, type
    """

    hors_found = [ (h, s, t, h+s-1) for h, s, t in read.hors if t in ["~", "D39", "D28", "22U", "D1", "D12"] ]
    reg = []

    for i in range(len(hors_found) - 1):
        h, s, t, e = hors_found[i]
        nh, ns, nt, ne = hors_found[i+1]

        if read.ori == '+':
            gap = read.mons[nh].begin - read.mons[e].end
        if read.ori == '-':
            gap = read.mons[e].begin - read.mons[nh].end

        if verbose:
            print(f"{hors_found[i]}\t{hors_found[i+1]}\t{gap}\t{read.ori}")

        reg += [(h, s, t)]
        for j in range(int((gap + 1000) / 2057)):
            reg += [(h, s, "*")]

    # add the last HOR
    if len(hors_found) > 0:
        h, s, t, e = hors_found[len(hors_found)-1]
        reg += [(h, s, t)]

    if verbose and reg:
        print(reg)
        print(f"reglen\t{len(reg)}")

    return reg

def var(reads, units = None, hor_type = "~", skips = [],
        err_rate = 0.03, fq_upper_bound = 0.75, comprehensive = False):
    """ Define SNVs over `units` in `reads`. Currently, SNV is dict with following keys: k, p, b, f, c, binom_p
        If `units` is not specified, it look for all unit of `hor_type`.
        If `units` nor `hor_type` is specified, it defaults to use "~".
        It assumes error rate of `err_rate` for calculation of FDR.
        Too frequent SNVs (>`fq_upper_bound`) are ignored as noninformative. """

    counter = Counter()

    # default to ~
    if not units:
        units = [ (ri, h) for ri, er in enumerate(reads) for h, _, t in er.hors if t == hor_type ]

    if not units:
        return []

    size = [ s for h, s, _ in reads[units[0][0]].hors if h == units[0][1] ][0]
    n_tests = 171 * 3 * (size - len(skips)) # approximate

    for ri, i in units:
        for j in [ j for j in range(size) if j not in skips ]:
            counter.update([ (j, s.pos, s.base) for s in reads[ri].mons[i+j].monomer.snvs ])

    # Here I'll return every detected variants positions!
    if comprehensive:
        return [ dict(k = k, p = p, b = b, f = c/len(units), c = c, binom_p = 1 - binom.cdf(c, len(units), err_rate)) for (k, p, b), c in counter.most_common() ]

    # remove too frequent ones, and allow <1 false positives
    if not comprehensive:
        _tmp = [ (k, p, b, c / len(units), c) for (k, p, b), c in counter.most_common() if c / len(units) < fq_upper_bound ]
        _nt_vars = [ dict(k = k, p = p, b = b, f = c/len(units), c = c, binom_p = 1 - binom.cdf(c, len(units), err_rate)) for k, p, b, f, c in _tmp ]
        return [ s for s in _nt_vars if s["binom_p"] * n_tests < 1.0 ]

def vec(read, h, snvs):
    """ bit array for a single HOR unit, starting at h-th monomers of a read """
    def has_s(s):
        # i.e., if s in t
        return any([ (t.pos, t.base) == (int(s["p"]), s["b"]) for t in read.mons[h+int(s["k"])].monomer.snvs ])
    return [ 1 if has_s(s) else -1 for s in snvs ]

def ba(read, arr, snvs):
    """ bit array for a single array """
    zeros, v = [0] * len(snvs), []

    #print(f"arr = {arr}\nlen(mons) = {len(read.mons)}")
    #print([ (i, m.monomer.name) for i, m in enumerate(read.mons) ])

    for h, s, t in arr:
        #print(h)
        if t == "~":
            v += vec(read, h, snvs)
        else:
            v += zeros
    return np.array(v).reshape(len(arr), len(snvs))

def ucomp(rd1, h1, rd2, h2, snvs = None):
    """ compare two units on `snvs`. use all columns if `snvs = None`.
        `rd1`, `rd2`: HOR encoded read
        `h1`, `h2`: start midx of HOR unit to be compared.
        NOTE: specialized for chrX
        """

    v1 = { (k, t.pos, t.base) for k in range(12) for t in rd1.mons[h1+k].monomer.snvs }
    v2 = { (k, t.pos, t.base) for k in range(12) for t in rd2.mons[h2+k].monomer.snvs }
    d = v1 ^ v2

    if snvs:
        _snv = { (int(s["k"]), int(s["p"]), s["b"]) for s in snvs }
        return len(d & _snv) / len(snvs)
    else:
        return len(d) / 2057.0

def acomp(rd1, a1, rd2, a2, snvs = None):
    """ compare two (arrays in) reads to generate matrix ready for dotplot.
        `rd1`, `rd2`: HOR encoded read
        `a1`, `a2`: arrays to be compared. returns from `fillx()`.
        NOTE: specialized for chrX
    """

    m = []
    for i, (h1, _, t1) in enumerate(a1):
        for j, (h2, _, t2) in enumerate(a2):
            if t1 == "~" and t2 == "~":
                m += [1.0 - ucomp(rd1, h1, rd2, h2, snvs)]
            elif t1 == t2: # variant unit match
                m += [-1]
            else: # variant unit inconsistency
                m += [-2]

    return np.array(m).reshape(len(a1), len(a2))

# TODO: use a generalized metric g # NOTE: reworking from impl in Alignment
# def mismatX(self, i, ai, j, aj, k, bits_dict):

def mismatX(i, ai, j, aj, b1, b2, k, g = None):
    """ aux func returns with some metadata a (3, len(snv)) matrix, which summarizes #11, #01/10, #00 of each column,
        in comparison of bit-represented read `b1` and `b2` with displacement of k units.
        g is the metric.
    """

    # li, lj = bits_dict[(i, ai)].shape[0], bits_dict[(j, aj)].shape[0]

    assert b1.shape[1] == b2.shape[1], "wrong shape"
    """
    if g == None:
        c11 = 1.0
        c00 = 0.1
        c01 = -2.0
        g = np.array(([c11] * b1.shape[1]) + ([c00] * b1.shape[1]) + ([c01] * b1.shape[1])).reshape(3, b1.shape[1])
    else:
        assert g.shape[0] == 3 and b1.shape[1] == g.shape[1], "wrong shape"
    """

    li, lj = b1.shape[0], b2.shape[0]
    lov = min(lj, -k+li) if k > 0 else min(li, k+lj) # length of overlap
    xi = b1[k:(k+lov),] if k > 0 else b1[0:lov,]
    xj = b2[0:lov,] if k > 0 else b2[-k:(-k+lov),]

    m = np.multiply(xi, xj) # 1: match, 0: masked, -1: mismatch
    match = np.multiply((m + 1), m) / 2 # 1: match, 0: masked/mismatch
    match_11 = np.multiply((xi + 1), match) / 2
    match_00 = np.multiply((1 - xi), match) / 2
    mismatch = np.multiply(m, m) - match # 1: mismatch (01/10), 0: masked/match

    # calculate score
    score = np.sum(np.multiply(g[0,], np.sum(match_11, axis=0)))
    score += np.sum(np.multiply(g[1,], np.sum(match_00, axis=0)))
    score += np.sum(np.multiply(g[2,], np.sum(mismatch, axis=0)))

    # effective length of overlap
    eov = lov - sum([ 1 if x == 0 else 0 for x in m[:,0] ])

    #return dict(li = li, lj = lj, lov = lov, eov = eov,
    #        n11 = n11, n00 = n00, match = np.sum(match), mismatch = np.sum(mismatch))
    #X = self.mismatX(i, ai, j, aj, k, bits_dict = bits_dict)
    #score = 1000 * (X["n11"] + self.c_00 * X["n00"]) / (0.01 + X["n11"] + self.c_00 * X["n00"] + self.c_ms * X["mismatch"]) # per-mille !

    return Aln(
        i = i, ai = ai, li = li, j = j, aj = aj, lj = lj, k = k,
        score = score, lov = lov, eov = eov, f_ext = max(0, k+lj-li), r_ext = max(0, -k),
        n11 = np.sum(match_11), n00 = np.sum(match_00), nmis = np.sum(mismatch))


def print_snvs(snvs, sort = "freq", n_tests = 2057, alt_snvs = None, innum = False):

    """ Show SNVs """

    lines = "midx\tpos\tbas\tcount\tfreq\t   p-val\t cfd."
    lines += "\t f.gl\n" if alt_snvs else "\n"

    if sort == "freq":
        key = lambda x: -x["c"]
    else:
        key = lambda x: (x["k"], x["p"])

    for s in sorted(snvs, key = key):
        k, p, b, c, f, pv = s["k"], s["p"], s["b"], s["c"], s["f"], s["binom_p"]

        if innum:
            b = {'A':0, 'C':1, 'G':2, 'T':3}[b]

        if alt_snvs:
            alt_f = [ a["f"] for a in alt_snvs if (k, p, b) == (a["k"], a["p"], a["b"]) ]
            alt_info = f"{100*alt_f[0]:>5.2f}" if alt_f else "    *"
            lines += f"{k:>3}\t{p:>3}\t{b}\t{c:>5}\t{100*f:.2f}\t{pv:.2e}\t{pv*n_tests:>7.2f}\t{alt_info}\n"
        else:
            lines += f"{k:>3}\t{p:>3}\t{b}\t{c:>5}\t{100*f:.2f}\t{pv:.2e}\t{pv*n_tests:>7.2f}\n"

    print(lines)

def print_align(aln, bits, arrs):
    """ Visualize an alignment in text to clarify the overlap. """
    pair_to_char = {
            ( 1, 1): "#", ( 1, 0): "0", ( 1, -1): ".",
            ( 0, 1): "0", ( 0, 0): "@", ( 0, -1): "0",
            (-1, 1): ",", (-1, 0): "0", (-1, -1): "_"}

    bit_to_char = { 1: "x", 0: "o", -1: " "}
    r, q = bits[aln.i], bits[aln.j]
    rs, qs = max(0, aln.k), max(-aln.k, 0)
    re, qe = rs + aln.lov, qs + aln.lov

    ra, qa = arrs[aln.i], arrs[aln.j]

    #if is_top:
    print(f"\n*\tri\tqi\trs\tre\trl\tqs\tqe\tql\tscr.\tovlp\teff.ovlp")
    print(f"ALIGN\t{aln.i}\t{aln.j}\t{rs}\t{re}\t{aln.li}\t{qs}\t{qe}\t{aln.lj}\t{aln.score:.3f}\t{aln.lov}\t{aln.eov}\n")

    n_snvs = r.shape[1]
    _l = f"== {n_snvs} SNVs =="
    lines = ("   \t  i\t{0:^" + f"{n_snvs}" + "}\tj  ").format(_l)

    # dangling part if any
    for i in range(rs):
        lines += f"\n{ra[i][2]:>3}\t{i:>3}\t" + "".join([ bit_to_char[e] for e in r[i,:] ]) + f"\t***\t***"
    for i in range(qs):
        lines += f"\n***\t***\t" + "".join([ bit_to_char[e] for e in q[i,:] ]) + f"\t{i:<3}\t{qa[i][2]:<3}"

    # the aligned part
    for i in range(re-rs):
        lines += f"\n{ra[rs+i][2]:>3}\t{rs+i:>3}\t" + "".join([ pair_to_char[(e, d)] for e, d in zip(r[rs+i,:], q[qs+i,:]) ]) + f"\t{qs+i:<3}\t{qa[qs+i][2]:<3}"

    # the remaining dangling part if any
    for i in range(re, aln.li):
        lines += f"\n{ra[i][2]:>3}\t{i:>3}\t" + "".join([ bit_to_char[e] for e in r[i,:] ]) + f"\t***\t***"
    for i in range(qe, aln.lj):
        lines += f"\n***\t***\t" + "".join([ bit_to_char[e] for e in q[i,:] ]) + f"\t{i:<3}\t{qa[i][2]:<3}"

    print(lines)

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='Analyze HOR encoded read.')
    parser.add_argument('action', metavar='action', type=str,
            help='action to perform: print-var, print-gap, print-snv-evolution...')
    parser.add_argument('--hor-reads', dest='hors', help='pickled hor-encoded long reads')

    # For print-var
    parser.add_argument('--hor-type', dest='hor_type', help='HOR unit on which variants will be reported.')
    parser.add_argument('--skips', dest='skips', help='idx of monomers to be ignored.')
    parser.add_argument('--all', dest='allvars', action='store_const', const=True, help='report all mismatches')
    parser.add_argument('--innum', dest='innum', action='store_const', const=True, help='bases are encoded into int (for plotting)')

    #parser.add_argument('-o', dest='outfile', help='the file to be output (required for align)')
    args = parser.parse_args()
    assert args.hors, "specify HOR-encoded reads"
    hers = pickle.load(open(args.hors, "rb"))
    hor_type = args.hor_type if args.hor_type else "~"
    skips = [ int(i) for i in args.skips.split(",") ] if args.skips else []

    # Print SNVs detected
    if args.action == "print-var":

        v = var(hers, hor_type = hor_type,
            fq_upper_bound = 1.1, skips = skips,
            comprehensive = True if args.allvars else False)

        print(f"# Variants on HOR units of type: {hor_type}")
        print_snvs(v, sort = "freq", innum = True if args.innum else False)

    elif args.action == "print-gap":

        for her in hers:
            print(her.name)
            fillx(her, verbose = True)

    elif args.action == "print-snv-evolution":

        v_major = var(hers, hor_type = hor_type,
            fq_upper_bound = 1.1, skips = skips,
            comprehensive = False)

        print(f"print-snv-evolution : {len(v_major)} SNVs\nd\t%MM-All\t%MM-SNV")

        for her in hers:

            arr = fillx(her)
            if not arr:
                continue

            # print(her.name)
            for d in range(1, 21):
                for i in range(len(arr) - d):
                    h, s, t = arr[i]
                    _h, _s, _t = arr[i+d]
                    if (t == "~") & (_t == "~"):
                        div = ucomp(her, h, her, _h)
                        divm = ucomp(her, h, her, _h, snvs = v_major)
                        #print(f"{d}\t{h}\t{_h}\t{100*div:.2f}\t{100*divm:.2f}")
                        print(f"{d}\t{100*div:.3f}\t{100*divm:.3f}")

    elif args.action == "test-align":

        v_major = var(hers, hor_type = "~", err_rate = 0.05,
            fq_upper_bound = 1.1, skips = skips,
            comprehensive = False)

        g_row_0 = [ s["f"] * (1.0 - s["f"]) for s in v_major ]
        g_row_1 = [ 0.1 * s["f"] * (1.0 - s["f"]) for s in v_major ]
        g_row_2 = [ -1.5 * s["f"] * (1.0 - s["f"]) for s in v_major ]
        g = np.array(g_row_0 + g_row_1 + g_row_2).reshape(3, len(v_major))

        # c11, c00, c01 = 1.0, 0.1, -2.0
        # g = np.array(([c11] * b1.shape[1]) + ([c00] * b1.shape[1]) + ([c01] * b1.shape[1])).reshape(3, b1.shape[1])

        arrs = [ fillx(her) for her in hers ]
        bits = { i: ba(her, arrs[i], v_major) for i, her in enumerate(hers) if arrs[i] and len(arrs[i]) > 9 }

        bits_keys_longer = sorted(bits.keys(), key = lambda x: -bits[x].shape[0])

        print(f"{len(bits)} bb long arrs in {len(hers)} reads", flush=True)

            # alns_dict[(i, ai)] = dict()

        def look_for_ext(i):
            # TODO: add description!!
            """ returns a dict of possible extension from read i.
                an entry (for a single read) is a list of `aln` sorted by the score. """

            ad0 = {}
            for j in [ j for j in bits.keys() if j != i ]:
                li, lj = bits[i].shape[0], bits[j].shape[0]
                # with at least 4 units in overlap
                res = [ mismatX(i, 0, j, 0, bits[i], bits[j], k, g) for k in range(-lj+4, li-4) ]
                res = [ r for r in res if r.score > -100 and r.eov > 4 ]
                ad0[j] = sorted(res, key=lambda x: -x.score)[:2]

            print(i, flush=True)

            return ad0

        ad = [ (i, look_for_ext(i)) for i in bits_keys_longer[:10] ]

        for i, ad0 in ad:
            print(f"# alignment for  {i}")
            
            for j in sorted([ j for j in ad0.keys() if ad0[j] ], key=lambda x: -ad0[x][0].score):
                print(f"- alignment {i} with {j}")
                for aln in ad0[j][:2]:
                    print_align(aln, bits, arrs)

            # sorted(res, key = lambda x: x.score
            # ad[0] |= { j : [ mismatX(0, 0, j, 0, bits[0], bits[j], k) for k in range(-lj+2, li-2) ] }

                #alns = [ self.calc_align(i, ai, j, aj, k, bits_dict = bits_dict) for k in range(-lj+2, li-2) ]
                #alns = [ aln for aln in alns if aln.score > T_dag - 100 and aln.eff_ovlp > 4 ] # NOTE: threshold depends on scoring scheme.

        """
        for i, ai in targets:
            if not quiet:
                print(".", end = "", flush = True)

            alns_dict[(i, ai)] = dict()
            for j, aj in [ (j, aj) for j, aj in queries if not (j, aj) == (i, ai)]:
                li, lj = bits_dict[(i, ai)].shape[0], bits_dict[(j, aj)].shape[0]
                # with at least 2 units in overlap
                alns = [ self.calc_align(i, ai, j, aj, k, bits_dict = bits_dict) for k in range(-lj+2, li-2) ]
                alns = [ aln for aln in alns if aln.score > T_dag - 100 and aln.eff_ovlp > 4 ] # NOTE: threshold depends on scoring scheme.
        """


        """
    if args.action == "align":
        assert args.hors, "specify HOR-encoded reads"
        assert args.outfile, "specify output file"

        hers = pickle.load(open(args.hors, "rb"))
        alignment = Alignment(hers)
        store = alignment.get_all_vs_all_aln()
        with open(args.outfile, "wb") as f:
            pickle.dump(store, f)
        print("done.")

    if args.action == "print":
        assert args.alns, "specify pickled alignment file"

        with open(args.alns, "rb") as f:
            store = pickle.load(f)
            alignment = Alignment(store.reads, arrs = store.arrays, variants = store.variants)
            store.describe(alignment.bits)
        """

    else:
        assert False, "invalid action."
