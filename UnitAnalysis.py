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

# TODO: I'll work on this.
# TODO: can this be abstracted to work with other chromosomes?

def fillx(read, verbose = False):

    """
    this is more like `fill`ing appropriate number of `*` HOR units.
    get an array from a read.
    """

    # TODO: check this impl would work
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
        print( hors_found[len(hors_found)-1] )
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
    v = [ [0] * len(snvs) if h < 0 else vec(read, h, snvs) for h in arr ]
    return np.array(v).reshape(len(arr), len(snvs))

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

    # Print SNVs detected
    if args.action == "print-var":
        hor_type = args.hor_type if args.hor_type else "~"
        skips = [ int(i) for i in args.skips.split(",") ] if args.skips else []

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

        print("print-snv-evolution")
        for her in hers:
            print(her.name)
            arr = fillx(her)

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
