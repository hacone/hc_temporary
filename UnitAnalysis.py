# Utilities to analyze HOR encoded reads before alignment.
# As a standalone script, this outputs relevant data and figures available at this stage.

from scipy.stats import binom
from collections import Counter
import hashlib
from collections import namedtuple
import numpy as np
import pandas as pd
import networkx as nx
import pickle

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns

import svgwrite

import random
random.seed(42)
import sys

# TODO: I just need datatype definition. That can be separated from other codes.
from EncodedRead import *
from HOR_segregation import *

t2c = {"*": "*", "~": "", "D1":"Y", "D12":"X", "22U":"U", "D39":"V", "D28":"W"}
t2col = {"*": "black", "~": "white",
         "D1":"yellow", "D12":"green",
         "22U":"brown", "D39":"orange", "D28":"red"}

def squarify(M,val):
    a , b = M.shape
    if a > b:
        padding = ((0,0),(0,a-b))
    else:
        padding = ((0,b-a),(0,0))
    return np.pad(M,padding,mode='constant',constant_values=val)

# this is an atomic alignment
Aln = namedtuple("Aln", ("koff", "score", "eov", "fext", "rext"))
# for each read pair i, j, calc one BestAln out of Aln's
BestAln = namedtuple("BestAln", ("i", "j", "aln", "gap", "nround"))
# plus, minus, embed are list of BestAln
Extension = namedtuple("Extension", ("i", "plus", "minus", "embed", "vf_history"))


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
        err_rate = 0.05, fq_upper_bound = 0.75, comprehensive = False):
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
        # no SNV filter
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
                #m += [-1]
                m += [np.nan]
            else: # variant unit inconsistency
                #m += [-2]
                m += [np.nan]

    return np.array(m).reshape(len(a1), len(a2))

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
    parser.add_argument('--vars', dest='vars', help='pickled variant sites (disable auto detection)')

    # For print-snv-evolution
    parser.add_argument('--rand', dest='random', action='store_const', const=True, help='compare 50k random pairs of units')
    parser.add_argument('--nvars', dest='nvars', help='number of variants to be used')

    # For print-var
    parser.add_argument('--hor-type', dest='hor_type', help='HOR unit on which variants will be reported.')
    parser.add_argument('--skips', dest='skips', help='idx of monomers to be ignored.')
    parser.add_argument('--all', dest='allvars', action='store_const', const=True, help='report all mismatches')
    parser.add_argument('--err-rate', dest='err_rate', help='error rate assumed in variants detection')
    # and --nvars as well
    parser.add_argument('--innum', dest='innum', action='store_const', const=True, help='bases are encoded into int (for plotting)')
    parser.add_argument('--save-to', dest='save', help='pickle variants for later use')

    # For layout, layout-2
    parser.add_argument('--reverse', dest='backwards', action='store_const', const=True, help='backwards search')
    parser.add_argument('--range', dest='sizeranges', help='size ranges to be used')
    # parser.add_argument('--parallel', dest='n_parallel', help='# of cores to be used') # deprecated
    parser.add_argument('--park', dest='park', help='for parallelly processing read id = k mod 24')
    parser.add_argument('--edges', dest='edgefile', help='edge list file')
    parser.add_argument('--layouts', dest='layouts', help='precomputed layouts to be analysed')


    #parser.add_argument('-o', dest='outfile', help='the file to be output (required for align)')
    args = parser.parse_args()
    assert args.hors, "specify HOR-encoded reads"
    hers = pickle.load(open(args.hors, "rb"))
    arrs = [ fillx(her) for her in hers ]
    hor_type = args.hor_type if args.hor_type else "~"
    units = [ (her, h) for her in hers for h, s, t in fillx(her) if t == hor_type ]
    skips = [ int(i) for i in args.skips.split(",") ] if args.skips else []
    import datetime
    datetimestr = datetime.datetime.now().strftime("%Y.%m%d.%H%M")

    # Print SNVs detected
    if args.action == "print-var":

        v = var(hers, hor_type = hor_type,
            fq_upper_bound = 1.1, skips = skips,
            err_rate = float(args.err_rate) if args.err_rate else 0.03,
            comprehensive = True if args.allvars else False)

        print(f"# Variants on HOR units of type: {hor_type}")
        print_snvs(v, sort = "freq", innum = True if args.innum else False)

        if args.save:
            if args.nvars:
                v = v[:int(args.nvars)]

            print(f"\n\n{len(v)} of these variants are to be saved to {args.save}")
            with open(args.save, "wb") as f:
                pickle.dump(v, f)

    elif args.action == "print-gap":

        for her in hers:
            print(her.name)
            fillx(her, verbose = True)

    elif args.action == "draw-units-pca":
        """ generate a figure of all units' PCA """

        if args.vars:
            v_major = pickle.load(open(args.vars, "rb"))
        else:
            v_major = var(hers, hor_type = hor_type,
                fq_upper_bound = 1.1, skips = skips,
                comprehensive = True)[:100]
                #comprehensive = False)

        bits = { i: ba(her, arrs[i], v_major) for i, her in enumerate(hers) if arrs[i] }

        a, n, ridx = [], 0, []
        for i in range(len(hers)):
            if not arrs[i] or (len(arrs[i]) < 10):
                ridx.append([n, n])
                continue

            rs = n
            for j, (h, s, t) in enumerate(arrs[i]):
                if t != hor_type:
                    continue
                n += 1
                a.extend([i, j])
                a.extend(list(bits[i][j,]))
            ridx.append([rs, n])

        X = np.array(a).reshape(n, 2 + len(v_major))
        X = X[0:10000,:]
        print("load data")

        from sklearn.cluster import KMeans
        from sklearn.decomposition import PCA

        for nvars in [25, 50, 75, 100]:
            Xs = X[:,2:2+nvars]
            km_model = KMeans(n_clusters=3, random_state=0, n_jobs=-3).fit(Xs)
            cls = km_model.predict(Xs)
            print("kmeans done")

            pca = PCA(n_components=2).fit(Xs)
            Xrd = pca.fit_transform(Xs)
            print("pca done")

            xlim = min(Xrd[:,0]), max(Xrd[:,0])
            ylim = min(Xrd[:,1]), max(Xrd[:,1])

            g = sns.scatterplot(x = Xrd[:,0], y = Xrd[:,1], hue=cls, s = 5, alpha = 0.8, palette = "Accent")
            plt.xlim(xlim)
            plt.ylim(ylim)
            plt.savefig(f"PCA-{X.shape[0]}-units-{nvars}-vars.png")
            plt.close()

        # for each read !
        nvars = min([int(args.nvars), 100]) if int(args.nvars) else 50
        Xs = X[:,2:2+nvars]
        km_model = KMeans(n_clusters=3, random_state=0, n_jobs=-3).fit(Xs)
        print("k-means done")
        cls = km_model.predict(Xs)
        pca = PCA(n_components=2).fit(Xs)
        print("pca done")
        Xrd = pca.fit_transform(Xs)

        xlim = min(Xrd[:,0]), max(Xrd[:,0])
        ylim = min(Xrd[:,1]), max(Xrd[:,1])

        # take reads upto 5000 units, sorting by x-coord in PCA.
        gi = [ i for i in range(len(hers)) if ridx[i][1] -ridx[i][0] > 9 and ridx[i][0] < Xrd.shape[0] ][:200]
        gi = sorted(gi, key = lambda i: Xrd[ridx[i][0],0])

        for j, i in enumerate(gi):
            g = sns.lineplot(
                    x = Xrd[ridx[i][0]:ridx[i][1],0],
                    y = Xrd[ridx[i][0]:ridx[i][1],1],
                    lw = 1, sort = False)
            plt.xlim(xlim)
            plt.ylim(ylim)
            plt.text(xlim[0], ylim[0], f"Read-{i:04}", horizontalalignment='left', size='medium', color='black', weight='semibold')
            plt.savefig(f"PCA-read-{nvars}-vars-{j:04}.png")
            plt.close()
            print("|" if j%10 else ".", flush = True)


    elif args.action == "print-snv-evolution":

        if args.vars:
            # if vars are specified.
            v_major = pickle.load(open(args.vars, "rb"))
        else:
            v_major = var(hers, hor_type = hor_type,
                fq_upper_bound = 1.1, skips = skips,
                comprehensive = False)
        if args.nvars:
            # if the number of SNVs is specified.
            v_major = var(hers, hor_type = hor_type,
                fq_upper_bound = 1.1, skips = skips,
                comprehensive = True)
            v_major = v_major[0:int(args.nvars)]

        n = 10000
        print(f"print-snv-evolution ({n} of {len(units)} units): {len(v_major)} SNVs")
        a = [] # accumulator for plot
        for i in range(n):
            j, k = random.randrange(0, len(units)), random.randrange(0, len(units))
            r, x = units[j]
            s, y = units[k]
            div = ucomp(r, x, s, y)
            divm = ucomp(r, x, s, y, snvs = v_major)
            a += [-1, 100*div, 100*divm]

        df = pd.DataFrame(
            data = np.array(a).reshape(n, 3),
            columns =["d", "All", "Snv"])

        a = [] # accumulator for plot
        for her in hers:
            arr = fillx(her)
            if not arr:
                continue
            for d in range(1, 11):
                for i in range(len(arr) - d):
                    h, s, t = arr[i]
                    _h, _s, _t = arr[i+d]
                    if (t == "~") & (_t == "~"):
                        div = ucomp(her, h, her, _h)
                        divm = ucomp(her, h, her, _h, snvs = v_major)
                        a += [d, 100*div, 100*divm]

        df = df.append(pd.DataFrame(
            data = np.array(a).reshape(int(len(a)/3), 3),
            columns = ["d", "All", "Snv"]))


        df = pd.melt(df, id_vars = ["d"], var_name = "class", value_name = "pDiv") #, value_vars = data.columns[:-2].tolist())

        g = sns.boxplot(x='d', y="pDiv", hue='class', data=df, palette="PRGn")
        g.set_xticklabels(["Random"] + [ f"{i+1}" for i in range(10) ], rotation = -45)
        plt.savefig(f"units-divergence-{len(v_major)}vars.svg")
        plt.close()

    elif args.action == "layout":
        """ This action calculates list of viable edges using proper measure
            along with some nice figures. """

        ## common data structures and helper function to calculate proper similarity scores
        class error_tensor:

            def __init__(self, e0 = 0.01, e1 = 0.05):
                self.e0 = e0
                self.e1 = e1
                self.td = np.array([1-e0, e0, e1, 1-e1]).reshape(2, 2)
                self.dsq = np.array(
                        [ e0*e0 + (1-e0)*(1-e0), e0*(1-e1) + e1*(1-e0), 
                          e0*(1-e1) + e1*(1-e0), e1*e1 + (1-e1)*(1-e1) ]).reshape(2, 2)

        def b_to_plup(b, f):
            # TODO: this will do?
            #return np.stack([0.5*(1-b) - (1-f), 0.5*(1+b) - f])
            return np.stack([0.5*(1-b) - (1-f), 0.5*(1+b) - f]).reshape(2, b.shape[0])

        def denomi(plup, errt):
            # first; pluppp = plup[n,k] * dsq[l,n]
            # then; denomi = plup[l,k] * pluppp[l,k]
            pluppp = np.tensordot(errt.dsq, plup, [1, 0])
            return np.tensordot(plup, pluppp, [[0,1],[0,1]])

        def numeri(plup, bx, errt, t2_ve):
            x = (1 + np.outer(np.array([-1, 1]), bx)) / 2
            t1 = np.tensordot(x, errt.td, [0, 1])
            return np.tensordot(plup, t1 - t2_ve, [[0,1],[1,0]])

        ## entry point of this action.

        if args.vars:
            v_major = pickle.load(open(args.vars, "rb"))
        else:
            v_major = var(hers, hor_type = "~", err_rate = 0.05,
                fq_upper_bound = 1.1, skips = skips,
                comprehensive = False)
        vf = np.array([ v["f"] for v in v_major ]).reshape(len(v_major))

        v_all = var(hers, hor_type = "~", err_rate = 0.05,
            fq_upper_bound = 1.1, skips = skips,
            comprehensive = True)
        vf_all = np.array([ v["f"] for v in v_all ]).reshape(len(v_all))

        errt = error_tensor(e0 = 0.03, e1 = 0.10) # for PacBio
        # errt = error_tensor(e0 = 0.06, e1 = 0.20) # for ONT?
        # errt_lev = error_tensor(e0 = 0.10, e1 = 0.10)
        # vf_cst = np.array([ 0.1 for v in v_major ]).reshape(len(v_major))

        print("\ntd =")
        print(errt.td)

        print("\ndsq =")
        print(errt.dsq)

        print(f"\nvf = {len(v_major)} global vars:")
        print(vf)
        print(f"\nvf_all = {len(v_all)} comprehensive vars: suppressed...")

        koff = 0
        origin = 0

        bits = { i: ba(her, arrs[i], v_major) for i, her in enumerate(hers) if arrs[i] and len(arrs[i]) > 1 }

        keys_sorted = sorted(
                [ i for i in range(len(hers)) if arrs[i] and len(arrs[i]) > 1 ],
                key = lambda x: -1 * len(arrs[x]))
        sizeranges = [ int(n) for n in args.sizeranges.split(",") ] if args.sizeranges else [8, 0]
        keys_long = [ i for i in keys_sorted if len(arrs[i]) >= sizeranges[0] ]
        keys_short = [ i for i in keys_sorted
                if len(arrs[i]) < sizeranges[0] and len(arrs[i]) >= sizeranges[1] ]

        print(f"{len(keys_long)} long (>={sizeranges[0]}) arrs, {len(keys_short)} short (>={sizeranges[1]}) arrs, of total {len(hers)} reads")

        print("\n#\tid\tmons\tunits\treadname")
        print("\n".join([ f"{n+1}\t{i}\t{len(hers[i].mons)}\t{bits[i].shape[0]}\t{hers[i].name}" for n, i in enumerate(keys_sorted) if len(arrs[i]) > 7 ]))

        def dotplot(i, j, errt = errt, vf = vf, bits = bits,
                naive = False, vs = []):
            """ calculate each cell of dotplot of prop. sim.
                required context: arrs """

            li = bits[i].shape[0]
            lj = bits[j].shape[0]
            result = np.zeros(li*lj).reshape(li, lj)
            result[:] = np.nan

            if naive:
                ## Naive comparison; just counting the mismatched SNVs
                si = { n: { (k, snv.pos, snv.base)
                            for k in range(s) for snv in hers[i].mons[h+k].monomer.snvs }
                       for n in range(li) for h, s, t in [arrs[i][n]] if t == "~" }

                sj = { n: { (k, snv.pos, snv.base)
                            for k in range(s) for snv in hers[j].mons[h+k].monomer.snvs }
                       for n in range(lj) for h, s, t in [arrs[j][n]] if t == "~" }

                vs_set = { (v["k"], v["p"], v["b"]) for v in vs }

                for n in range(li):
                    if arrs[i][n][2] != "~":
                        continue
                    for m in range(lj):
                        if arrs[j][m][2] != "~":
                            continue
                        if vs:
                            result[n,m] = 1 - (len((si[n] ^ sj[m]) & vs_set) / len(vs_set))
                        else:
                            result[n,m] = 1 - (len(si[n] ^ sj[m]) / 2057.0)
                return result

            t2_ve = np.tensordot(np.stack([1-vf, vf]), errt.dsq, [0, 1])
            plup = [ b_to_plup(bits[i][n,], vf) for n in range(li) ]
            denomis = [ denomi(plup[n], errt) for n in range(li) ]

            for n in range(li):
                if arrs[i][n][2] != "~":
                    continue
                for m in range(lj):
                    if arrs[j][m][2] != "~":
                        continue
                    result[n,m] = numeri(plup[n], bits[j][m,], errt, t2_ve) / (denomis[n] + 0.0001)

            return result

        # TODO: temporary naming, obviously
        def get_result_i(i, vf = vf, bits = bits, targets = None, verbose = False):
            """ assumes errt in context """

            # aln for i
            li = bits[i].shape[0]
            plup = [ b_to_plup(bits[i][n,], vf) for n in range(li) ]
            denomis = [ denomi(plup[n], errt) for n in range(li) ]
            t2_ve = np.tensordot(np.stack([1-vf, vf]), errt.dsq, [0, 1])

            result_i = {}
            targets = targets if targets else keys_long
            print(f"aligning for {i} : {len(targets)} targets")

            for j in [ j for j in targets if i != j ]:
                # aln for i, j
                lj = bits[j].shape[0]

                result_j = []
                maxs_j = -10000

                for k in range(-lj+2, li-1): ## minimum overlap is set to be 2
                    # aln for i, j, k
                    s = 0 # score
                    e = 0 # units effectively overlapped

                    if k > 0:
                        lov = min(lj, -k+li) # length of overlap
                        # NOTE: xi = bits[i][k:(k+lov),] is aligned to xj = bits[j][0:lov,]
                        for n in [ n for n in range(lov) if arrs[i][k+n][2] == "~" and arrs[j][n][2] == "~"]:
                            e += 1
                            s += numeri(plup[k+n], bits[j][n,], errt, t2_ve) / (denomis[k+n] + 0.0001)

                    else:
                        lov =  min(li, k+lj) # length of overlap
                        # NOTE: xi = bits[i][0:lov,] is aligned to xj = bits[j][-k:(-k+lov),]
                        for n in [ n for n in range(lov) if arrs[i][n][2] == "~" and arrs[j][n-k][2] == "~"]:
                            e += 1
                            s += numeri(plup[n], bits[j][n-k,], errt, t2_ve) / (denomis[n] + 0.0001)

                    if e >= 2:
                        s = s / e
                        maxs_j = max(maxs_j, s)
                        f_ext, r_ext = max(0, k+lj-li), max(0, -k)
                        result_j += [ Aln(koff = k, score = s, eov = e, fext = f_ext, rext = r_ext) ] # offset, score, eov, fext, rext

                if result_j:
                    result_i[j] = sorted(result_j, key = lambda x: -x.score)

                if verbose:
                    print(".", end = "", flush = True)

            return result_i

        # layout is like [(j0, k0), (j1, k1), (j2, k2)...]
        # extension into the right
        def extension(i, backwards = False):

            vf_history = []
            bits_history = []

            # initialize for round 0
            nround = 0
            best_long_js = keys_long
            best_short_js = keys_short

            # TODO: make this deprecated
            n_long_target = len(keys_long)
            n_short_target = len(keys_short)

            vs_local = v_major
            vf_local = vf
            bits_local = bits
            vf_history += [vf_local]
            bits_history += [bits_local]

            plus = []
            minus = []
            embed = []

            #global fig_idx
            #global ctg_idx

            # NOTE: n_units must be >100
            while nround < 15 and \
                  best_long_js and best_short_js and \
                  vs_local and (vs_local[0]["c"] / vs_local[0]["f"]) > 100:

                result_i_long = get_result_i(i, vf_local, bits_local, targets = best_long_js) # ~7 us
                best_long_js = sorted(result_i_long.keys(), key = lambda j: (-1) * result_i_long[j][0].score)

                result_i_short = get_result_i(i, vf_local, bits_local, targets = best_short_js) # 2~7 us
                best_short_js = sorted(result_i_short.keys(), key = lambda j: (-1) * result_i_short[j][0].score)

                def gap(alns):
                    if not [ a for a in alns if a.score > 0.6 ]:
                        return 0.0
                    elif len(alns) < 2:
                        return 1.0
                    else: 
                        return (alns[0].score - alns[1].score) / alns[0].score

                uniques = sorted([ j for j in best_long_js if gap(result_i_long[j]) > 0.30 ], key = lambda x: -result_i_long[x][0].eov)

                def pickbest(l):
                    # l is aln
                    _l, li = [], list(range(len(l)))
                    while li:
                        _l += [ l[li[0]] ]
                        li = [ i for i in li if l[i].j not in [ e.j for e in _l ]]
                    return _l

                ## j, aln, gap, round
                plus += [ BestAln(i = i, j = j, aln = result_i_long[j][0], gap = gap(result_i_long[j]), nround = nround)
                          for j in uniques if result_i_long[j][0].fext > 0 and result_i_long[j][0].eov > 4 ]
                plus = pickbest(sorted(plus, key = lambda x: (-x.gap, -x.aln.eov)))

                minus += [ BestAln(i = i, j = j, aln = result_i_long[j][0], gap = gap(result_i_long[j]), nround = nround)
                           for j in uniques if result_i_long[j][0].rext > 0 and result_i_long[j][0].eov > 4 ]
                minus = pickbest(sorted(minus, key = lambda x: (-x.gap, -x.aln.eov)))

                embed += [ BestAln(i = i, j = j, aln = result_i_long[j][0], gap = gap(result_i_long[j]), nround = nround)
                           for j in uniques if result_i_long[j][0].fext == 0 and result_i_long[j][0].rext == 0 and result_i_long[j][0].eov > 4 ]
                embed = pickbest(sorted(embed, key = lambda x: (-x.gap, -x.aln.eov)))

                line = f"\n(+{len(plus)}, -{len(minus)}, ^{len(embed)}) uniques for {i} upto round {nround}."
                line += f"\n   i\t   j\t  k\t  e\tfex\trex\tscr\t%gap\tround"
                for balns in plus[:10] + minus[:10] + embed[:10]:
                    line += f"\n{balns.i:4d}\t{balns.j:4d}\t{balns.aln.koff:3d}\t{balns.aln.eov:3d}\t" + \
                            f"{balns.aln.fext:3d}+\t{balns.aln.rext:3d}-\t{balns.aln.score:.2f}\t{100*balns.gap:.1f} %\t{balns.nround}"

                print(line, flush=True)

                # NOTE: update target, vf, bits
                n_long_target = int(n_long_target*0.6)
                n_short_target = int(n_short_target*0.6)

                # take ones with plus extension
                # TODO: also for shorter reads!
                if backwards:
                    best_long_js = [ j for j in best_long_js if result_i_long[j][0].rext > 0 ][:n_long_target]
                    best_short_js = [ j for j in best_short_js if result_i_short[j][0].fext == 0 ][:n_short_target]
                else:
                    best_long_js = [ j for j in best_long_js if result_i_long[j][0].fext > 0 ][:n_long_target]
                    best_short_js = [ j for j in best_short_js if result_i_short[j][0].koff >= 0 ][:n_short_target]

                nround += 1
                vs_local = var([ hers[j] for j in best_long_js + best_short_js ],
                    hor_type = "~", err_rate = 0.05,
                    fq_upper_bound = 1.1, skips = skips,
                    comprehensive = False)

                vf_local = np.array([ v["f"] for v in vs_local ]).reshape(len(vs_local))
                bits_local = { j: ba(hers[j], arrs[j], vs_local) for j in [i] + best_long_js + best_short_js }
                vf_history += [vf_local]
                bits_history += [bits_local]

                if vs_local:
                    est_units = int(vs_local[0]["c"] / vs_local[0]["f"])
                    print(f"Next for {i}: {len(vs_local)} local SNVs found for ~{est_units} units from {len(best_long_js)} long / {len(best_short_js)} short reads.")


            print("\n".join([ f"{i}\t{aln.j}\t{aln.aln.koff}\t{aln.aln.eov}\t{aln.aln.fext}\t{aln.aln.rext}\t" +\
                              f"{100*aln.aln.score:.2f}\t{100*aln.gap:.2f}\t{aln.nround}" for aln in plus + minus + embed ]),
                  file=open(f"read-{i}.{'rev' if backwards else 'fwd'}.ext", "w"), flush=True)

            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            sns.set(font_scale=1.5)

            if backwards:
                alns_to_plot = [(embed, "embed-rev")]
                alns_to_plot += [(minus, "minus")]
            else:
                alns_to_plot = [(embed, "embed")]
                alns_to_plot += [(plus, "plus")]

            ylabels = [ t2c[t] for h, s, t in arrs[i] ]
            for alns, name in alns_to_plot:
                for aln in alns[:3]:
                    # j, eov+fex, %gap
                    for nr in range(-1, aln.nround + 1):
                        fig = plt.figure(figsize=(15, 12))
                        # plt.axis("off")
                        ax1 = fig.add_subplot(1, 1, 1)
                        if nr == -1:
                            dots = squarify(dotplot(i, aln.j, vf = vf_history[aln.nround],
                                bits = bits_history[aln.nround], naive = True), np.nan)
                        elif nr == 0:
                            dots = squarify(dotplot(i, aln.j, vf = vf_history[aln.nround],
                                bits = bits_history[aln.nround], naive = True, vs = v_major), np.nan)
                        else:
                            dots = squarify(dotplot(i, aln.j, vf = vf_history[nr],
                                bits = bits_history[nr]), np.nan)
                        xlabels = [ t2c[t] for h, s, t in arrs[aln.j] ]
                        g1 = sns.heatmap(dots,
                                vmin=np.nanmin(dots), vmax=np.nanmax(dots),
                                cmap="YlGn" if nr > 0 else "Reds", ax = ax1,
                                xticklabels = xlabels, yticklabels = ylabels)
                                #cmap="coolwarm" if nr > 0 else , ax = ax1)

                        ax1.set_title(f"{i}-{aln.j}; @{aln.aln.koff}~{aln.aln.eov}+{aln.aln.fext}-{aln.aln.rext}; R{nr if nr > -1 else 'x'}/{aln.nround}; Score={100*aln.aln.score:.1f}% SG={100*aln.gap:.1f}%.")
                        plt.savefig(f"{i}-to-{aln.j}-{name}-R{nr if nr > -1 else 'x'}.png")
                        plt.close()

            exts = Extension(i = i, plus = plus, minus = minus, embed = embed, vf_history = vf_history)
            return exts

        ## NOTE: always parallelised by 24 
        exts = [ extension(i, backwards = args.backwards) for i in keys_long if i % 24 == int(args.park) ]
        outfile = f"exts-rev-mod{int(args.park)}.pickle" if args.backwards else f"exts-fwd-mod{int(args.park)}.pickle"
        pickle.dump(exts, open(outfile, "wb"))

    elif args.action == "layout-2":

        v_major = var(hers, hor_type = "~", err_rate = 0.05,
            fq_upper_bound = 1.1, skips = skips,
            comprehensive = False)
        v_all = var(hers, hor_type = "~", err_rate = 0.05,
            fq_upper_bound = 1.1, skips = skips,
            comprehensive = True)

        nonstd = { i:  [ (ii, t) for ii, (h, s, t) in enumerate(a) if t != "~" ]
                   for i, a in enumerate(arrs) }

        Edge = namedtuple("Edge", ("From", "To", "K", "Eov", "Fext", "Rext", "Score", "Gap", "Round"))
        df = pd.read_csv(args.edgefile, sep="\t", names = ["From", "To", "K", "Eov", "Fext", "Rext", "Score", "Gap", "Round"])
        G = nx.MultiDiGraph()
        #G = nx.MultiGraph()

        for _k, row in df.iterrows():
            if row.Gap > 30 or (row.Round > 4 and row.Gap > 25):
                G.add_edge(int(row.From), int(row.To), key = int(row.K))
                for k, v in row.items():
                        G.edges[int(row.From), int(row.To), int(row.K)][k] = v

        print(f"Loaded {len(G.nodes)} nodes.")

        def merge_k(l1, l2, k):
            """ merge 2 layouts by offset k """
            return sorted(l1 + [ (n, k+l) for n, l in l2 ], key = lambda x: x[1])

        def bridge(l1, l2, ins = None):
            """ estimate offset from l1 to l2 ...
                here, ins is the edges within/between l1/l2 
            """
            ls1 = { li: lk for li, lk in l1 }
            ls2 = { li: lk for li, lk in l2 }
            _edges = [ (gi, gj, k) for gi, gj, k in (ins if ins else G.edges(keys = True)) if gi in ls1 or gj in ls2 ]

            bridges = Counter()
            bridges.update([ ls1[gi] - ls2[gj] + k for gi, gj, k in _edges if gi in ls1 and gj in ls2 ])
            bridges.update([ ls1[gj] - ls2[gi] - k for gi, gj, k in _edges if gj in ls1 and gi in ls2 ])
            print(f"{bridges.most_common()[:3]} of {sum(bridges.values())}")

            return bridges

        def sanitize(l, ins = None, dryrun = False):
            """ sanity of the layout, returning list of layouts """

            if len(l) < 2:
                return [l]

            print(f"sanitizing a layout with {l[0]} of length {len(l)}")
            ret = []

            ins = ins if ins else G.edges(keys = True)
            # # global blacklist
            ls = { li: lk for li, lk in l } # set representation
            self_edges = [ (gi, gj, k) for gi, gj, k in ins if gi in ls and gj in ls ]
            # for each node in layout, list and number of (in)consystencies
            ics = { li[0]: [ (ls[gj] - ls[gi] - k, gj if gi == li[0] else gi)
                          for gi, gj, k in self_edges if li[0] in [gi, gj]]
                    for li in l }

            nics = { li[0]: len([ cp for cp in ics[li[0]] if cp[0] != 0 ]) for li in l }
            bads = sorted([ li[0] for li in l if ics[li[0]]], key = lambda x: -(nics[x] / len(ics[x])))

            # print(ics)
            print(nics)
            print(bads)

            while bads and nics[bads[0]]:
                l = [ li for li in l if li[0] != bads[0] ] # remove the worst
                ret += [[ (bads[0], 0) ]]

                ls = { li: lk for li, lk in l } # set representation
                self_edges = [ (gi, gj, k) for gi, gj, k in (ins if ins else G.edges(keys = True)) if gi in ls and gj in ls ]
                # for each node in layout, list and number of (in)consystencies
                ics = { li[0]: [ (ls[gj] - ls[gi] - k, gj if gi == li[0] else gi)
                              for gi, gj, k in self_edges if li[0] in [gi, gj]]
                        for li in l }
                nics = { li[0]: len([ cp for cp in ics[li[0]] if cp[0] != 0]) for li in l }
                bads = sorted([ li[0] for li in l if ics[li[0]]], key = lambda x: -(nics[x] / len(ics[x])))
                print(",".join([ f"{bi}:{nics[bi]}/{len(ics[bi])}" for bi in bads]))

            if ret:
                print("fluctured:" + f"{len(ret)} + {len(l)}")
                print(ret)
            return [l] + ret

        def consensus(layout):
            """ take a consensus of the layout """
            m, M = 10000, -10000
            types = { i : Counter() for i in range(-10000, 10000) }

            for li, lk in layout:
                #print(f"arr {li}@{lk} = {arrs[li]}")
                for ii, (h, s, t) in enumerate(arrs[li]):
                    M, m = max(M, lk + ii), min(m, lk + ii)
                    types[lk+ii].update([t])
                    #print(f"{ii}, {(h, s, t)} => {ii} + {lk} = {ii+lk}")

            units = { i: { k : Counter() for k in range(12) } for i in range(m, M+1) }
            depth = { i: 0 for i in range(m, M+1) }

            for li, lk in layout:
                if all([ True or (t == types[lk+ii].most_common()[0][0]) or t == "*"
                         for ii, (h, s, t) in enumerate(arrs[li]) ]):
                    for ii, (h, s, t) in enumerate(arrs[li]):
                        if t == "~":
                            depth[lk+ii] += 1
                            for mk in range(0, 12):
                                units[lk+ii][mk].update(hers[li].mons[h + mk].monomer.snvs)

            cons_hors = [ (i * 12, 12, types[i+m].most_common()[0][0]) for i in range(0, M-m+1) ] 
            cons_mons = [ AssignedMonomer(begin=lk*2057 + k*171, end=lk*2057 + (k+1)*171, ori="+",
                monomer = Monomer(
                    name = f"Consensus-{lk}-{k}",
                    snvs = [ s for s, sc in units[lk][k].most_common() if sc / depth[lk] > 0.4 ]))
                for lk in range(m, M+1) for k in range(12) ]

            print("recall, specificity?, fpr, fnr")
            print("\ni\ttps(%)\ttns(%)\tfps(%)\tfns(%)\tnvars\tdepth\ttypes observed.")

            ## total (default) units 
            coverage = sum([ depth[lk] for lk in range(m, M+1) ])

            vfs = { lk: 
                [ sc for k in range(12) for s, sc in units[lk][k].most_common() ]
                if depth[lk] > 0 and cons_hors[lk-m][2] == "~" else []
                for lk in range(m, M+1) }

            #print(vfs)
            name = layout[0][0]
            for lk in range(m, M+1):
                if depth[lk] > 0 and vfs[lk]:
                    print(f"VARS\t{name}\t{lk}\t{depth[lk]} => {Counter(vfs[lk])}")
                else:
                    print(f"VARS\t{name}\t{lk}\t* => *")

            print(f"name\tlk\tsc\td\tfq")
            for lk in range(m, M+1):
                if depth[lk] > 0 and vfs[lk]:
                    tps = [ sc for sc in vfs[lk] if sc / depth[lk] > 0.5  ]
                    fps = [ sc for sc in vfs[lk] if sc / depth[lk] <= 0.5 ]
                    sensitivity = sum(tps) / (depth[lk]*len(tps)) if tps else 0
                    specificity = 1 - (sum(fps) / (depth[lk]*(2057-len(tps)))) if fps else 0
                    #print("\n".join([
                    #    f"VARS\t{name}\t{lk}\t{sc}\t{depth[lk]}\t{100*sc/depth[lk]:.2f}"
                    #    for sc in sorted(vfs[lk], key = lambda x: -x/depth[lk])]))
                    print(f"{lk}" +\
                          f"\t{sum(tps)}/{len(tps)*depth[lk]} ({100*sensitivity:.2f})" +\
                          f"\t{sum(fps)}/{(2057-len(tps))*depth[lk]} ({100*specificity:.2f})" +\
                          f"\t{len(tps)}\t{depth[lk]}\t" +\
                          "\t".join([ f"{t}:{c}" for t, c in types[lk].most_common()]))
                else:
                    print(f"{lk}" +\
                          f"\t*\t*\t*\t*\t" +\
                          "\t".join([ f"{t}:{c}" for t, c in types[lk].most_common()]))


            hor_read = HOR_Read(name = "Consensus",
                            mons = cons_mons, hors = cons_hors,
                            length = (M-m)*2057, ori="+")

            print("readname\tbegin\tend\tidx\tsize\telem\tgap\tvars")
            print_HOR_read(hor_read)
            return hor_read

        def show_layout(layout):
            """ visualize consensus read """

            cons_read = consensus(layout)
            cons_arr = fillx(cons_read)
            snvs = var([ hers[li] for li, lk in layout ])
            print_snvs(snvs)

            snvs_read = var([cons_read], err_rate = 0.01, comprehensive = False)
            print(f"Total {len(snvs_read)} SNVs on this layout.")
            print_snvs(snvs_read)

            # t2c = {"*": "*", "~": "", "D1":"Y", "D12":"X", "22U":"U", "D39":"V", "D28":"W"}

            labels = [ t2c[t] for h, s, t in cons_read.hors ]

            # using reads SNVs
            for vs, name in [(v_all, "global"), (v_major, "gl-freq"), (snvs, "local"), (snvs_read, "consensus")]:
            # for vs, name in [(v_major, "gl-freq"), (snvs, "local"), (snvs_read, "consensus")]:
                dots = acomp(cons_read, cons_arr, cons_read, cons_arr, snvs = vs)
                fig = plt.figure(figsize=(20, 16))
                sns.set(font_scale=2)
                # plt.axis("off")
                ax1 = fig.add_subplot(1, 1, 1)
                g1 = sns.heatmap(dots,
                        vmin=np.nanmin(dots), vmax=np.nanmax(dots),
                        cmap="coolwarm", ax = ax1,
                        xticklabels = labels, yticklabels = labels)
                ax1.set_title(f"Layout-{layout[0][0]}; {len(layout)} rds; {len(cons_arr)} units; with {len(vs)} of {name} vars")
                plt.savefig(f"Layout-{layout[0][0]}-{name}.png")
                plt.close()

            dwg = svgwrite.Drawing(filename=f"Layout-{layout[0][0]}-str.svg")
            lkmin = min([ lk for li, lk in layout ])

            for i, (li, lk) in enumerate(layout):
                read = dwg.add(dwg.g(id=f"rd{li}", stroke='green'))
                #fig = plt.figure(figsize=(10, 8))
                #ax1 = fig.add_subplot(1, 1, 1)
                #ax1.set_xlim(-1, 1500)
                #ax1.set_ylim(-1, 300)
                for n, (h, s, t) in enumerate(arrs[li]):
                    read.add(
                        dwg.rect(insert=((n + lk - lkmin)*5,i*5), size=(4,4),
                        fill=t2col[t], stroke='black', stroke_width=0.5))
                    # ax1.text((lk + n) * 1, i * 1, t2c[t], fontsize=9) 
            dwg.save()

            #plt.savefig(f"Layout-{layout[0][0]}-str.svg")
            #plt.close()

        # if layouts is given
        if args.layouts:

            # say something about global mismatch distribution.
            layouts = pickle.load(open(args.layouts, "rb"))
            layouts = sorted(layouts, key = lambda x: -len(x))

            for i in range(len(layouts)): 
                if len(layouts[i]) > 4:
                    print(f"consensus for {i}")
                    show_layout(layouts[i])

            sys.exit()

        # Initialize
        layouts = [ [(n, 0)] for n in list(G.nodes) ]

        ins = [ (gi, gj, k)
            for l in layouts
            for gi, gj, k in G.edges(keys = True)
            if len(l) > 1 and {gi, gj} <= { li[0] for li in l } ]

        outs = [ (gi, gj, k) for gi, gj, k in G.edges(keys = True) if (gi, gj, k) not in set(ins) ]
        # eov sort, filter
        outs = sorted(outs, key = lambda x: -G.edges[x]["Gap"])

        for eov_t in [9, 8, 7, 6, 5, 4]:

            layouts = list(chain.from_iterable([ sanitize(l) for l in layouts ]))

            ins = [ (gi, gj, k)
                for l in layouts
                for gi, gj, k in G.edges(keys = True)
                if len(l) > 1 and {gi, gj} <= { li[0] for li in l } ]
            outs = [ (gi, gj, k) for gi, gj, k in G.edges(keys = True) if (gi, gj, k) not in set(ins) ]
            outs = sorted(outs, key = lambda x: -G.edges[x]["Gap"])
            # NOTE: for ONT
            outs = [ o for o in outs if G.edges[o]["Eov"] > eov_t * 5 and G.edges[o]["Gap"] > eov_t * 5 ]
            #outs = [ o for o in outs if G.edges[o]["Eov"] > eov_t and G.edges[o]["Gap"] > eov_t * 5 ]

            while outs:

                print(f"{len([l for l in layouts if len(l) > 1])} non-trivial layouts, {len(ins)} ins, {len(outs)} outs;")
                print([ (o, G.edges[o]["Gap"], G.edges[o]["Eov"]) for o in outs[:10] ])

                gi, gj, k = outs[0]
                gil = [ l for l in layouts if gi in {li[0] for li in l} ]
                gjl = [ l for l in layouts if gj in {li[0] for li in l} ]
                others = [ l for l in layouts if not len({gj, gi} & {li[0] for li in l}) ]

                assert len(gil) == 1, f"len(gil) = {len(gil)} is not 1"
                assert len(gjl) == 1, f"len(gjl) = {len(gjl)} is not 1"

                bs = bridge(gil[0], gjl[0])
                bsmc = bs.most_common()[0]

                print(f"Trying to merge L{gil[0][0][0]} = {len(gil[0])} and L{gjl[0][0][0]} = {len(gjl[0])}")
                print(f"{bsmc} / {sum(bs.values())} = {100*bsmc[1] / sum(bs.values()):.2f}% support.")
                print(f"{bsmc} / {sum([ x[1] for x in bs.most_common()[:2] ])}" +\
                      f"= {100*bsmc[1] / sum([ x[1] for x in bs.most_common()[:2] ]):.2f}% mergin.", flush = True)
                
                #if bsmc[1] / sum(bs.values()) < 0.8:
                if bsmc[1] / sum([ x[1] for x in bs.most_common()[:2] ]) < 0.9:
                    print("do not reject alternative; continue...", flush = True)
                    outs = outs[1:]
                    continue
                #if (len(gil[0]) + len(gjl[0])) > 2 and bsmc[1] < 2:
                #    print("not enough support continue...", flush = True)
                #    outs = outs[1:]
                #    continue

                merged = merge_k(gil[0], gjl[0], bsmc[0])

                # Check if merged is OK
                m, M = 10000, -10000
                for li, lk in merged:
                    for ii, (h, s, t) in enumerate(arrs[li]):
                        M, m = max(M, lk + ii), min(m, lk + ii)
                types = { i : Counter() for i in range(m, M+1) }

                for li, lk in merged:
                    for ii, (h, s, t) in enumerate(arrs[li]):
                        if t in ["D39", "D28"]:
                            types[lk+ii].update(["5mer"])
                        elif t in ["D1", "D12"]:
                            types[lk+ii].update(["10/11mer"])
                        elif t in ["22U"]:
                            types[lk+ii].update(["22mer"])
                        elif t == "~":
                            types[lk+ii].update(["def"])

                if any([ len(types[i].most_common(n=2)) > 1 for i in range(m, M+1)]):
                    print("Conflict! continue...", flush = True)
                    outs = outs[1:]
                    continue

                layouts = others + [merged]

                print(f"layout extended; {len(merged)} rds = {merged}", flush = True)
                # update ins, outs 
                ins = list(set(ins + [ (gi, gj, k)
                    for gi, gj, k in G.edges(keys = True)
                    if {gi, gj} <= { li[0] for li in merged } ]))

                outs = [ (gi, gj, k) for gi, gj, k in outs if (gi, gj, k) not in set(ins) ]
                outs = sorted(outs, key = lambda x: -G.edges[x]["Gap"])
                outs = [ o for o in outs if G.edges[o]["Eov"] > eov_t and G.edges[o]["Gap"] > eov_t * 5 ]

                print(f"{len([l for l in layouts if len(l) > 1])} non-trivial layouts, {len(ins)} ins, {len(outs)} outs.", flush = True)

            print(f"{len([l for l in layouts if len(l) > 1])} non-trivial layouts, {len(ins)} ins, {len(outs)} outs at the end of the round for {eov_t}.", flush = True)
            print("nreads\tlen_in_units\tcontents")
            print("\n".join([ f"{len(l)}\t{l[-1][1]-l[0][1]}\t{l}" for l in sorted(layouts, key=lambda x: len(x)) if len(l) > 1 ]))

            pickle.dump(layouts, open(f"layouts-{datetimestr}-noconflict-round-eov{eov_t}.pkl", "wb"))
        
        pickle.dump(layouts, open("layouts-{datetimestr}.pkl", "wb"))
        print("done")

    else:
        assert False, "invalid action."
