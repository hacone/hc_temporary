# Perform alignment among hor encoded reads (from X array).
from scipy.stats import binom
from collections import Counter
import hashlib
from collections import namedtuple
import numpy as np
import pickle
from underscore import _ as us

# TODO: I just need datatype definition. That can be separated from other codes.
from EncodedRead import *
from HOR_segregation import *

# a record for an alignment (for reads, and for layouts).
# eff_ovlp = effective length of ovlp, f_ext/r_ext = extention forward/reverse
Aln = namedtuple("Aln", ("i", "ai", "li", "j", "aj", "lj", "k", "len_ovlp", "eff_ovlp", "f_ext", "r_ext", "n00", "nmis", "n11", "score"))
LoAln = namedtuple("LoAln", ("l1", "l2", "k", "len_ovlp", "eff_ovlp", "n00", "nmis", "n11", "score")) # NOTE: do I have to calc extension stats?

# I need total length for calc overlap length
Layout = namedtuple("Layout", ("reads", "begin", "end"))

T_agr = 700 # agree score threshold. alignment with score above this can be ...
T_gap = 100 # required score gap between best one vs 2nd best one.
T_dag = 600 # disagree score threshold. alignments with scores below this are considered false.

# TODO: can this be abstracted to work with other chromosomes?
def chopx(reads):
    """
        chop HOR-encoded reads from X, return `{ rid : [[mids, ...], ...] }`
        Specifically, they're valid consecutive regions after 22/23/24-mons or 33/34/35-mons gaps are filled.
        the format is { ri : [[mi, mi, mi, ...], [mi, mi, ...]] } (mi < 0 for masked)
    """
    regs = {}

    def _chopx(hors):
        """ chop one read """
        ret, ai, lh = [], -1, 0.5
        for h, s, t in hors:
            if t not in ["~", "@LO", "D39", "D28", "22U", "D1", "D12"]:
                continue

            if t in ["@LO"] and h == lh: # special element to facilitate non-canonical units in layout!
                ret[ai] += [-1]
            elif t == "D39" and h == lh:
                ret[ai] += [-2]
            elif t == "D28" and h == lh:
                ret[ai] += [-3]
            elif t == "D1" and h == lh:
                ret[ai] += [-4]
            elif t == "D12" and h == lh:
                ret[ai] += [-5]
            elif t == "22U" and h == lh:
                ret[ai] += [-6,-6]

            elif t == "~" and h == lh:
                ret[ai] += [h] # NOTE: positivity matters
            elif t == "~" and h in [lh + 10, lh + 11, lh + 12]:
                ret[ai] += [-0.5, h]
            elif t == "~" and h in [lh + 21, lh + 22, lh + 23]:
                ret[ai] += [-0.5, -0.5, h]
            else:
                ret += [[h]] if t == "~" else [[-0.5]] # TODO logic!!
                ai += 1
            lh = h + s

        # return with removing empty items
        return [ l for l in ret if l ]

    # return with removing empty items
    regs = { i : _chopx(her.hors) for i, her in enumerate(reads) }
    return { k : v for k, v in regs.items() if v }

def var(reads, units = None, err_rate = 0.03, fq_upper_bound = 0.75):
    """ Define SNVs over `units` in `reads`. Currently, SNV is dict with following keys: k, p, b, f, c, binom_p
        If `units` is not specified, it uses all default "~" units. """

    counter = Counter()
    if not units:
        units = [ (ri, h) for ri, er in enumerate(reads) for h, _, t in er.hors if t == "~" ]

    for ri, i in units:
        for j in range(2, 12): # NOTE: specialied for X. skipping the first 2 monomers, where assignment is always wrong
            counter.update([ (j, s.pos, s.base) for s in reads[ri].mons[i+j].monomer.snvs ])

    # remove too frequent ones, and allow <1 false positives
    _tmp = [ (k, p, b, c / len(units), c) for (k, p, b), c in counter.most_common() if c / len(units) < fq_upper_bound ]
    _nt_vars = [ dict(k = k, p = p, b = b, f = c/len(units), c = c, binom_p = 1 - binom.cdf(c, len(units), err_rate)) for k, p, b, f, c in _tmp ]
    return [ s for s in _nt_vars if s["binom_p"]*(171*10*3) < 1.0 ]

class Alignment: # TODO: rename this!!
    """
        The Alignment class, initialized with a set of HOR encoded reads,
        holds relevant data such as induced (chopped) array, default global variants, matrix representation,
        and other utility methods to perform alignment (all-vs-all or specific).
    """

    def __init__(self, hers, arrs = None, variants = None):
        self.hers = hers # HOR encoded reads
        self.arrs = chopx(self.hers) if not arrs else arrs
        self.variants = var(hers) if not variants else variants

        # NOTE: self.bits can be reconstructed from vars and arrs, so depend on them essentially.
        self.bits = self.get_bits_dict(self.variants, self.longer_than(0))

        # NOTE: I won't change them never?
        self.c_00 = 0.1
        self.c_ms = 1.5

    def longer_than(self, tl):
        """ helper function to return region indices for later use. """
        return [ (i, ai) for i, a in self.arrs.items() for ai, l in enumerate(a) if len(l) > tl ]

    def get_bits_dict(self, snvs, regs):
        """ construct bit array representation as dict """

        reads, arrs = self.hers, self.arrs
        # regs = self.longer_than(0)

        def vec(i, h):
            """ bit array for a single unit, starting at h-th monomers of read i """
            def has_s(s):
                # i.e., if s in t
                #return any([ (t.pos, t.base) == (s["p"], s["b"]) for t in reads[i].mons[h+s["k"]].monomer.snvs ])
                return any([ (t.pos, t.base) == (int(s["p"]), s["b"]) for t in reads[i].mons[h+int(s["k"])].monomer.snvs ])
            return [ 1 if has_s(s) else -1 for s in snvs ]

        def _ba(i, ai):
            """ bit array for a single array """
            l = arrs[i][ai]
            v = [ [0] * len(snvs) if h < 0 else vec(i, h) for h in l ]
            return np.array(v).reshape(len(l), len(snvs))

        return { (i, ai) : _ba(i, ai) for i, ai in regs } 

    def mismatX(self, i, ai, j, aj, k, bits_dict):
        """ aux func returning stats of matches and mismatches
            eventually, this could be a mismatch matrix X
        """

        li, lj = bits_dict[(i, ai)].shape[0], bits_dict[(j, aj)].shape[0]
        if k < 0:
            len_ovlp = min(li, k+lj)
            if len_ovlp < 1:
                return dict(li = li, lj = lj, len_ovlp = 0, eff_ovlp = 0, n11 = 0, n00 = 0, match = 0, mismatch = 0)
            xi = bits_dict[(i, ai)][0:len_ovlp,:]
            xj = bits_dict[(j, aj)][-k:-k+len_ovlp,:]
        else:
            len_ovlp = min(lj, -k+li) 
            if len_ovlp < 1:
                return dict(li = li, lj = lj, len_ovlp = 0, eff_ovlp = 0, n11 = 0, n00 = 0, match = 0, mismatch = 0)
            xi = bits_dict[(i, ai)][k:k+len_ovlp,:]
            xj = bits_dict[(j, aj)][0:len_ovlp,]

        m = np.multiply(xi, xj) # 1: match, 0: masked, -1: mismatch
        match = np.multiply((m + 1), m) / 2 # 1: match, 0: masked/mismatch
        mismatch = np.sum(np.multiply(m, m) - match) # 1: mismatch, 0: masked/match # n01 or n10
        n11 = np.sum(np.multiply(match, (xi + 1) / 2)) # n11
        n00 = np.sum(match) - n11
        eff_ovlp = len_ovlp - sum([ 1 if x == 0 else 0 for x in m[:,0] ])

        return dict(li = li, lj = lj, len_ovlp = len_ovlp, eff_ovlp = eff_ovlp,
                n11 = n11, n00 = n00, match = np.sum(match), mismatch = mismatch)

    def calc_align(self, i, ai, j, aj, k, bits_dict):
        """ calculate alignment for a pair; in a specified configuration """

        X = self.mismatX(i, ai, j, aj, k, bits_dict = bits_dict)
        score = 1000 * (X["n11"] + self.c_00 * X["n00"]) / (0.01 + X["n11"] + self.c_00 * X["n00"] + self.c_ms * X["mismatch"]) # per-mille !

        return Aln(i = i, ai = ai, li = X["li"], j = j, aj = aj, lj = X["lj"], k = k,
            score = score, len_ovlp = X["len_ovlp"], eff_ovlp = X["eff_ovlp"],
            f_ext = max(0, k+X["lj"]-X["li"]), r_ext = max(0, -k),
            n00 = X["n00"], nmis = X["mismatch"], n11 = X["n11"])

    def calc_align_layout(self, l1, l2, k, bits_dict):
        """ calculate alignment between two layouts; in a specified (by displacement k) configuration """

        t_ovlp, t_eff_ovlp = 0, 0
        t_n11, t_n00, t_mat, t_mis = 0, 0, 0, 0

        for (i, ai), di in l1.reads:
            for (j, aj), dj in l2.reads:
                X = self.mismatX(i, ai, j, aj, dj-di+k, bits_dict = bits_dict)
                t_ovlp += X["len_ovlp"]; t_eff_ovlp += X["eff_ovlp"]
                t_n11 += X["n11"]
                t_n00 += X["n00"]
                t_mis += X["mismatch"]

        score = 1000 * (t_n11 + self.c_00 * t_n00) / (0.01 + t_n11 + self.c_00 * t_n00 + self.c_ms * t_mis)

        return LoAln(l1 = l1, l2 = l2, k = k, len_ovlp = t_ovlp, eff_ovlp = t_eff_ovlp,
                n00 = t_n00, nmis = t_mis, n11 = t_n11, score = score)

    def some_vs_some_alignment(self, targets, queries, bits_dict): # TODO: bits_dict might be adaptive (local w.r.t. targets?)
        """ calculate some-vs-some alignments among reads(regions), returning dict of alns.
            also you have more control over the parameters (variants to be used) thru bits """
        # TODO: i won't need all_vs_all_aln anymore?

        alns_dict = dict()
        n = 0

        print(f"aligning for {len(targets)}")
        print("[P]", end = "", flush = True)
        for i, ai in targets:
            print(".", end = "", flush = True)
            #n += 1
            #print(f"aligning {i} - {ai}. {n} / {len(targets)}")
            alns_dict[(i, ai)] = dict()
            for j, aj in [ (j, aj) for j, aj in queries if not (j, aj) == (i, ai)]:
                li, lj = bits_dict[(i, ai)].shape[0], bits_dict[(j, aj)].shape[0]
                # with at least 2 units in overlap
                alns = [ self.calc_align(i, ai, j, aj, k, bits_dict = bits_dict) for k in range(-lj+2, li-2) ]
                alns = [ aln for aln in alns if aln.score > T_dag - 100 and aln.eff_ovlp > 4 ] # NOTE: threshold depends on scoring scheme.
                if alns:
                    alns_dict[(i, ai)][(j, aj)] = sorted(alns, key = lambda x: -1 * x.score)
            #l = f"{ len([ 0 for t, alns in alns_dict[(i, ai)].items() for aln in alns if aln.score > T_dag -100 and aln.eff_ovlp > 4 ]) } targets found for saving. "
            #l += f"{ len([ 0 for t, alns in alns_dict[(i, ai)].items() for aln in alns if aln.score > T_dag and aln.eff_ovlp > 4 ]) } targets above T_dag = {T_dag}."
            #print(l, flush = True)

        return alns_dict

    def get_all_vs_all_aln(self, regs = None, variants = None):
        """ calculate all-vs-all alignments among reads(regions), returning dict of alns.
            parameters are default, induced from the contextual object. (cf. some_vs_some_alignment) """

        regs = regs if regs else self.longer_than(6)

        if variants:
            bits = self.get_bits_dict(variants, regs)
        else: 
            variants = self.variants
            bits = self.bits

        # NOTE: following 2 are equivalent
        # regs = [ (i, ai) for i, a in self.arrs.items() for ai, l in enumerate(a) if len(l) > 6 ]

        return AlignmentStore(
            reads = self.hers,
            arrays = self.arrs,
            variants = variants,
            alignments = self.some_vs_some_alignment(regs, regs, bits),
            c_00 = self.c_00,
            c_ms = self.c_ms)

def print_snvs(snvs, alt_snvs = None):
    """ Show SNVs """
    lines = "  k\t  p\tb\t    c\t   f\t   p-val\t c.f.d.\t f.gl\n"
    for s in sorted(snvs, key = lambda x: -x["c"]):
        k, p, b, c, f, pv = s["k"], s["p"], s["b"], s["c"], s["f"], s["binom_p"]
        if alt_snvs:
            alt_f = [ a["f"] for a in alt_snvs if (k, p, b) == (a["k"], a["p"], a["b"]) ]
            alt_info = f"{100*alt_f[0]:>5.2f}" if alt_f else "    *"
            lines += f"{k:>3}\t{p:>3}\t{b}\t{c:>5}\t{100*f:.2f}\t{pv:.2e}\t{pv*171*10*3:>7.2f}\t{alt_info}\n"
        else:
            lines += f"{k:>3}\t{p:>3}\t{b}\t{c:>5}\t{100*f:.2f}\t{pv:.2e}\t{pv*171*10*3:>7.2f}\n"
    print(lines)

def print_align(aln, bits_dict):
    """ Visualize an alignment in text to clarify the overlap. """
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

    n_snvs = r.shape[1]
    _l = f"== {n_snvs} SNVs =="
    lines = ("  i\t{0:^" + f"{n_snvs}" + "}\tj  ").format(_l)

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

class AlignmentStore:
    """
    storing the results of all-vs-all alignments with parameters being used,
    which include: reads, arrays, variants. 
    NOTE: usually, arrays doesn't change for a set of reads (but it could).
    NOTE: bit matrix representation is (must be) fully determined by these, thus is not stored explicitly.
    """
    def __init__(self, reads, arrays, variants, alignments, c_00, c_ms):
        self.reads = reads
        self.arrays = arrays
        self.variants = variants
        self.alignments = alignments
        self.c_00 = c_00
        self.c_ms = c_ms

    def describe(self, bits_dict, alt_alns_dict = None, alt_bits_dict = None):
        """ visualize alns_dict in text """

        T_print = 650

        for (i, ai), d in self.alignments.items():
            targets = sorted([ (j, aj, alns) for (j, aj), alns in d.items() if alns[0].score > T_print ], key = lambda x: -x[2][0].score)
            if targets:
                print("\n--------------------------------------------------")
                n11, n00 = targets[0][2][0].n11, targets[0][2][0].n00
                print(f"Alignments for read {i} region {ai}...: {100*n11/(n11+n00):.2f} " + " ".join([f"{alns[0].score/10:.3f}" for j, aj, alns in targets[:10]]))

            for j, aj, alns in targets[:10]: # for the first 10 reads
                print(f"\nFrom read {i}, {ai} To read {j}, {aj}...")
                print("alt.confs: " + " ".join([ f"{s.score/10:.2f}~({s.k})" for s in alns[:10] ]))
                print_align(alns[0], bits_dict)

                """
                if alt_alns_dict and alt_bits_dict and (i, ai) in alt_alns_dict and (j, aj) in alt_alns_dict[(i, ai)]:
                    alt_alns = alt_alns_dict[(i, ai)][(j, aj)]
                    print("\noriginal alt.confs: " + " ".join([ f"{s.score/10:.2f}~({s.k})" for s in alt_alns[:10] ]))
                    alt_alns = [ alt_aln for alt_aln in alt_alns if alt_aln.k == alns[0].k ]
                    if alt_alns:
                        print_align(alt_alns[0], alt_bits_dict)
                    else:
                        print("\nNo corresponding alignment in original setting.")
                """

def stats_alns_dict(alns_dict): # NOTE: I'm not using this, but leave it here for later reference
    """ empirical score distribution for plotting (deprecated?) """
    # (score, gap, n11/n11+n00, eff_ovlp, rank)
    print("i\tj\tscore\tscoreGap\tvars_frac\teff_ovlp\trank")
    for (i, ai), d in alns_dict.items():
        for rank, (j, aj, alns) in enumerate(
                sorted([ (j, aj, alns) for (j, aj), alns in d.items() if (j, aj) != (i, ai) ], key = lambda x: -1.0 * x[2][0].score)[:10]):
            scoreGap = alns[0].score - alns[1].score if len(alns) > 1 else alns[0].score - (T_dag-100)
            vars_frac = alns[0].n11 / (alns[0].n11 + alns[0].n00) # n11/(n11+n00)
            line = f"{i}\t{j}\t{alns[0].score/10:.2f}\t{scoreGap/10:.2f}\t{100*vars_frac:.2f}\t{alns[0].eff_ovlp}\t{rank}"
            print(line)

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='perform pair-wise alignment among HOR encoded reads on SNV data, again')
    parser.add_argument('action', metavar='action', type=str, help='action to perform: align, ...')
    parser.add_argument('--hor-reads', dest='hors', help='pickled hor-encoded long reads (required for align)')
    parser.add_argument('--alignments', dest='alns', help='pickled AlignmentStore (required for print)')
    parser.add_argument('-o', dest='outfile', help='the file to be output (required for align)')
    args = parser.parse_args()

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

    else:
        assert False, "invalid action."
