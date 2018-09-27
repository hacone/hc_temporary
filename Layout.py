# This is evolved from Alignments.py, which is then evolved from SNV_distribution.py
from scipy.stats import binom
from collections import Counter
import hashlib
from collections import namedtuple
import numpy as np
import pickle

from EncodedRead import *
from HOR_segregation import *

# Pool = namedtuple("Pool", ("hers", "arrs")) # simple class for the context; hor encoded reads, layouts, ...
# TODO: this is not what i wanted. for now, let us accept overriding (hers, arrs) wherever possible, falling back to default pool if it's not specified.

# a record for an alignment (for reads, and for layouts).
# eff_ovlp = effective length of ovlp, f_ext/r_ext = extention forward/reverse
Aln = namedtuple("Aln", ("i", "ai", "li", "j", "aj", "lj", "k", "len_ovlp", "eff_ovlp", "f_ext", "r_ext", "n00", "nmis", "n11", "score"))
LoAln = namedtuple("LoAln", ("l1", "l2", "k", "len_ovlp", "eff_ovlp", "n00", "nmis", "n11", "score")) # NOTE: do I have to calc extension stats?

# TODO: may I enrich this data? TODO: make this class...
# I need total length for calc overlap length
Layout = namedtuple("Layout", ("reads", "begin", "end"))
# class Layout():

T_agr = 700 # agree score threshold. alignment with score above this can be ...
T_gap = 100 # required score gap between best one vs 2nd best one.
T_dag = 600 # disagree score threshold. alignments with scores below this are considered false.

class Assembly:
    """ The assembly class holds all data you need to perform assembly.
        Currently, those would be HOR encoded reads, consecutive units in them, and tentative layouts. """

    def __init__(self, hers, arrs = None):

        self.hers = hers # HOR encoded reads

        # initialize self.arrs
        """ `arrs` is valid consecutive regions after filling up 22/23/24-mons or 33/34/35-mons gaps:
            the format is { ri : [[mi, mi, mi, ...], [mi, mi, ...]] } (mi = -1 for masked) """
        # TODO -0.1 for just masked, -n for n-monomer non--canonical units!
        if not arrs:
            self.arrs = self.get_consecutive_units(self.hers)
        else:
            self.arrs = arrs

        self.c_00 = 0.1
        self.c_ms = 1.5

    def get_consecutive_units(self, hers):
        regs = {}
        c_5u, c_22u, c_d1, c_d12, c_lo = 0, 0, 0, 0, 0

        for i, her in enumerate(hers):
            regs[i], ai, lh = [], -1, 0.5

            for h, s, t in her.hors:

                if t not in ["~", "@LO", "D39", "D28", "22U", "D1", "D12"]:
                    continue

                if t in ["@LO"] and h == lh: # special element to facilitate non-canonical units in layout!
                    regs[i][ai] += [-1]
                    c_lo += 1
                elif t in ["D39", "D28"] and h == lh:
                    regs[i][ai] += [-5]
                    c_5u += 1
                elif t == "D1" and h == lh:
                    regs[i][ai] += [-11]
                    c_d1 += 1
                elif t == "D12" and h == lh:
                    regs[i][ai] += [-10]
                    c_d12 += 1
                elif t == "22U" and h == lh:
                    regs[i][ai] += [-11,-11]
                    c_22u += 1
                elif t == "~" and h == lh:
                    regs[i][ai] += [h]
                elif t == "~" and h in [lh + 10, lh + 11, lh + 12]:
                    regs[i][ai] += [-0.1, h]
                elif t == "~" and h in [lh + 21, lh + 22, lh + 23]:
                    regs[i][ai] += [-0.1, -0.1, h]
                else:
                    regs[i] += [[h]]
                    ai += 1
                lh = h + s

            regs[i] = [ l for l in regs[i] if l ] # remove if empty

        print(f"Spanned variant units: 5u={c_5u}, 22u={c_22u}, d1={c_d1}, d12={c_d12}, and @LO = {c_lo}.")
        return { k : v for k, v in regs.items() if v } # remove if empty

    def get_variants(self, units, err_rate = 0.03, fq_upper_bound = 0.75):
        """ define SNVs over `units`. Currently, SNV is dict with following keys: k, p, b, f, c, binom_p """

        assert units, "trying to call variants over nothing. abort."

        counter = Counter()
        for ri, i in units:
            for j in range(2, 12): # skipping the first 2 monomers, where assignment is always wrong
                counter.update([ (j, s.pos, s.base) for s in self.hers[ri].mons[i+j].monomer.snvs ])

        # remove too frequent ones, and allow <1 false positives
        snvs = [ (k, p, b, c / len(units), c) for (k, p, b), c in counter.most_common() if c / len(units) < fq_upper_bound ]
        snvs = [ dict(k = k, p = p, b = b, f = c/len(units), c = c, binom_p = 1 - binom.cdf(c, len(units), err_rate)) for k, p, b, f, c in snvs ]
        snvs = [ s for s in snvs if s["binom_p"]*(171*10*3) < 1.0 ]
        return snvs

    def get_bits_dict(self, snvs, regs):
        """ construct bit vectors in dict """

        # NOTE: don't you think this is too operational?
        bits_dict = {}
        for i, ai, l in [ (i, ai, self.arrs[i][ai]) for i, ai in regs ]:
            v = []
            for h in l:
                if h < 0: # TODO: fix
                    v += [0] * len(snvs)
                else:
                    v += [ 1 if any([ (sp.pos, sp.base) == (s["p"], s["b"]) for sp in self.hers[i].mons[h+s["k"]].monomer.snvs ]) else -1 for s in snvs ]
            bits_dict[(i, ai)] = np.array(v).reshape(len(l), len(snvs))
        return bits_dict

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
        """ calculate some-vs-some alignments among reads(regions), returning dict of alns. """
        # TODO: i won't need all_vs_all_aln anymore?

        alns_dict = dict()
        n = 0

        for i, ai in targets:
            n += 1
            print(f"aligning {i} - {ai}. {n} / {len(targets)}")
            alns_dict[(i, ai)] = dict()
            for j, aj in [ (j, aj) for j, aj in queries if not (j, aj) == (i, ai)]:
                li, lj = bits_dict[(i, ai)].shape[0], bits_dict[(j, aj)].shape[0]
                # with at least 2 units in overlap
                alns = [ self.calc_align(i, ai, j, aj, k, bits_dict = bits_dict) for k in range(-lj+2, li-2) ]
                alns = [ aln for aln in alns if aln.score > T_dag - 100 and aln.eff_ovlp > 4 ] # NOTE: threshold depends on scoring scheme.
                if alns:
                    alns_dict[(i, ai)][(j, aj)] = sorted(alns, key = lambda x: -1 * x.score)
            l = f"{ len([ 0 for t, alns in alns_dict[(i, ai)].items() for aln in alns if aln.score > T_dag -100 and aln.eff_ovlp > 4 ]) } targets found for saving. "
            l += f"{ len([ 0 for t, alns in alns_dict[(i, ai)].items() for aln in alns if aln.score > T_dag and aln.eff_ovlp > 4 ]) } targets above T_dag = {T_dag}."
            print(l, flush = True)

        return alns_dict

    def get_all_vs_all_aln(self):

        units = [ (ri, h) for ri, er in enumerate(self.hers) for h, _, t in er.hors if t == "~" ]
        n_units = len(units) # NOTE: do I need this here?
        snvs = self.get_variants(units)

        regs = [ (i, ai) for i, a in self.arrs.items() for ai, l in enumerate(a) if len(l) > 6 ]
        bits_dict = self.get_bits_dict(snvs, regs)

        alns_dict = self.some_vs_some_alignment(regs, regs, bits_dict)

        # FIXME how I realize params consistency!
        return Alignment(alignments = alns_dict, ctx_hash = 0, snv_hash = 0, c_00 = self.c_00, c_ms = self.c_ms)

class Alignment:

    def __init__(self, alignments, ctx_hash, snv_hash, c_00, c_ms):
        # alignments :: alns_dict
        self.alignments = alignments
        self.ctx_hash = ctx_hash
        self.snv_hash = snv_hash
        self.c_00 = c_00
        self.c_ms = c_ms
        # c_00, c_ms = 0.4735, 1.5 # with good theory 
        # c_00, c_ms = 0.0, 1.0 # ignoring 00

    def __hash__(self):
        return 1234 # TODO


def iterative_layout(alns_dict, bits_dict, nodes, asm):
    """ for nodes and alns_dict, perform iterative layout process and returns list of obtained layout
        nodes :: [(i, ai)]
    """

    result_layouts = [] # this accumulates the results.
    n_iter = 0

    while len(nodes) > 1:

        n_iter += 1; print(f"\n{len(nodes)} nodes are to be layout. Layout No.{n_iter}.")

        edges = [ (i, ai, j, aj, alns_dict[(i, ai)][(j, aj)])
                for i, ai in nodes for j, aj in nodes
                if i < j and (j, aj) in alns_dict[(i, ai)] and alns_dict[(i, ai)][(j, aj)][0].score > T_dag ]

        # remove slippy pairs
        edges = [ e for e in edges if len(e[4]) < 2 or e[4][0].score - e[4][1].score > T_gap ] 

        edges = sorted(edges, key = lambda x: -x[4][0].score) # sorted with the best score for the pair

        if not edges:
            #print("No Good edges any more.")
            break
        else:
            #print(f"{len(edges)} edges available.") 
            #print("\n".join([ f"{x[4][0]}" for x in edges[:10] ]))
            pass


        # TODO: maybe check for internal consistency, ... then final iterative assembly procedure.
        # TODO: assert the size of edges

        # make a seed
        best = edges[0][4][0]
        seed_layout = Layout(
            reads = [ ((best.i, best.ai), 0), ((best.j, best.aj), best.k) ],
            begin = -best.r_ext, end = best.li + best.f_ext)
        edges = edges[1:]
        visited = [ (i, ai) for (i, ai), d in seed_layout.reads ]

        #print(f"\ncurrent layout has {len(seed_layout.reads)} reads:", flush = True) ; print(seed_layout)

        next_e = [ e for e in edges if bool((e[4][0].i, e[4][0].ai) in visited) ^ bool((e[4][0].j, e[4][0].aj) in visited) ]

        if next_e:
            #print("nexts:"); print("\n".join([ f"{x[4][0]}" for x in next_e[:5] ]))
            best = next_e[0][4][0]
        else:
            #print(f"\nNo initial extention: FINAL LAYOUT with {len(seed_layout.reads)} reads, {seed_layout.end - seed_layout.begin} units:")
            #print(seed_layout)
            result_layouts += [seed_layout]
            nodes -= set(visited)
            continue

        while True:

            if (best.i, best.ai) in visited:
                ll = best.lj
                lo_alns = [ asm.calc_align_layout(seed_layout,
                    Layout(reads = [((best.j, best.aj), 0)], begin = 0, end = best.lj), k, bits_dict)
                    for k in range(seed_layout.begin - best.lj + 2, seed_layout.end - 2 + 1) ]

            elif (best.j, best.aj) in visited:
                ll = best.li
                lo_alns = [ asm.calc_align_layout(seed_layout,
                    Layout(reads = [((best.i, best.ai), 0)], begin = 0, end = best.li), k, bits_dict)
                    for k in range(seed_layout.begin - best.li + 2, seed_layout.end - 2 + 1) ]

            lo_alns = sorted(lo_alns, key = lambda x: -x.score)
            lo_alns = [ la for la in lo_alns if la.eff_ovlp > 4 ]

            best_ext_aln = lo_alns[0]
            best_ext_i, best_ext_ai = best_ext_aln.l2.reads[0][0]

            #print(f"\na read {best_ext_aln.l2.reads[0][0][0]} aligned to the current layout.\nscr\t  k\teol")
            #print("\n".join([ f"{lo_aln.score}\t{lo_aln.k}\t{lo_aln.eff_ovlp}" for lo_aln in lo_alns[:3] ]))

            if best_ext_aln.score > T_dag and (best_ext_aln.score - lo_alns[1].score) > T_gap:
                # add new node into layout
                visited += [(best_ext_i, best_ext_ai)]
                seed_layout = Layout(
                        reads = seed_layout.reads + [((best_ext_i, best_ext_ai), best_ext_aln.k)],
                        begin = min(seed_layout.begin, best_ext_aln.k), end = max(seed_layout.end, best_ext_aln.k + ll))

                #next_e = [ e for e in edges if bool((e[4][0].i, e[4][0].ai) in visited) ^ bool((e[4][0].j, e[4][0].aj) in visited) ]

                # NOTE: this should be faster
                next_e = next_e + [ e for e in edges if 
                           (((e[4][0].i, e[4][0].ai) == (best_ext_i, best_ext_ai)) & ((e[4][0].j, e[4][0].aj) not in visited))
                         | (((e[4][0].j, e[4][0].aj) == (best_ext_i, best_ext_ai)) & ((e[4][0].i, e[4][0].ai) not in visited))]
                next_e = [ e for e in next_e if bool((e[4][0].i, e[4][0].ai) in visited) ^ bool((e[4][0].j, e[4][0].aj) in visited) ]

                #print(f"\ncurrent layout has {len(seed_layout.reads)} reads:", flush = True); print(seed_layout)
                # layout_consensus(seed_layout) # TODO: temporary?

                if next_e:
                    #print("nexts:"); print("\n".join([ f"{x[4][0]}" for x in next_e[:5] ]))
                    best = next_e[0][4][0]
                else:
                    break
            else:
                #print(f"\nBad alignment; current layout has {len(seed_layout.reads)} reads:", flush = True)# ; print(seed_layout)
                if len(next_e) > 1:
                    next_e = next_e[1:]
                    #print("nexts:"); print("\n".join([ f"{x[4][0]}" for x in next_e[:5] ]))
                    best = next_e[0][4][0]
                else:
                    break

        #print(f"\nNo additional extention: FINAL LAYOUT with {len(seed_layout.reads)} reads, {seed_layout.end - seed_layout.begin} units:")
        result_layouts += [seed_layout]
        #print(seed_layout)
        #print_layout(seed_layout, bits_dict)
        nodes -= set(visited)

        # layout_consensus(seed_layout) # TODO: temporary?
        seed_layout.describe(asm)
        print_layout(seed_layout, bits_dict)

    return result_layouts

def double_edge_component(alns_dict):
    """
    obtain so called slippy component, whose elements are connected by multi-edges.
    """

    def is_double(alns):
        return 1 < len([ aln for aln in alns if aln.score > T_dag and aln.len_ovlp > 4 ])

    nodes = set()
    for (i, ai), d in alns_dict.items():
        nodes |= { (i, ai) } | { (j, aj) for (j, aj), alns in d.items() }

    components = []
    while len(nodes) > 0:
        k, ak = nodes.pop()
        # print(f"\nchecking component for {k}-{ak} / {len(nodes)} remain to be checked.")

        curr_component = { (k, ak) } | { (j, aj) for (j, aj), alns in alns_dict[(k, ak)].items() if is_double(alns) }
        curr_checked = set()

        while len(curr_checked) < len(curr_component):
            l, al = (curr_component - curr_checked).pop()
            curr_checked |= { (l, al) }
            curr_component |= { (j, aj) for (j, aj), alns in alns_dict[(l, al)].items() if is_double(alns) }

        components += [[ e for e in curr_component ]]
        nodes -= curr_component
        # print(f"{len(curr_component)} nodes in the component")

    return components


class Layout:
    """ The layout is essentially a set of reads with their relative configurations specified.
        Reads are represented as ((rid, aid), rpos) tuple, so the layout must point to the context where it was defined. """

    def __init__(self, reads = None, begin = None, end = None):
        # TODO: check the consistency of ctx
        #self.ctx_hash = ctx_hash
        self.reads = reads
        self.begin = begin # NOTE: could be calculated on the fly?
        self.end = end
        self.consensus = None # lazy consensus

    def get_consensus(self, ctx): # TODO: make tidy here.
        """ construct a consensus read out of a layout / compute minimal-mismatches to evaluate layout
            returns a new HOR encoded read and an associated array. """

        min_variant_freq = 0.3 # freq of mismatches required to be considered as true variants

        #if self.consensus:
        #    return self.consensus

        lb, le = self.begin, self.end
        cons_snvs = dict() # NOTE: i can be negative, thus it must be dict.
        depth = dict()
        for i in range(lb, le): # NOTE: this is not really efficient, but hope it's not bottleneck.
            cons_snvs[i] = [ Counter() for ki in range(12) ]
            depth[i] = 0
            # i-th unit.
            for (ri, rai), k in self.reads:
                l = ctx.arrs[ri][rai]
                if k > i or k + len(l) <= i:
                    continue
                if l[i-k] < 0: # TODO fix
                    continue
                else: # if i[i-k] is canonical unit
                    for ki in range(2, 12): # NOTE: the first two mons should be ignored?
                        cons_snvs[i][ki].update([ (s.pos, s.base) for s in ctx.hers[ri].mons[l[i-k]+ki].monomer.snvs ])
                    depth[i] += 1

        nm = 0 # NOTE: for now, this excludes the first 2 (wrong) monomers
        for i in range(lb, le):
            for ki in range(2,12):
                for s, fq in cons_snvs[i][ki].items():
                    # FIXME: it allows >1 snvs at the same position !! 
                    nm += (depth[i] - fq) if (fq / depth[i] > min_variant_freq) else fq

        tot_units = sum(depth.values())
        len_cons = le - lb
        print(f"consensus taken: {len(self.reads)}\t{tot_units} => {len_cons}\t{nm}\t{nm/(tot_units-len_cons):.2f}")
        
        self.consensus = HOR_Read(
            name = "anything", length = (le-lb)*12*171, ori = "+",
            mons = [ AssignedMonomer(begin = 12*171*i + 171*ki, end = 12*171*i+171, ori = "+",
                monomer = Monomer(name = "*", snvs = [ SNV(pos = p, base = b) for (p, b), fq in cons_snvs[i+lb][ki].items() if fq/depth[i+lb] > min_variant_freq ])) # TODO: fix threshold
                for i in range(le-lb) for ki in range(12) ],
            hors = [ (i*12, 12, "~") if depth[lb+i] > 0 else (i*12, 12, "@LO") for i in range(le-lb) ])
        return self.consensus


    def define_local_variants(self, ctx):
        """ returns layout-wide locally-defined variants. """
        newpool = Assembly(hers = [ ctx.hers[i] for (i, ai), k in self.reads ])
        regs = [ (i, ai) for i, a in newpool.arrs.items() for ai, l in enumerate(a) if len(l) > 4 ]
        units_in_layout = [ (i, mi) for i, ai in regs for mi in newpool.arrs[i][ai] if mi >= 0 ]
        print(f"{len(units_in_layout)} units found in desc layout")
        return newpool.get_variants(units_in_layout)

    # TODO test
    def define_ends_variants(self, ctx, plusend = True):
        """ returns variants locally-defined at ends of the layout"""
        newpool = Assembly(hers = [ ctx.hers[i] for (i, ai), k in self.reads ])
        regs = [ (i, ai)
            for i, a in newpool.arrs.items()
            for ai, l in enumerate(a)
            if len(l) > 4 ]

        units_in_layout = [ ((i, mi), self.end - k + p)
            for (i, ai), k in self.reads
            for p, mi in enumerate(newpool.arrs[i][ai])
            if mi >= 0 ]
        units_in_end = sorted(units_in_layout,
            key = lambda x: (1 if plusend else -1) * x[1])

        print(f"{len(units_in_layout)} units found in desc layout")
        return newpool.get_variants(units_in_layout)

    def prune(self, ctx):

        """ remove any reads slippery w.r.t. consensus. """
        # rd = ((rid, aid), rpos)
        def leave_one_out_consensus(rd):
            reads_except_me = [ r for r in self.reads if r != rd ]
            b = min([ rpos for (rid, aid), rpos in reads_except_me ])
            e = max([ rpos + len(ctx.arrs[rid][aid]) for (rid, aid), rpos in reads_except_me ])
            loo_layout = Layout(reads_except_me, b, e)
            loo_consensus = loo_layout.get_consensus(ctx)
            return loo_consensus

        llsnvs = self.define_local_variants(ctx)

        def if_uniquely_align(ccs, read):
            np = Assembly(hers = [ccs] + ctx.hers[read[0]])
            np_regs = [ (i, ai) 
                for i, a in np.arrs.items()
                for ai, l in enumerate(a) if len(l) > 4 ]
            ccs_parts = len([ r for r in np_regs if r[0] == 0 ]) 
            if ccs_parts > 1:
                print(f"consensus is fluctured into {ccs_parts}.")
                return True # TODO: can I do that?
            # TODO: alignment actually.
            llbd = np.get_bits_dict(llsnvs, np_regs)
            i, ai = np_regs[0]
            j, aj = read[0]
            li = llbd[(i, ai)].shape[0]
            lj = llbd[(j, aj)].shape[0]
            alns = [ np.calc_align(i, ai, j, aj, k, bits_dict = llbd) for k in range(-lj-5, li+5) ]
            alns = [ aln for aln in sorted(alns, key = lambda x: -x.score) if aln.eff_ovlp > 4 ]
            return (aln[0].score - aln[1].score) > 15.0


        goods = [ rd for rd in self.reads
            if if_uniquely_align(leave_one_out_consensus(rd), rd) ]
        print(f"Pruned; {len(goods)} passed out of {len(self.reads)} reads.")

        # making layout!
        b = min([ rpos for (rid, aid), rpos in goods ])
        e = max([ rpos + len(ctx.arrs[rid][aid]) for (rid, aid), rpos in goods ])
        pruned_layout = Layout(goods, b, e)
        return pruned_layout


    def describe(self, ctx):
        """ describe layout's consistency (stability with local mask?). This should work fine for medium-size layout. """

        #newpool = Assembly(hers = [self.get_consensus(ctx)] + [ ctx.hers[i] for (i, ai), k in self.reads ])
        newpool = Assembly(hers = [self.get_consensus(ctx)] + [ ctx.hers[i] for (i, ai), k in sorted(self.reads, key=lambda x: x[1]) ])
        regs = [ (i, ai) for i, a in newpool.arrs.items() for ai, l in enumerate(a) if len(l) > 4 ]
        layout_local_snvs = self.define_local_variants(ctx)
        print_snvs(layout_local_snvs)

        layout_local_bits_dict = newpool.get_bits_dict(layout_local_snvs, regs)
        l0 = layout_local_bits_dict[regs[0]].shape[0] # layout's length
        print("The structure of the layout:")
        print(newpool.hers[0].hors)

        i, ai = regs[0]

        for j, aj in regs[1:]:
            lj = layout_local_bits_dict[(j, aj)].shape[0]
            k_in_l = self.reads[j-1][1] - self.begin # ?

            # with at least 2 units in overlap
            #alns = [ newpool.calc_align(i, ai, j, aj, k, bits_dict = layout_local_bits_dict) for k in range(-lj+2, l0-2) ]
            alns = [ newpool.calc_align(i, ai, j, aj, k, bits_dict = layout_local_bits_dict) for k in range(-lj-5, l0+5) ]
            alns = [ aln for aln in sorted(alns, key = lambda x: -x.score) if aln.eff_ovlp > 4 ]

            print("\n")
            print(f"- structure of the read {j}-{aj} ({lj} units):")
            print(f"{newpool.hers[j].hors}")
            print("- alignments found ... [scr]/[ovlp]@[pos]")
            print(f" ".join([f"*{aln.score/10:.1f}/{aln.eff_ovlp}@{aln.k}"
                if aln.k == k_in_l else f"{aln.score/10:.1f}/{aln.eff_ovlp}@{aln.k}" for aln in alns[:10] ]))

            if len(alns) > 1:
                if (alns[0].score - alns[1].score)/10 > 15:
                    print(f"\nScore Gap = {(alns[0].score - alns[1].score)/10:.1f} - PASS")
                else:
                    print(f"\nScore Gap = {(alns[0].score - alns[1].score)/10:.1f} - FAIL")
            elif alns:
                print(f"\nScore Gap = {(alns[0].score - 0):.1f} - PASS (ONLY_ALN)")
            if alns:
                print_align(alns[0], layout_local_bits_dict)

def print_snvs(snvs, alt_snvs = None):
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

def describe_alns_dict(alns_dict, bits_dict, alt_alns_dict = None, alt_bits_dict = None):
    """ visualize alns_dict in text """

    T_print = 650
    
    for (i, ai), d in alns_dict.items(): 
        targets = sorted([ (j, aj, alns) for (j, aj), alns in d.items() if alns[0].score > T_print ], key = lambda x: -x[2][0].score)
        if targets:
            print("\n--------------------------------------------------")
            n11, n00 = targets[0][2][0].n11, targets[0][2][0].n00
            print(f"Alignments for read {i} region {ai}...: {100*n11/(n11+n00):.2f} " + " ".join([f"{alns[0].score/10:.3f}" for j, aj, alns in targets[:10]]))
        for j, aj, alns in targets[:10]: # for the first 10 reads
            print(f"\nFrom read {i}, {ai} To read {j}, {aj}...")
            print("alt.confs: " + " ".join([ f"{s.score/10:.2f}~({s.k})" for s in alns[:10] ]))
            print_align(alns[0], bits_dict)

            if alt_alns_dict and alt_bits_dict and (i, ai) in alt_alns_dict and (j, aj) in alt_alns_dict[(i, ai)]:
                alt_alns = alt_alns_dict[(i, ai)][(j, aj)]
                print("\noriginal alt.confs: " + " ".join([ f"{s.score/10:.2f}~({s.k})" for s in alt_alns[:10] ]))
                alt_alns = [ alt_aln for alt_aln in alt_alns if alt_aln.k == alns[0].k ]
                if alt_alns:
                    print_align(alt_alns[0], alt_bits_dict)
                else:
                    print("\nNo corresponding alignment in original setting.")

def stats_alns_dict(alns_dict): # NOTE: I'm not using this, but leave it here for later reference
    """ empirical score distribution for plotting """
    # (score, gap, n11/n11+n00, eff_ovlp, rank)
    print("i\tj\tscore\tscoreGap\tvars_frac\teff_ovlp\trank")
    for (i, ai), d in alns_dict.items():
        for rank, (j, aj, alns) in enumerate(
                sorted([ (j, aj, alns) for (j, aj), alns in d.items() if (j, aj) != (i, ai) ], key = lambda x: -1.0 * x[2][0].score)[:10]):
            scoreGap = alns[0].score - alns[1].score if len(alns) > 1 else alns[0].score - (T_dag-100)
            vars_frac = alns[0].n11 / (alns[0].n11 + alns[0].n00) # n11/(n11+n00)
            line = f"{i}\t{j}\t{alns[0].score/10:.2f}\t{scoreGap/10:.2f}\t{100*vars_frac:.2f}\t{alns[0].eff_ovlp}\t{rank}"
            print(line)

# TODO: def plot_scoreGap_distr(), or any other equivalent
# TODO: some function to map shorter reads into layouts/components to enrich data.

def print_layout(layout, bits_dict):

    def is_in(i, j, aj, d):
        lj = bits_dict[(j, aj)].shape[0]
        if d <= i and i < d + lj:
            if bits_dict[(j, aj)][i-d,0] == 0:
                return "*"
            else:
                return "|"
        else:
            return " "

    print("\n   i:  i2\t>  >  >  reads  <  <  <")
    for i in range(layout.begin, layout.end):
        # NOTE: you want to sort by starting coordinate?
        #line = f"{i-layout.begin:<4}\t" + "".join([ "|" if is_in(i, j, aj, d) else " " for (j, aj), d in sorted(layout.reads, key=lambda x: x[1]) ])
        #line = f"{i:<4}:{i-layout.begin:<4}\t" + "".join([ is_in(i, j, aj, d) for (j, aj), d in layout.reads ])
        line = f"{i:<4}:{i-layout.begin:<4}\t" + "".join([ is_in(i, j, aj, d) for (j, aj), d in sorted(layout.reads, key=lambda x: x[1]) ])
        print(line)
    print("----:----\t" + "-" * len(layout.reads) + ";")

    return None

def layout(alns_dict = None, ctx = None, layouts = None):
    """ implementing layout idea """

    assert alns_dict, "needs Alignment."
    assert ctx, "needs Assembly as a context"

    bag_of_units = [ (ri, h) for ri, er in enumerate(ctx.hers) for h, _, t in er.hors if t == "~" ]
    n_units = len(bag_of_units)
    print(f"{n_units} units found in {len(ctx.hers)} reads.")

    snv_sites = ctx.get_variants(bag_of_units)
    #snv_sites = detect_snvs(bag_of_units)[:40] # force 40?
    print(f"{len(snv_sites)} SNV sites defined globally.")
    print_snvs(snv_sites)

    n_raw_hers = len(ctx.hers)

    # bits vectors are constructed for all regions, on SNVs defined globally.
    regs = [ (i, ai) for i, a in ctx.arrs.items() for ai, l in enumerate(a) if len(l) > 4 ]
    bits_dict = ctx.get_bits_dict(snv_sites, regs)

    # calculate alignments globally
    long_arrays = [ (i, ai) for i, a in ctx.arrs.items() for ai, l in enumerate(a) if len(l) > 6 ]

    #describe_alns_dict(alns_dict, bits_dict) 
    #stats_alns_dict(alns_dict)

    # NOTE: I'm trying to see what if layout with all nodes
    if not layouts:
        layouts = iterative_layout(alns_dict, bits_dict, set(long_arrays), ctx)
        with open("layouts.pickle", "wb") as f:
            pickle.dump(layouts, f)

    print(f"{len(layouts)} layouts obtained from all {len(long_arrays)} long (>6u) arrays")

    for il, l in enumerate(sorted(layouts, key = lambda x: -len(x.reads))):
        if len(l.reads) < 3:
            continue
        #print_layout(l, bits_dict)
        #l.describe(ctx)
        pass

    slippies = double_edge_component(alns_dict)
    print(f"{ len(slippies) } components found. CHECK1")
    print(f"Of these, { len([ s for s in slippies if len(s) > 1 ]) } comprises multiple reads (slippery components).")

    good_nodes = { comp[0] for comp in slippies if len(comp) == 1 } # nodes = [(i, ai)]
    nonslip_snv_sites = ctx.get_variants([ (i, h) for i, ai in good_nodes for h in ctx.arrs[i][ai] if h >= 0 ])
    nonslip_bits_dict = ctx.get_bits_dict(nonslip_snv_sites, regs)

    if True:
        nonslip_layouts = pickle.load(open("./non-slip-layouts.pickle", "rb"))
    else:
        nonslip_alns_dict = ctx.some_vs_some_alignment(good_nodes, good_nodes, nonslip_bits_dict)
        nonslip_layouts = iterative_layout(nonslip_alns_dict, nonslip_bits_dict, good_nodes, ctx)

    for il, l in enumerate(sorted(nonslip_layouts, key = lambda x: -len(x.reads))):
        if len(l.reads) < 3:
            continue

        print(f"a Layout no.{il} for non-slippery nodes")
        print_layout(l, nonslip_bits_dict)
        l.describe(ctx)

        continue

        # TODO : remove rest? already handled in the layout.describe
        cons, arr = layout_consensus(l)
        pool = Pool(hers = pool.hers + [cons], arrs = pool.arrs)
        pool.arrs.update({(n_raw_hers+il):arr})
        print(arr)
        print(f"{len(pool.hers)} reads, {len(pool.arrs)} arrays.  k = {n_raw_hers + il}")

        # using layout-local snv and recalculate alignment dict.
        comp = [ (i, ai) for (i, ai), k in l.reads ]
        units_in_comp = [ (i, mi) for (i, ai), k in l.reads for mi in arrs[i][ai] if mi >= 0 ]
        snvs = detect_snvs(units_in_comp) # local
        print(f"\n{len(snvs)} SNVs for the layout {il}: {len(l.reads)} reads; {len(units_in_comp)} units.")
        print_snvs(snvs, snv_sites)

        bits_dict_local = get_bits_dict(snvs, comp) # local
        alns_dict_original = all_vs_all_aln(comp, bits_dict)
        alns_dict_local = all_vs_all_aln(comp, bits_dict_local)
        describe_alns_dict(alns_dict_local, bits_dict_local, alns_dict, bits_dict)

        # print_layout for old and updated?

    print(f"{len(nonslip_layouts)} layouts obtained for non-slippery parts: " +\
          f"{sum([ len(l.reads) for l in nonslip_layouts ])} nodes, "+\
          f"{sum([ l.end - l.begin for l in nonslip_layouts ])} units long.")

    if True:
        pass
    else:
        with open("non-slip-layouts.pickle", "wb") as f:
            pickle.dump(nonslip_layouts, f)

    for ic, comp in enumerate(slippies): # NOTE: suppressed just for the experiment below.

        if len(comp) < 2:
            continue

        units_in_comp = [ (i, mi) for i, ai in comp for mi in ctx.arrs[i][ai] if mi >= 0 ]
        slip_snv_site = ctx.get_variants(units_in_comp)
        slip_bits_dict = ctx.get_bits_dict(slip_snv_site, regs) # NOTE: a bit superflous to calc dict for all regs
        slip_alns_dict = ctx.some_vs_some_alignment(comp, comp, slip_bits_dict)
        slip_layouts = iterative_layout(slip_alns_dict, slip_bits_dict, set(comp), ctx)

        #snvs = detect_snvs(units_in_comp) # local
        print(f"\n{len(slip_snv_site)} SNVs for the component {ic}: {len(comp)} reads; {len(units_in_comp)} units.")
        #print_snvs(snvs, snv_sites)
        #bits_dict_local = get_bits_dict(snvs, comp) # local

        # alignment locally, using global SNVs 
        #alns_dict_original = all_vs_all_aln(comp, bits_dict)
        #print("alns in original SNVs")
        #stats_alns_dict(alns_dict_original)

        # alignment locally in the component, using new SNVs
        #alns_dict_local = all_vs_all_aln(comp, bits_dict_local)
        #print("alns in updated SNVs")
        #stats_alns_dict(alns_dict_local)

        """
        describe_alns_dict(alns_dict_local, bits_dict_local, alns_dict, bits_dict)
        """

        #layouts = iterative_layout(alns_dict, bits_dict, set(comp))
        print(f"\n\nUsing SNVs, {len(slip_layouts)} layouts obtained for slippery component No.{ic}: " +\
              f"{sum([ len(l.reads) for l in slip_layouts ])} (/ {len(comp)}) nodes, "+\
              f"{sum([ l.end - l.begin for l in slip_layouts ])} (/ {len(units_in_comp)}) units long.")
 
        for il, l in enumerate(sorted(slip_layouts, key = lambda x: -len(x.reads))):
            if len(l.reads) < 3:
                continue
            print(f"a Layout no.{il} for a slippery component no.{ic}")
            print_layout(l, slip_bits_dict)
            l.describe(ctx)

        pickle.dump(slip_layouts, open(f"./slip-layouts/slip-layouts-{ic}.pickle", "wb"))

        continue

        import sys
        sys.exit()

        layouts = iterative_layout(alns_dict_local, bits_dict_local, set(comp))
        print(f"Using  updated SNVs, {len(layouts)} layouts obtained for slippery component No.{ic}: " +\
              f"{sum([ len(l.reads) for l in layouts ])} (/ {len(comp)}) nodes, "+\
              f"{sum([ l.end - l.begin for l in layouts ])} (/ {len(units_in_comp)}) units long.")

        for l in layouts:
            layout_consensus(l)
            describe_layout(l)

        for l in layouts:
            cons, arr = layout_consensus(l)
            pool.arrs.update({len(pool.hers):arr})
            pool = Pool(hers = pool.hers + [cons], arrs = pool.arrs)
            print(arr)
            print(f"{len(pool.hers)} reads, {len(pool.arrs)} arrays.  k = {n_raw_hers + il}")

    import sys
    sys.exit()

    # all-layout v all-layout alignment
    regs = [ (i, ai) for i, a in pool.arrs.items() for ai, l in enumerate(a) if len(l) > 4 and i >= n_raw_hers ]
    units_in_comp = [ (i, mi) for i, ai in regs for mi in arrs[i][ai] if mi >= 0 ]
    snvs = detect_snvs(units_in_comp) # local
    print(f"\n{len(snvs)} SNVs for LvL: {len(regs)} layouts; {len(units_in_comp)} units.")
    print_snvs(snvs, snv_sites)

    bits_dict_local = get_bits_dict(snvs, regs) # local
    #bits_dict = get_bits_dict(snv_sites, regs)
    layouts_alns_dict = all_vs_all_aln(regs, bits_dict_local)
    describe_alns_dict(layouts_alns_dict, bits_dict_local)
    
    layout_of_layouts = iterative_layout(layouts_alns_dict, bits_dict_local, set(regs))
    for l in layout_of_layouts:
        print_layout(l, bits_dict_local)

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='perform pair-wise alignment among HOR encoded reads on SNV data, again') # TODO: explain
    parser.add_argument('action', metavar='action', type=str, help='action to perform: layout, ...')
    parser.add_argument('--hor-reads', dest='hors', help='pickled hor-encoded long reads')
    parser.add_argument('--alns', dest='alns', help='pickled alignments')
    parser.add_argument('--layouts', dest='layouts', help='pickled layouts; it\'s your job to take care of consistency with alns...')
    args = parser.parse_args()

    # global pool # NOTE: to be deprecated soon

    if args.action == "align":

        assert args.hors, "need HOR-encoded reads"

        hers = pickle.load(open(args.hors, "rb"))
        assembly = Assembly(hers)
        alns = assembly.get_all_vs_all_aln()
        with open("calign.pickle", "wb") as f:
            pickle.dump(alns, f)

    elif args.action == "layout":

        assert args.hors, "need HOR-encoded reads"
        assert args.alns, "need alignments"
        hers = pickle.load(open(args.hors, "rb"))
        assembly = Assembly(hers)
        alns = pickle.load(open(args.alns, "rb"))

        # TODO: check consistency here
        if args.layouts:
            layouts = pickle.load(open(args.layouts, "rb"))
            layout(alns.alignments, ctx = assembly, layouts = layouts)
        else:
            layout(alns.alignments, ctx = assembly)

    else:
        assert None, "invalid action."
