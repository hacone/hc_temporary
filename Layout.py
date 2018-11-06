# Generates and refines layouts of hor encoded reads.
from Alignment import *

from scipy.stats import binom
from collections import Counter
from collections import namedtuple
import numpy as np
import pickle
from underscore import _ as us
from EncodedRead import *
from HOR_segregation import *
from Alignment import *

# I need total length for calc overlap length
Layout = namedtuple("Layout", ("reads", "begin", "end"))

T_agr = 700 # agree score threshold. alignment with score above this can be ...
T_gap = 100 # required score gap between best one vs 2nd best one.
T_dag = 600 # disagree score threshold. alignments with scores below this are considered false.

def read_cover(alignmentstore, radius = 20, thres = 650):
    """ overlapping cover of reads to define locally valid SNVs.
        returns :: [[reg]] """

    ast = alignmentstore # alias
    reads = ast.reads
    all_nodes = set(ast.alignments.keys())
    remaining_nodes = set(ast.alignments.keys())

    def one_cover(centroid, nodes):
        """ generage a possible cover. """

        assert len(nodes) > 0, "calculating a cover for an empty set of nodes"

        next_n = set([ (centroid, 0) ])
        seen = set()

        while len(next_n) > 0:
            reg, k = next_n.pop()
            seen |= set([reg])
            alns = sorted(ast.alignments[reg].items(), key = lambda x: -x[1][0].score)[:5]

            # print(f"seen = {len(seen)} : next_n = {len(next_n)} : k = {k}")

            next_n |= set([ ((j, aj), k + aln[0].k) for (j, aj), aln in alns \
                if (j, aj) not in seen and aln[0].score > thres and k + aln[0].k < radius and k + aln[0].k > -radius ])

        return seen

    result_covers = []
    print(f"\n{len(remaining_nodes)} nodes remaining...") 
    while len(remaining_nodes) > 1:
        cent = remaining_nodes.pop()
        s = one_cover(cent, all_nodes - set([cent]))
        if not s:
            break
        remaining_nodes -= s
        if len(s) > 0:
            result_covers += [s]
            print(f"\n{cent} covered with {len(s)} nodes; {len(remaining_nodes)} nodes remaining...")
        else:
            print(".", end="", flush=True)

    print(f"{len(result_covers)} good covers found.")
    return result_covers

# TODO: I'll get AlignmentStore as an input. check all the caller.
# TODO: make this readable.
def iterative_layout(alignmentstore, nodes = None):
    """
        Perform iterative layout process and returns list of obtained layout.
        Alignments and other contexts are given as `alignmentstore`. Only `nodes` specified are used.
        nodes :: [(i, ai)]
    """

    # TODO: this is for backwards compat.
    ast = Alignment(alignmentstore.reads, alignmentstore.arrays, alignmentstore.variants)
    alns_dict = alignmentstore.alignments
    bits_dict = ast.bits
    remaining_nodes = nodes if nodes else alignmentstore.alignments.keys()

    result_layouts = [] # this accumulates the results.

    def one_layout(nodes):
        """ generage a possible layout. """

        def valid_edges(nodes):
            """ return non-slippy, valid edges (alignments) among nodes. """
            def is_aligned(d, i, ai, j, aj, threshold):
                return (j, aj) in d[(i, ai)] and d[(i, ai)][(j, aj)][0].score > threshold
            # For removing slippy pairs
            def is_ambiguous(d, i, ai, j, aj, threshold):
                alns = d[(i, ai)][(j, aj)]
                return len(alns) > 1 and (alns[0].score -alns[1].score) < threshold

            edges = [
                (i, ai, j, aj, alns_dict[(i, ai)][(j, aj)]) for i, ai in nodes for j, aj in nodes
                if i < j and is_aligned(alns_dict, i, ai, j, aj, T_dag)
                and not is_ambiguous(alns_dict, i, ai, j, aj, T_gap) ]
            # sorted with the best score for the pair
            return sorted(edges, key = lambda x: -x[4][0].score)

        edges = valid_edges(nodes)

        if not edges:
            return None

        # make a seed
        best = edges[0][4][0]
        seed_layout = Layout(
            reads = [ ((best.i, best.ai), 0), ((best.j, best.aj), best.k) ],
            begin = -best.r_ext, end = best.li + best.f_ext)
        edges = edges[1:]
        visited = [ (i, ai) for (i, ai), d in seed_layout.reads ]

        #print(f"\ncurrent layout has {len(seed_layout.reads)} reads:", flush = True) ; print(seed_layout)
        next_e = [ e for e in edges if bool((e[4][0].i, e[4][0].ai) in visited) ^ bool((e[4][0].j, e[4][0].aj) in visited) ]

        #print("nexts:"); print("\n".join([ f"{x[4][0]}" for x in next_e[:5] ]))

        while next_e:

            def align_to_current_layout(best, origin):
                """ just returns the set of alignments from best edges to the current layout collectively. """

                if (best.i, best.ai) in origin:
                    i, ai, li = best.j, best.aj, best.lj
                else:
                    i, ai, li = best.i, best.ai, best.li
                target = Layout(reads = [((i, ai), 0)], begin = 0, end = li)
                return us([ k for k in range(seed_layout.begin - li + 2, seed_layout.end - 2 + 1) ]).chain() \
                    .map(lambda k, *a: ast.calc_align_layout(seed_layout, target, k, ast.bits)) \
                    .filter(lambda al, *a: al.eff_ovlp > 4) \
                    .sortBy(lambda al, *a: -al.score).value()

            def admittable_pairwise_align(alns):
                if len(alns) == 1:
                    return True # TODO: can I do this?
                return alns[0].score > T_dag and (alns[0].score - alns[1].score) > T_gap

            best = next_e[0][4][0]
            alns = align_to_current_layout(best, visited)

            if not admittable_pairwise_align(alns):
                next_e = next_e[1:]
                continue

            if (best.i, best.ai) in visited:
                j, aj, lj = best.j, best.aj, best.lj
            else:
                j, aj, lj = best.i, best.ai, best.li
            print(".", end="", flush=True)

            # add new node into layout
            visited += [(j, aj)]
            seed_layout = Layout(
                    reads = seed_layout.reads + [((j, aj), alns[0].k)],
                    begin = min(seed_layout.begin, alns[0].k), end = max(seed_layout.end, alns[0].k + lj)) # TODO

            # update the set of extending edges
            # this is a slower implementation
            #next_e = [ e for e in edges if bool((e[4][0].i, e[4][0].ai) in visited) ^ bool((e[4][0].j, e[4][0].aj) in visited) ]
            next_e = next_e + [ e for e in edges if 
                       (((e[4][0].i, e[4][0].ai) == (j, aj)) & ((e[4][0].j, e[4][0].aj) not in visited))
                     | (((e[4][0].j, e[4][0].aj) == (j, aj)) & ((e[4][0].i, e[4][0].ai) not in visited))]
            next_e = [ e for e in next_e if bool((e[4][0].i, e[4][0].ai) in visited) ^ bool((e[4][0].j, e[4][0].aj) in visited) ]

            #print(f"\ncurrent layout has {len(seed_layout.reads)} reads:", flush = True); print(seed_layout)
            # layout_consensus(seed_layout) # TODO: temporary?
            #print("nexts:"); print("\n".join([ f"{x[4][0]}" for x in next_e[:5] ]))

        return dict(layout = seed_layout, remaining_nodes = nodes - set(visited))

    n = 0
    while len(remaining_nodes) > 1:
        print(f"\n{len(remaining_nodes)} nodes remaining...") 
        r = one_layout(remaining_nodes)
        if not r:
            break
        result_layouts += [ r["layout"] ]
        remaining_nodes = r["remaining_nodes"]

        if len(r["layout"].reads) > 2:
            print("\n", end="")
            print_layout(r["layout"], bits_dict)
        print(f"{sum([ len(l.reads) for l in result_layouts ])} nodes in layout.")
        print(f"{sum([ l.end - l.begin for l in result_layouts ])} units long in total.")


    return result_layouts

        #print(f"\nNo additional extention: FINAL LAYOUT with {len(seed_layout.reads)} reads, {seed_layout.end - seed_layout.begin} units:")
        #print(seed_layout)
        #print_layout(seed_layout, bits_dict)
        #layout_consensus(seed_layout) # TODO: temporary?
        #seed_layout.describe(asm)
        #print_layout(seed_layout, asm.bits)



# TODO: change signature to use alignmentstore obj
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
        newpool = Alignment(hers = [ ctx.hers[i] for (i, ai), k in self.reads ])
        regs = [ (i, ai) for i, a in newpool.arrs.items() for ai, l in enumerate(a) if len(l) > 4 ]
        units_in_layout = [ (i, mi) for i, ai in regs for mi in newpool.arrs[i][ai] if mi >= 0 ]
        #return newpool.get_variants(units_in_layout)
        return var(newpool.hers, units_in_layout)

    # TODO test
    def define_ends_variants(self, ctx, plusend = True):
        """ returns variants locally-defined at ends of the layout"""
        newpool = Alignment(hers = [ ctx.hers[i] for (i, ai), k in self.reads ])
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

        return var(newpool.hers, units_in_layout)

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
            np = Alignment(hers = [ccs] + ctx.hers[read[0]])
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

        #newpool = Alignment(hers = [self.get_consensus(ctx)] + [ ctx.hers[i] for (i, ai), k in self.reads ])
        newpool = Alignment(hers = [self.get_consensus(ctx)] + [ ctx.hers[i] for (i, ai), k in sorted(self.reads, key=lambda x: x[1]) ])
        regs = [ (i, ai) for i, a in newpool.arrs.items() for ai, l in enumerate(a) if len(l) > 4 ]
        layout_local_snvs = self.define_local_variants(ctx)
        print_snvs(layout_local_snvs)

        layout_local_bits_dict = newpool.get_bits_dict(layout_local_snvs, regs)
        l0 = layout_local_bits_dict[regs[0]].shape[0] # layout's length
        print("The structure of the layout:")
        print(newpool.hers[0].hors)

        def er_to_svs_dp(asm, i, ai, bd):
            # encoded read to self-vs-self dot plot in hor resolution, given snv mask.
            m = bd[(i, ai)]
            l = m.shape[0]
            def sc_xy(x, y):
                mx, my = m[x,:], m[y,:]
                mch = np.sum(np.multiply(mx, my) == 1)
                n11 = np.sum((mx + my) == 2)
                n00 = np.sum((mx + my) == -2)
                mis = np.sum(np.multiply(mx, my) == -1)
                assert mch == n11 + n00, "matrix calculation wrong"
                score = (n11 + ctx.c_00 * n00) / (0.01 + n11 + ctx.c_00 * n00 + ctx.c_ms * mis)
                return score
            # return [ (x, y, sc_xy(x, y)) for x in range(l) for y in range(x, l) ]
            return f"<svg width=\"{1.43*7*l+30}\" height=\"{7*l+30}\" transform=\"translate=(100, 0)\">" +\
                    "\n".join([ f'<rect x="{10+7*x:.2f}" y="{10+7*y:.2f}" width="8" height="8" fill-opacity="{sc_xy(x,y):.3f}" transform="rotate(-45, {10+7*x:.2f}, {10+7*y:.2f})"/>' for x in range(l) for y in range(x, l) ]) + "</svg>"


        i, ai = regs[0]
        dotplot = er_to_svs_dp(newpool, i, ai, layout_local_bits_dict)
        print("</pre>\n" + dotplot + "\n<pre>")

        bag_of_units = [ (ri, h) for ri, er in enumerate(ctx.hers) for h, _, t in er.hors if t == "~" ]
        glbd = newpool.get_bits_dict(ctx.get_variants(bag_of_units), regs)
        dotplot2 = er_to_svs_dp(newpool, i, ai, glbd)
        print("<br></pre>\n" + dotplot2 + "\n<pre>")

        return 

        for j, aj in regs[1:]:
            lj = layout_local_bits_dict[(j, aj)].shape[0]
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
    assert ctx, "needs Alignment as a context"

    bag_of_units = [ (ri, h) for ri, er in enumerate(ctx.hers) for h, _, t in er.hors if t == "~" ]
    n_units = len(bag_of_units)
    print(f"{n_units} units found in {len(ctx.hers)} reads.")

    #snv_sites = ctx.get_variants(bag_of_units)
    snv_sites = var(ctx.hers)

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
        print_layout(l, bits_dict)
        l.describe(ctx)
        pass

    import sys
    sys.exit()

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

def layout_to_json(alns_dict = None, ctx = None, layouts = None):
    """ converting the layout pickle into json; temporary """
    import json

    assert alns_dict, "needs Alignment."
    assert ctx, "needs Alignment as a context"

    def horUnits(hers, target):
        return [ (ri, h) for ri, er in enumerate(hers) for h, _, t in er.hors if t == target ]

    bag_of_units = horUnits(ctx.hers, "~")

    n_units = len(bag_of_units)
    print(f"{n_units} units found in {len(ctx.hers)} reads.")

    snv_sites = ctx.get_variants(bag_of_units)
    # NOTE: this amount to... snv_sites = ctx.get_variants(horUnits(ctx.hers, "~"))
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

    d = dict(reads = reads, layouts = [ l.reads.items() for l in layouts ])


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='generate layout from alignment results.')
    parser.add_argument('action', metavar='action', type=str, help='action to perform: layout, ...')
    # parser.add_argument('--hor-reads', dest='hors', help='pickled hor-encoded long reads')
    parser.add_argument('--alignments', dest='alns', help='pickled alignments')
    parser.add_argument('--layouts', dest='layouts', help='pickled layouts; it\'s your job to take care of consistency with alns...')
    args = parser.parse_args()

    if args.action == "layout":

        # assert args.hors, "need HOR-encoded reads"
        # hers = pickle.load(open(args.hors, "rb"))
        # assembly = Alignment(hers)

        assert args.alns, "need alignments"
        ast = pickle.load(open(args.alns, "rb"))
        context = Alignment(ast.reads, arrs = ast.arrays, variants = ast.variants)

        # TODO: this is working well
        #layouts = iterative_layout(alignments)
        #print(f"{len(layouts)} layouts found. done.")

        covers = read_cover(ast)

        n = 0
        layouts = []

        fluct = set([ r for c in covers for r in c if len(c) <= 9 ])

        print(f"{len(fluct)} reads in fluctured covers.")

        for cover in [ c for c in covers if len(c) > 9 ]:

            print(f"\n### {len(cover)} nodes in Cover {n} ###")
            print(cover)

            rds = [ ast.reads[i] for i in set([ ri for ri, rai in cover ]) ]
            snvs = var(rds)
            print(f"\n--- {len(snvs)} variants. ---")
            print_snvs(snvs)

            local_ast = context.get_all_vs_all_aln(cover | fluct, snvs)
            los = iterative_layout(local_ast)
            print(f"\n--- inside {len(los)} layouts for Cover {n} ---")
            for il, lo in enumerate(los):
                print("\n" + "\n".join([ f"{n}\t{il}\t{i}:{ai}\t{k}" for (i, ai), k in lo.reads ]))

            n += 1
            layouts += los

            #cover_local_ctx =  Alignment(ast.reads, arrs = ast.arrays, variants = snvs)

            # Alignment(reads = [ alignments.reads[i] for i in set() # NOTE: do not spawn new variatble space!!!!!

        with open("layouts-for-covers.pickle", "wb") as f:
            pickle.dump(layouts, f)

        # layout(alignments)

        # TODO: deprecated
        """
        if args.layouts:
            layouts = pickle.load(open(args.layouts, "rb"))
            layout(alns.alignments, ctx = assembly, layouts = layouts)
        else:
            layout(alns.alignments, ctx = assembly)
        """

        """
        elif args.action == "layout_to_json":
            assert args.hors, "need HOR-encoded reads"
            assert args.alns, "need alignments"
            assert args.layouts, "need layouts"
            hers = pickle.load(open(args.hors, "rb"))
            assembly = Alignment(hers)
            alns = pickle.load(open(args.alns, "rb"))
            layouts = pickle.load(open(args.layouts, "rb"))
            layout_to_json(alns.alignments, ctx = assembly, layouts = layouts)
        """

    else:
        assert False, "invalid action."
