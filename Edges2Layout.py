
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
from UnitAnalysis import *

def loadEdges(edges, score_t = 100, prom_t = 30, round_t = -1,
        no_lateral = False, only_dovetail = False):

    # colnames for edge list
    df = pd.read_csv(edges, sep="\t",
            names = ["I", "J", "Li", "Lj", "K", "Eov", "Fext", "Rext", "Score", "Prom", "nVars", "Round"])

    G = nx.MultiDiGraph()
    # TODO: parametrize initial filter
    for _k, row in [ (_k, r) for _k, r in df.iterrows()
            if float(r.Score) > score_t and float(r.Prom) > prom_t and int(r.Round) > round_t]:

        if no_lateral and int(row.K) == 0:
            continue

        if only_dovetail:
            if 0 <= int(row.K) and int(row.K) + int(row.Lj) <= int(row.Li):
                continue
            if int(row.K) <= 0 and int(row.K) + int(row.Lj) >= int(row.Li):
                continue

        if int(row.K) == 0:
            if int(row.I) < int(row.J):
                e = ( int(row.I), int(row.J), 0 )
            else:
                e = ( int(row.J), int(row.I), 0 )
        elif int(row.K) > 0:
            e = ( int(row.I), int(row.J), int(row.K) )
        else:
            e = ( int(row.J), int(row.I), -int(row.K) )

        # Takes only the edge with best score
        if e in G.edges and G.edges[e]["Score"] >= row.Score:
            continue

        if int(row.K) == 0:
            if int(row.I) < int(row.J):
                G.add_edge(int(row.I), int(row.J), key = 0)
                for k, v in row.items():
                    G.edges[int(row.I), int(row.J), 0][k] = v
            else:
                G.add_edge(int(row.J), int(row.I), key = 0)
                for k, v in row.items():
                    G.edges[int(row.J), int(row.I), 0][k] = v

        elif int(row.K) > 0:
            G.add_edge(int(row.I), int(row.J), key = int(row.K))
            for k, v in row.items():
                G.edges[int(row.I), int(row.J), int(row.K)][k] = v
        else:
            G.add_edge(int(row.J), int(row.I), key = -int(row.K))
            for k, v in row.items():
                G.edges[int(row.J), int(row.I), -int(row.K)][k] = v

    print(f"Loaded {len(G.nodes)} nodes and {len(G.edges)} edges.")
    # print(f"describe the graph here.")

    return G

def transitive_cycles(G, verbose = True):
    """ returns transitively reducible triples & contradictory cycles """

    # see farther edges earlier
    edges = sorted(list(G.edges), key=lambda x: (-x[2], x[0], x[1]))
    transitive_cycles = []
    contradictory_cycles = []
    n_bads = 0

    for e in edges:
        for _e, f, fk in list(G.out_edges(e[0], keys = True)):
            if (f, e[1], e[2] - fk) in G.edges:
                if verbose:
                    print(f"GOOD\t{e} = ({_e}, {f}, {fk}) + ({f}, {e[1]}, {e[2]-fk})",
                        flush = True)
                transitive_cycles += [(e[0], f, e[1])]

            for _f, g, gk in list(G.out_edges(f, keys = True)):
                if g == e[1] and fk + gk != e[2]:
                    if verbose:
                        print(f"BAD\t{e} != ({_e}, {f}, {fk}) + ({_f}, {g}, {gk})")
                    contradictory_cycles += [(e[0], f, e[1])]
                    n_bads += 1

    if verbose:
        print(f"{n_bads} intransitive; {len(transitive_cycles)} transitive cycles; " +\
              f"{len(set([ (a,c) for a,b,c in transitive_cycles ]))} reducible edges." )
    return (transitive_cycles, contradictory_cycles)

def dang_nodes(G, plus = True):
    """ pick nodes with degree 1, which are then reduced """
    if plus:
        return [ n for n in list(G.nodes) if G.in_degree(n) == 1 and G.out_degree(n) == 0 ]
    else:
        return [ n for n in list(G.nodes) if G.in_degree(n) == 0 and G.out_degree(n) == 1 ]

def chop_nodes(G, ns, plus = True):
    """ trim nodes found by dang_nodes """
    for n in ns:
        if plus:
            b = list(G.in_edges(n, keys = True, data = True))[0]
            if "dang" in G.nodes[b[0]]:
                G.nodes[b[0]]["dang"] += [(n, b[2])]
            else:
                G.nodes[b[0]]["dang"] = [(n, b[2])]
            if "dang" in G.nodes[n]:
                G.nodes[b[0]]["dang"] += [ (m, b[2] + k) for m, k in G.nodes[n]["dang"] ]
            G.remove_edge(b[0], n)

        else:
            b = list(G.out_edges(n, keys = True, data = True))[0]
            if "dang" in G.nodes[b[1]]:
                G.nodes[b[1]]["dang"] += [(n, -b[2])]
            else:
                G.nodes[b[1]]["dang"] = [(n, -b[2])]
            if "dang" in G.nodes[n]:
                G.nodes[b[1]]["dang"] += [ (m, -b[2] + k) for m, k in G.nodes[n]["dang"] ]
            G.remove_edge(n, b[1])
        G.remove_node(n)

def simple_nodes(G):
    return [ n for n in list(G.nodes) if G.in_degree(n) == 1 and G.out_degree(n) == 1 ]

def pop_nodes(G, ns):
    """ pop out simple nodes, unless simple nodes form 3-clique.
        3-clique may be formed when other simple nodes are popped out, thus it must be checked
        just before trying to pop them out. """

    for n in sorted(ns):
        b = list(G.in_edges(n, keys = True, data = True))[0]
        e = list(G.out_edges(n, keys = True, data = True))[0]

        if (b[0], e[1]) in G.edges:
            continue
        G.add_edge(b[0], e[1], key = b[2] + e[2])
        # NOTE: too high ???
        G.edges[b[0], e[1], b[2]+e[2]]["Score"] = b[3]["Score"] + e[3]["Score"]
        G.remove_edge(b[0], n)
        G.remove_edge(n, e[1])
        if "dang" in G.nodes[b[0]]:
            G.nodes[b[0]]["dang"] += [(n, b[2])]
        else:
            G.nodes[b[0]]["dang"] = [(n, b[2])]
        if "dang" in G.nodes[n]:
            G.nodes[b[0]]["dang"] += [ (m, b[2] + k) for m, k in G.nodes[n]["dang"] ]
        G.remove_node(n)

def remove_transitives(G, cycles):

    # remove goods (longer edges will be gone)
    for i, j, k in cycles:
        if (i, k) in G.edges(keys = False):
            G.remove_edge(i, k)

def remove_intransitives(G, cycles):
    # remove bad edges based on its score

    def pick_edge(i, j):
        f = list(filter(lambda x: x[1] == j,
                        list(G.edges(i, keys = True, data = True))))
        return f[0] if f else None

    for i, j, k in cycles:
        # print(f"bads; {i} {j} {k}")
        if ((i, j) not in G.edges(keys = False)) \
                or ((j, k) not in G.edges(keys = False)) \
                or ((i, k) not in G.edges(keys = False)):
            continue

        ij = pick_edge(i, j)[3]["Score"]
        jk = pick_edge(j, k)[3]["Score"]
        ik = pick_edge(i, k)[3]["Score"]

        if ij < jk:
            if ij < ik:
                G.remove_edge(i, j) # remove ij
            else:
                G.remove_edge(i, k) # remove ik 
        else: # jk <= ij
            if jk < ik:
                G.remove_edge(j, k) # remove jk
            else:
                G.remove_edge(i, k) # remove ik 


# pop_nodes(G, simple_nodes(G))
# remove_transitives(G):

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='Analyze edges list obtained.')
    parser.add_argument('action', metavar='action', type=str, help='layout, ...')

    parser.add_argument('--hor-reads', dest='hors', help='pickled hor-encoded long reads')
    parser.add_argument('--vars', dest='vars', help='pickled variant sites (disable auto detection)')
    parser.add_argument('--err-rate', dest='err_rate', help='error rate used in variant detecion')

    # NOTE it also accepts --err-rate in action layout?
    parser.add_argument('--edges', dest='edges', help='edge list file')
    parser.add_argument('--layouts', dest='layouts', help='precomputed layouts to be analysed')

    parser.add_argument('--params', dest='params', help='params for filtering edges')

    parser.add_argument('-o', dest='outfile', help='the file to be output (consensus reads pickle)')

    args = parser.parse_args()

    assert args.hors, "required --hor-reads; specify HOR-encoded reads"
    assert args.edges, "required --edges"

    hers = pickle.load(open(args.hors, "rb"))
    arrs = [ fillx(her) for her in hers ]
    units = [ (her, h) for her in hers for h, s, t in fillx(her) if t == "~" ]

    # Time when invoked.
    import datetime
    datetimestr = datetime.datetime.now().strftime("%Y.%m%d.%H%M")

    if args.action == "transitivity": # check transitivity


        if args.params:
            ap = args.params.split(",")
            score_t, prom_t, round_t = int(ap[0]), int(ap[1]), int(ap[2])
        else:
            score_t, prom_t, round_t = 100, 30, -1

        print(f"loading edges with params = s>{score_t}, p>{prom_t}, r>{round_t}")
        G = loadEdges(args.edges, score_t, prom_t, round_t)

        print(f"name           \t#node\t#edge\t#simp\t#+dang\t#-dang\t#t_cyc\t#i_cyc")

        def describe(G, name = "Graph"):
            tcycles, icycles = transitive_cycles(G, verbose = False)
            simples = simple_nodes(G)
            plus_dang, minus_dang = dang_nodes(G, True), dang_nodes(G, False)
            print(f"{name:15}\t{len(G.nodes)}\t{len(G.edges)}\t{len(simples)}\t" +\
                  f"{len(plus_dang)}\t{len(minus_dang)}\t{len(tcycles)}\t{len(icycles)}")

        def clusters(G):
            """ collapsed nodes stored as "dang" attribute represents one consistent layout. """
            clusters = [ [(n, 0)] + G.nodes[n]["dang"] for n in G.nodes if "dang" in G.nodes[n] ]
            print(",".join([ f"{len(c)}" for c in sorted(clusters, key = lambda x: -len(x)) ][:20]))

        # NOTE: remove icycles nodes from simples  # TODO
        describe(G, "initial")
        clusters(G)

        # remove isolated nodes # NOTE: unnecessarry
        #for n in list(G.nodes):
        #    if G.out_degree(n) == 0 and G.in_degree(n) == 0:
        #        G.remove_node(n)

        describe(G, "initial")
        clusters(G)

        iteration = 0
        while iteration < 9:

            iteration += 1
            tcycles, icycles = transitive_cycles(G, verbose = False)
            remove_transitives(G, tcycles)
            describe(G, f"no-tr-edge.{iteration}")
            clusters(G)
            
            last_nedges = 100000000
            while last_nedges > len(G.edges):
                last_nedges = len(G.edges)
                chop_nodes(G, dang_nodes(G, True), True)
                chop_nodes(G, dang_nodes(G, False), False)

            describe(G, f"no-dang.{iteration}")
            clusters(G)

            pop_nodes(G, simple_nodes(G))
            describe(G, f"pop-simple.{iteration}")
            clusters(G)

            tcycles, icycles = transitive_cycles(G, verbose = False)
            remove_intransitives(G, icycles)
            describe(G, f"!no-int-edge.{iteration}")

    elif args.action == "layout": # perform layout

        # NOTE: do i need to change err_rate (if it affect anything)?
        v_major = var(hers, hor_type = "~", err_rate = args.err_rate if args.err_rate else 0.05,
            fq_upper_bound = 1.1, comprehensive = False)
        v_all = var(hers, hor_type = "~", err_rate = 0.05, fq_upper_bound = 1.1, comprehensive = True)

        # positions of non-canonical units in read i (index and types)
        nonstd = { i:  [ (ii, t) for ii, (h, s, t) in enumerate(a) if t != "~" ] for i, a in enumerate(arrs) }
        G = loadEdges(args.edges, 90, 30, 0)

        sys.exit()

        def merge_k(l1, l2, k):
            """ merge 2 layouts by offset k; the very first position is set to k=0 """
            s = sorted(l1 + [ (n, k+l) for n, l in l2 ], key = lambda x: x[1])
            return [ (n, l - s[0][1]) for n, l in s ]

        def bridge(l1, l2, ins = None):
            """ estimate offset from l1 to l2 ...
                here, ins is the edges within/between l1/l2 
                returns Counter object.
            """
            ins = ins if ins else G.edges(keys = True)
            ls1 = { li: lk for li, lk in l1 }
            ls2 = { li: lk for li, lk in l2 }
            _edges = [ (gi, gj, k) for gi, gj, k in ins if gi in ls1 or gj in ls2 ]

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

        # NOTE: this is too specific to X
        def consensus(layout, verbose = False):
            """ take a consensus of the layout """
            m, M = 10000, -10000
            types = { i : Counter() for i in range(-10000, 10000) }

            for li, lk in layout:
                #print(f"arr {li}@{lk} = {arrs[li]}")
                for ii, (h, s, t) in enumerate(arrs[li]):
                    M, m = max(M, lk + ii), min(m, lk + ii)
                    types[lk+ii].update([t])

            units = { i: { k : Counter() for k in range(12) } for i in range(m, M+1) } # NOTE X
            depth = { i: 0 for i in range(m, M+1) }

            for li, lk in layout:
                if all([ True or (t == types[lk+ii].most_common()[0][0]) or t == "*"
                         for ii, (h, s, t) in enumerate(arrs[li]) ]):
                    for ii, (h, s, t) in enumerate(arrs[li]):
                        if t == "~":
                            depth[lk+ii] += 1
                            for mk in range(0, 12):# NOTE X
                                units[lk+ii][mk].update(hers[li].mons[h + mk].monomer.snvs)

            cons_hors = [ (i * 12, 12, types[i+m].most_common()[0][0]) for i in range(0, M-m+1) ] # NOTE X
            cons_mons = [ AssignedMonomer(begin=lk*2057 + k*171, end=lk*2057 + (k+1)*171, ori="+",
                monomer = Monomer(
                    name = f"Consensus-{lk}-{k}",
                    snvs = [ s for s, sc in units[lk][k].most_common() if sc / depth[lk] > 0.4 ]))
                for lk in range(m, M+1) for k in range(12) ]

            hor_read = HOR_Read(name = "Consensus",
                            mons = cons_mons, hors = cons_hors,
                            length = (M-m)*2057, ori="+")

            if verbose:

                ## total (default) units 
                coverage = sum([ depth[lk] for lk in range(m, M+1) ])

                print("recall, specificity?, fpr, fnr")
                print("\ni\ttps(%)\ttns(%)\tfps(%)\tfns(%)\tnvars\tdepth\ttypes observed.")

                vfs = { lk: 
                    [ sc for k in range(12) for s, sc in units[lk][k].most_common() ] # NOTE X
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

                print("readname\tbegin\tend\tidx\tsize\telem\tgap\tvars")
                print_HOR_read(hor_read)

            return hor_read

        def show_layout(layout):
            """ visualize consensus read """

            cons_read = consensus(layout)
            cons_arr = fillx(cons_read)
            # NOTE: do I need this?
            snvs = var([ hers[li] for li, lk in layout ])
            print_snvs(snvs)

            snvs_read = var([cons_read], err_rate = 0.01, comprehensive = True) # NOTE: this must be true!
            print(f"Total {len(snvs_read)} SNVs on this layout.")
            print_snvs(snvs_read)

            # t2c = {"*": "*", "~": "", "D1":"Y", "D12":"X", "22U":"U", "D39":"V", "D28":"W"}

            labels = [ t2c[t] for h, s, t in cons_read.hors ]

            # using reads SNVs
            for vs, name in [(v_all, "global"), (v_major, "gl-freq"), (snvs, "local"), (snvs_read, "consensus")]:
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
                for n, (h, s, t) in enumerate(arrs[li]):
                    read.add(
                        dwg.rect(insert=((n + lk - lkmin)*5,i*5), size=(4,4),
                        fill=t2col[t], stroke='black', stroke_width=0.5))
                    # ax1.text((lk + n) * 1, i * 1, t2c[t], fontsize=9) # not svgwrite
                    # NOTE: add label...

            dwg.save()

            return cons_read

        # if layouts is given
        if args.layouts:

            # say something about global mismatch distribution.
            layouts = pickle.load(open(args.layouts, "rb"))
            layouts = sorted(layouts, key = lambda x: -len(x))

            cons_reads = []
            for i in range(len(layouts)): 
                if len(layouts[i]) > 4:
                    print(f"consensus for {i}")
                    # cons_reads += [ show_layout(layouts[i]) ]
                    cons_reads += [ consensus(layouts[i]) ]

            # TODO
            if args.outfile:
                pickle.dump(cons_reads, open(args.outfile, "wb"))

            sys.exit()

        # Initialize
        layouts = [ [(n, 0)] for n in list(G.nodes) ]

        ins = [ (gi, gj, k)
            for l in layouts
            for gi, gj, k in G.edges(keys = True)
            if len(l) > 1 and {gi, gj} <= { li[0] for li in l } ]

        outs = [ (gi, gj, k) for gi, gj, k in G.edges(keys = True) if (gi, gj, k) not in set(ins) ]
        # eov sort, filter
        outs = sorted(outs, key = lambda x: -G.edges[x]["Prom"])

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
        
        pickle.dump(layouts, open(f"layouts-{datetimestr}.pkl", "wb"))
        print("done")

    else:
        assert False, "invalid action."
