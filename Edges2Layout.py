
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

def describe(G, name = "Graph"):
    tcycles, icycles = transitive_cycles(G, verbose = False)
    simples = simple_nodes(G)
    plus_dang, minus_dang = dang_nodes(G, True), dang_nodes(G, False)
    clens = sorted([ len(comp) for comp in nx.weakly_connected_components(G) ], reverse = True)
    print(f"{name:15}\t{len(G.nodes)}\t{clens[0]} ({100* clens[0] / sum(clens):.1f})\t" +\
          f"{len(G.edges)}\t{len(simples)}\t" +\
          f"{len(plus_dang)}\t{len(minus_dang)}\t{len(tcycles)}\t{len(icycles)}")

def largest_component(G):
    ln = sorted([ (len(comp), comp) for comp in nx.weakly_connected_components(G) ],
            reverse = True)[0][1]
    return G.subgraph(ln).copy()

def danglings(G, verbose = True):

    n_dangs = [ (n, G.nodes[n]["dang"]) for n in G.nodes if "dang" in G.nodes[n] ] #node_with_dangs
    n_dangs = sorted(n_dangs, key = lambda x: -len(x[1]))

    def gl_dang(d):
        """ genomic distance in units """
        return max([ dd[1] for dd in d ])

    sizes = [ len(n[1]) + 1 for n in n_dangs ]
    gls = [ gl_dang(n[1]) + 1 for n in n_dangs ]

    if verbose:
        print(f"{len(n_dangs)} nodes with {sum(sizes)} dangs, spanning {sum(gls)} units in total.")
        print("[(#nodes, gen.len)][:10] = ")
        print(f"{[ z for z in zip(sizes[:10], gls[:10]) ]}...")

    return n_dangs # list of tuple as in (parent node, list of dangling nodes = (nodeid, rel.pos))

def loadEdges(edges, score_t = 100, prom_t = 30, round_t = -1,
        no_lateral = False, only_dovetail = False):

    # colnames for edge list
    df = pd.read_csv(edges, sep="\t",
            names = ["I", "J", "Li", "Lj", "K", "Eov",
                     "Fext", "Rext", "Score", "Prom", "nVars", "Round"])

    valid_edges = [ (_k, r) for _k, r in df.iterrows()
        if float(r.Score) > score_t and \
           float(r.Prom) > prom_t and \
           int(r.Round) > round_t ]

    G = nx.MultiDiGraph()

    for _k, row in valid_edges:

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

# NOTE: deprecate
def dang_nodes(G, plus = True):
    """ pick nodes with degree 1, which are then reduced """
    if plus:
        return [ n for n in list(G.nodes) if G.in_degree(n) == 1 and G.out_degree(n) == 0 ]
    else:
        return [ n for n in list(G.nodes) if G.in_degree(n) == 0 and G.out_degree(n) == 1 ]

# NOTE: deprecate
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

# NOTE: deprecate
def simple_nodes(G):
    return [ n for n in list(G.nodes) if G.in_degree(n) == 1 and G.out_degree(n) == 1 ]

# NOTE: deprecate
def pop_nodes(G, ns):
    """ pop out simple nodes, unless simple nodes form 3-clique.
        3-clique may be formed when other simple nodes are popped out, thus it must be checked
        just before trying to pop them out. """

    for n in sorted(ns):
        b = list(G.in_edges(n, keys = True, data = True))[0]
        e = list(G.out_edges(n, keys = True, data = True))[0]

        if (b[0], e[1]) in G.edges:
            continue # triple

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

        # NOTE: equivalent to merge_nodes(G, b[0], n, b[2]) ?

def merge_nodes(G, n, m, k):
    """ m is merged into n """

    if "dang" in G.nodes[n]:
        G.nodes[n]["dang"] += [(m, k)]
    else:
        G.nodes[n]["dang"] = [(m, k)]
    # TODO: check why I can't reach here ?!
    if "dang" in G.nodes[m]:
        s = len(G.nodes[m]['dang']) + len(G.nodes[n]['dang']) + 1
        print(f"{s}", end = "")
        G.nodes[n]["dang"] += [ (l, k + lk) for l, lk in G.nodes[m]["dang"] ]
    G.remove_node(m)

def prioritise_pair(G):
    """ require: G has no 3-clique. or should I check that here? (or blacklist might be passed?) """

    def enumerate_triples(G): # e is a form of (s, t, k) or (s, t)
        triples = []
        for e in G.edges():
            for _e, f, fk in list(G.out_edges(e[0], keys = True)):
                for _f, g, gk in list(G.out_edges(f, keys = True)):
                    if g == e[1]:
                        triples += [(e[0], e[1]), (_e, f), (_f, g)]

        # print(f"len(triples) = {len(triples)}")
        return triples

    triples = enumerate_triples(G)

    el = []
    for n, m, k in G.edges(keys = True):
        if (n, m) in triples or (m, n) in triples:
            continue
        el += [ (G.edges[n, m, k]["Score"], n, m, k) ]

    return max(el) if el else None

def remove_nontriple(G, n, m, k, arrs = None):
    """ if no triple (transitively consistent or not) is found,
    then any pair of neighboring nodes can be reduced into single nodes,
    without affecting regularity of the graph. such reduction would mean we trust
    that edge (as we cannot inspect it anymore), so I look for best edges as candidates
    for reduction pair (n, m); see above. 
    """

    assert (n, m, k) in G.edges(keys = True), f"no edge {(n, m, k)} found."

    """ if necessarry, it checks proper alignment of special units using arrs data given. """
    # TODO: check mergeability here, before any update on graph

    if arrs:
        drydang = [(n, 0), (m, k)]
        if "dang" in G.nodes[n]:
            drydang += G.nodes[n]["dang"]
        if "dang" in G.nodes[m]:
            drydang += [ (di, dk + k) for di, dk in G.nodes[m]["dang"] ]

        types = { i : Counter() for i in range(-10000, 10000) }

        for li, lk in drydang:
            for ii, (h, s, t) in enumerate(arrs[li]):
                if t == "*":
                    continue # I allow wildcard units
                if t == "D28":
                    types[lk+ii].update(["D39"]) # D28 is assumed to be the same as D39
                else:
                    types[lk+ii].update([t])

        if any([ len(types[i]) > 1 for i in types.keys() ]):
            print("got conflict l.289")
            G.remove_edge(n, m, k) # I dont wanna see you again.
            return -1 # without doing anything
        # print(drydang)

    G.remove_edge(n, m, k)
    for mm, mj, mk in G.out_edges(m, keys = True):
        G.add_edge(n, mj, key = k + mk) # TODO add data
        G.edges[n, mj, k + mk]["Score"] = G.edges[mm, mj, mk]["Score"]
    for mi, mm, mk in G.in_edges(m, keys = True):
        if mk - k >= 0:
            G.add_edge(mi, n, key = mk - k) # TODO add data
            G.edges[mi, n, mk - k]["Score"] = G.edges[mi, mm, mk]["Score"]
        else:
            G.add_edge(n, mi, key = k - mk) # TODO add data
            G.edges[n, mi, k - mk]["Score"] = G.edges[mi, mm, mk]["Score"]

    # TODO: record as dang of n
    merge_nodes(G, n, m, k)

def remove_nontriple_all(G, arrs = None):
    # TODO make clear
    n, e = 0, prioritise_pair(G)
    while e:
        remove_nontriple(G, e[1], e[2], e[3], arrs = arrs)
        n, e = n + 1, prioritise_pair(G)
        print(".", end = "", flush = True)

    print(f"\nremoved {n} nodes")
    describe(G)
    return n


def remove_transitives(G, cycles = None):

    if not cycles:
        cycles, i = transitive_cycles(G, verbose = False)
    # remove goods (longer edges will be gone)
    n = 0
    for i, j, k in cycles:
        if (i, k) in G.edges(keys = False):
            G.remove_edge(i, k)
            n += 1
            # NOTE: to record data to existing edge, I need to consider the order of elimination
    describe(G)
    return n

def iter_tr(G, arrs = None):
    n = remove_transitives(G)
    while n > 0:
        n = remove_nontriple_all(G, arrs = arrs)
        n += remove_transitives(G)

def remove_intransitives(G, cycles = None):
    # remove bad edges based on its score

    if not cycles:
        i, cycles = transitive_cycles(G, verbose = False)

    def pick_edge(i, j):
        f = list(filter(lambda x: x[1] == j,
                        list(G.edges(i, keys = True, data = True))))
        return f[0] if f else None

    for i, j, k in cycles:
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
    describe(G)

def limit_outedges(G):
    """ is this a final resort?? """
    pass

def layout190512(G, prefix = "layouts", arrs = None):
    """ returns layout :: [(nodeid, pos)],
    which is then pickled and passed to visualization and stuff """

    describe(G) # description of initial graph loaded

    # let's focus on the largest components
    G = largest_component(G)
    describe(G)

    ## then I believe it's not so harmful to remove in-transitive node firstly;
    ## but let's make sure firstly this won't be harmful really.
    # remove_transitives(G)
    # remove_intransitives(G)

    # then, iteratively reduce, edges & nodes ...
    iter_tr(G, arrs = arrs)

    nd = sorted(danglings(G, verbose = True),
            key = lambda x: -x[1][-1][1]) # sort by the supposed genomic length
    layouts = [ [(d[0], 0)] + d[1] for d in nd ]
    pickle.dump(layouts, open(prefix + ".pkl", "wb"))

    # then, check the resultant is OK as a set of layout.

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
    parser.add_argument('--prefix', dest='prefix', help='prefix of consensus read name')
    parser.add_argument('--noc', dest='noc', action='store_const', const=True,
            help='not allow conflict of special HOR vars in layouting')
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

    if args.action == "layout190512": # make layout !

        if args.params:
            ap = args.params.split(",")
            score_t, prom_t, round_t = int(ap[0]), int(ap[1]), int(ap[2])
        else:
            score_t, prom_t, round_t = 100, 30, -1

        print(f"loading edges with params = s>{score_t}, p>{prom_t}, r>{round_t}")
        G = loadEdges(args.edges, score_t, prom_t, round_t)

        # TODO: change?
        print(f"name\t#node\t#in-largest\t#edge\t#simp\t#+dang\t#-dang\t#t_cyc\t#i_cyc")

        # NOTE: maybe I want to change its logic a bit
        #layout190512(G, prefix = datetimestr + f".{score_t}.{prom_t}.{round_t}")
        if args.noc:
            layout190512(G, prefix = datetimestr + f".{score_t}.{prom_t}.{round_t}-noc", arrs = arrs)
        else:
            layout190512(G, prefix = datetimestr + f".{score_t}.{prom_t}.{round_t}", arrs = None)

        print("All done.")

    else:
        assert False, "invalid action."
