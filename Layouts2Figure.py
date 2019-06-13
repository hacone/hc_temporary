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

t2c = {"*": "*", "~": "", "D1":"Y", "D12":"X", "22U":"U", "D39":"V", "D28":"W"}

# TODO: I just need datatype definition. That can be separated from other codes.
from EncodedRead import *
from HOR_segregation import *
from Edges2Layout import *
from UnitAnalysis import *

gid = 0

def draw_read(read, arr, snvs, name = "read"):

    t2col2 = {"*": "black", "~": "grey", "D1":"yellow", "D12":"green",
             "22U":"brown", "D39":"orange", "D28":"pink"}

    bits = ba(read, arr, snvs)
    r = bits[:, 0:100] # matrix representation of each read
    li = r.shape[0]
    ra = arr # arr data for HOR variants
    nvars = r.shape[1]
    dwg = svgwrite.Drawing(filename = f"read-{read.name}-{nvars}-{name}.svg")

    def draw_unit(r, i, x, y, col = "grey"):
        global gid
        bit2col = { -1: "white", 0: "#cccccc", 1: "black" }
        if i > -1: 
            unit_a = dwg.add(dwg.g(id=f"a-{i}-{gid}"))
            gid += 1
            unit_a.add(dwg.rect(insert = (x, y), size = (5 + nvars, 10),
                stroke_width = 1, fill = col))
            for n, e in enumerate(r[i,:]):
                unit_a.add(dwg.line(
                    start = (x + 2 + n, y + 1),
                    end   = (x + 2 + n, y + 9),
                    stroke = bit2col[e]))

    xcoord, ycoord, xskip = 0, 20, r.shape[1] + 10

    for i in range(li):
        draw_unit(r, i, xcoord * xskip, ycoord + 20, col=t2col2[ra[i][2]])
        dwg.add(dwg.text(f"{i+1}", insert=(xcoord * xskip, ycoord + 40)))
        xcoord += 1

    dwg.save()

    import scipy.cluster.hierarchy as shc
    bits_def = np.array([ bits[i,:] for i in range(li) if ra[i][2] == "~" ])

    #plt.figure(figsize=(10, 7))  
    fs = max(int( 10 * li / 20 ), 10)
    plt.figure(figsize=(int(fs*0.6), fs))
    plt.title(f"{read.name}-{nvars}-{name}-dendro")
    linkage = shc.linkage(bits_def, metric = "hamming",
            method='single', optimal_ordering = True)
    linkage = linkage + np.array([0, 0, 0.001, 0]*linkage.shape[0]).reshape(linkage.shape)
    dend = shc.dendrogram(linkage, orientation = "left",
            labels = [ f"{i+1}" for i in range(li) if ra[i][2] == "~"]) 
    plt.savefig(f"read-{read.name}-{nvars}-{name}-dendro." +\
            ("svg" if args.filetype == "svg" else "png"))
    plt.close()

# NOTE: this is too specific to X; TODO: expose me!
# TODO: correlate with given hors
def consensus(layout, hers, arrs, verbose = False, name = "Consensus"):
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

    hor_read = HOR_Read(name = f"{name}",
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

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='Visualize & make stats on layouts / contigs')
    parser.add_argument('action', metavar='action', type=str,
            help='action to perform: align, show, ...')

    # for action: show
    parser.add_argument('--hor-reads', dest='hors', help='pickled hor-encoded long reads')
    parser.add_argument('--layouts', dest='layouts', help='pickle of layouts to be visualized')
    parser.add_argument('--onlys', dest='onlys', help='ids of layouts to be visualized')
    parser.add_argument('--v-scale', dest='v_scale', help='scale of dot plot as in: 99.5,100')
    parser.add_argument('--err-rate', dest='err_rate', help='err. rate assumed in calc of gl. var. freq.')

    # for action: align
    parser.add_argument('--a-layouts', dest='alayouts', help='A layouts to be aligned')
    parser.add_argument('--b-layouts', dest='blayouts', help='B layouts to be aligned')
    parser.add_argument('--a-hor-reads', dest='ahors', help='reads for A layouts')
    parser.add_argument('--b-hor-reads', dest='bhors', help='reads for B layouts')
    parser.add_argument('--filetype', dest='filetype', help='png (default) or svg')

    # parser.add_argument('--consensi', dest='consensi', help='comma-separated list of pickled hor-encoded long reads (consensi)')
    #parser.add_argument('--vars', dest='vars', help='pickled variant sites (disable auto detection)')
    args = parser.parse_args()

    # NOTE: this assumes arrs in context
    def show_layout(layout, hers, arrs, cons_read = None,
            snvs_list = None, v_scale = None, force_denom = False):
        """ visualize consensus read """

        if not cons_read:
            cons_read = consensus(layout, hers, arrs)
        cons_arr = fillx(cons_read)

        # t2c = {"*": "*", "~": "", "D1":"Y", "D12":"X", "22U":"U", "D39":"V", "D28":"W"}
        labels = [ "." if (t == "~" and (i+1)%10 == 0) else t2c[t]
                for i, (h, s, t) in enumerate(cons_read.hors) ]

        if not snvs_list:

            # NOTE: do i need to change err_rate (if it affect anything)?
            v_major = var(hers, hor_type = "~",
                err_rate = args.err_rate if args.err_rate else 0.05,
                fq_upper_bound = 1.1, comprehensive = False)
            v_all = var(hers, hor_type = "~",
                err_rate = 0.05, fq_upper_bound = 1.1, comprehensive = True)
            snvs = var([ hers[li] for li, lk in layout ])
            snvs_read = var([cons_read], err_rate = 0.01, comprehensive = True) # NOTE: this must be true!
            print(f"Total {len(snvs_read)} SNVs on this layout.")

            snvs_list = [
                (v_all, "all"),
                (v_major, "major"),
                (snvs, "local"),
                (snvs_read, "consensus"),
                (None, "naive")]

        for snvs, name in snvs_list:
            if not snvs:
                continue
            draw_read(cons_read, cons_arr, snvs, name = f"{layout[0][0]}-" + name)

        # NOTE: calculate depth etc.
        lkmin = min([ lk for li, lk in layout ])
        depth = Counter()
        for i, (li, lk) in enumerate(sorted(layout, key = lambda x: x[1])):
            for n, (h, s, t) in enumerate(arrs[li]):
                depth[n+lk] += 1 if t == "~" else 0

        # using reads SNVs
        for vs, name in snvs_list:

            if not force_denom:
                dots = acomp(cons_read, cons_arr, cons_read, cons_arr, snvs = vs)
            else:
                dots = acomp(cons_read, cons_arr, cons_read, cons_arr, snvs = vs, force_denom = 2057)

            # Dotplot of layout
            fs = max(int( 20 * len(dots) / 100 ), 20)
            fig = plt.figure(figsize=(fs, int(fs*0.8)))
            sns.set(font_scale=2)
            ax1 = fig.add_subplot(1, 1, 1)
            if not v_scale:
                v_scale = (np.nanmin(dots), np.nanmax(dots))
            g1 = sns.heatmap(dots,
                    vmin=v_scale[0], vmax=v_scale[1],
                    cmap="coolwarm", ax = ax1,
                    xticklabels = labels, yticklabels = labels)
            ax1.set_title(f"Layout-{layout[0][0]}; {len(layout)} rds;" +\
                    f"{len(cons_arr)} units; with {len(vs) if vs else 'all'} of {name} vars")
            plt.savefig(f"Layout-{layout[0][0]}-{name}." +\
                    ("svg" if args.filetype == "svg" else "png"))
            plt.close()

            # TODO: dendrogram by similarity of units
            #import AgglomerativeClustering from sklearn.cluster
            #AgglomerativeClustering(affinity = "precomputed", linkage = "single")
            print(f"dots.shape = {dots.shape}")

            # level of local variations
            cps = [ dots[i,i+1] for i in range(len(dots)-1) if depth[i] > 2 and depth[i+1] > 2 ]
            nanmean = np.nanmean(cps)
            print(f"AVGSIM\t{name}\t{len(cps)}\t{nanmean*100:.3f}")


        ## print out its structure / CENP-B presense
        dwg = svgwrite.Drawing(filename=f"Layout-{layout[0][0]}-str.svg")

        for i, (li, lk) in enumerate(sorted(layout, key = lambda x: x[1])):
            read = dwg.add(dwg.g(id=f"rd{li}", stroke='green'))
            dwg.add(dwg.text(f"{li}", insert=(5 + (lk - lkmin)*25, 40 + i*25) ))
            for n, (h, s, t) in enumerate(arrs[li]):
                read.add(
                    dwg.rect(insert=(40 + (n + lk - lkmin)*25, 30 + i*25), size=(20, 20),
                    fill=t2col[t], stroke='black', stroke_width=2))

        dwg.add(dwg.text(f"#V", insert=(5, len(layout)*25 + 165) ))
        dwg.add(dwg.text(f"Dep.", insert=(5, len(layout)*25 + 175) ))
        dwg.add(dwg.text(f"#CBS", insert=(5, len(layout)*25 + 185) ))
        dwg.add(dwg.text(f"#CBM", insert=(5, len(layout)*25 + 195) ))

        for n, (h, s, t) in enumerate(cons_arr):

            if t != "~":
                continue
            intact_sites, intact_motifs = cenpb_vars(cons_read, h, summary = True)

            height = 3 * (int(intact_sites) - 40) # NOTE: baseline is 40
            dwg.add(dwg.rect(
                insert=((n - lkmin)*25 + 40, 25*len(layout) + 110 - height),
                size=(20, height))) # NOTE: baseline is 40

            height = int(intact_motifs * 6)
            dwg.add(dwg.rect(
                insert=((n - lkmin)*25 + 40, 25*len(layout) + 155 - height),
                size=(20, height)))

            nvars = sum([ len(cons_read.mons[h+k].monomer.snvs)  for k in range(12) ])
            dwg.add(dwg.text(
                f"{nvars}", insert=((n - lkmin)*25 + 40, 25*len(layout) + 165) ))

            dwg.add(dwg.text(
                f"{depth[n]}", insert=((n - lkmin)*25 + 40, 25*len(layout) + 175) ))

            dwg.add(dwg.text(
                f"{int(intact_sites)}", insert=((n - lkmin)*25 + 40, 25*len(layout) + 185) ))

            dwg.add(dwg.text(
                f"{int(intact_motifs)}", insert=((n - lkmin)*25 + 40, 25*len(layout) + 195) ))

            print(f"VVD\t{depth[n]}\t{nvars}")

        # line at max of #intact cenp b motifs
        dwg.add(dwg.line(
            start=((0 - lkmin)*25 + 30, 25*len(layout) + 110 - 69),
            end=((len(cons_arr) - lkmin)*25 + 30, 25*len(layout) + 110 - 69),
            stroke="black", stroke_width=2)),

        # line at max of #intact cenp b sites
        dwg.add(dwg.line(
            start=((0 - lkmin)*25 + 30, 25*len(layout) + 155 - 42),
            end=((len(cons_arr) - lkmin)*25 + 30, 25*len(layout) + 155 - 42),
            stroke="black", stroke_width=2)),

        # dummy
        dwg.add(dwg.rect(
            insert=(25*(len(cons_arr)-lkmin) + 30, 25*len(layout) + 200), size=(30,30),
            stroke = 'black', stroke_width = 5))

        dwg.save()

        return cons_read

    def dotplot_lvl(aread, bread, snvs_list, filetype = None, v_scale = None):
        # output png when filetype is None

        aarr = fillx(aread)
        barr = fillx(bread)

        a_labels = [ t2c[t] for h, s, t in aread.hors ]
        b_labels = [ t2c[t] for h, s, t in bread.hors ]

        # using reads SNVs
        for vs, name in snvs_list:
            dots = squarify(acomp(aread, aarr, bread, barr, snvs = vs))

            if not v_scale:
                v_scale = (np.nanmin(dots), np.nanmax(dots))

            fs = max(int( 20 * len(dots) / 100 ), 20)
            fig = plt.figure(figsize=(fs, int(fs*0.8)))
            sns.set(font_scale=2)
            ax1 = fig.add_subplot(1, 1, 1)
            g1 = sns.heatmap(dots,
                    vmin=v_scale[0], vmax=v_scale[1],
                    cmap="coolwarm", ax = ax1,
                    xticklabels = b_labels, yticklabels = a_labels)
            ax1.set_title(f"{aread.name} - {bread.name}; {len(aarr)} vs. {len(barr)} units; with {len(vs)} of {name} vars")

            plt.savefig(f"{aread.name}-{bread.name}-{name}." +\
                    ("svg" if filetype == "svg" else "png"))
            plt.close()

    if args.action == "align":
        # TODO: need output restriction here too.

        assert args.ahors, "specify HOR-encoded reads"
        assert args.bhors, "specify HOR-encoded reads"
        ahers = pickle.load(open(args.ahors, "rb"))
        aarrs = [ fillx(her) for her in ahers ]
        bhers = pickle.load(open(args.bhors, "rb"))
        barrs = [ fillx(her) for her in bhers ]

        # units = [ (her, h) for her in hers for h, s, t in fillx(her) if t == "~" ]
        assert args.alayouts, "no A layouts specified, abort."
        assert args.blayouts, "no B layouts specified, abort."
        alayouts = pickle.load(open(args.alayouts, "rb"))
        alayouts = sorted(alayouts, key = lambda x: -len(x))
        blayouts = pickle.load(open(args.blayouts, "rb"))
        blayouts = sorted(blayouts, key = lambda x: -len(x))

        def glen(l):
            # genomic length (in units) of layouts
            ul = [ p for i, p in l ]
            return max(ul) - min(ul) + 1

        print(f"aligns {len([l for l in alayouts if glen(l) > 9])} of {len(alayouts)} A layouts" +\
              f"against {len([l for l in blayouts if glen(l) > 9])} of {len(blayouts)} B layouts.")

        # NOTE: can these be precalculated?
        v_major = var(ahers + bhers, hor_type = "~",
            err_rate = args.err_rate if args.err_rate else 0.05,
            fq_upper_bound = 1.1, comprehensive = False)
        v_all = var(ahers + bhers, hor_type = "~",
            err_rate = 0.05, fq_upper_bound = 1.1, comprehensive = True)

        if args.v_scale:
            v_scale = [ float(s) for s in args.v_scale.split(",") ]
        else:
            v_scale = None

        # for ca in cons_a:
        for i in range(len(alayouts)):
            cons_a = consensus(alayouts[i], ahers, aarrs,
                    name = f"A-CS{i}:{alayouts[i][0][0]}")
            if glen(alayouts[i]) < 10:
                continue

            #for cb in cons_b:
            for j in range(len(blayouts)):
                cons_b = consensus(blayouts[j], bhers, barrs,
                        name = f"B-CS{j}:{blayouts[j][0][0]}")
                if glen(blayouts[j]) < 10:
                    continue

                # NOTE: assuming ahers and bhers are disjoint
                snvs = var([ ahers[li] for li, lk in alayouts[i] ] +\
                           [ bhers[li] for li, lk in blayouts[j] ])
                snvs_read = var([cons_a, cons_b], err_rate = 0.01, comprehensive = True)
                snvs_list = [ (v_all, "all"), (v_major, "major"),
                    (snvs, "local"), (snvs_read, "consensus"), (None, "naive") ]
                dotplot_lvl(cons_a, cons_b, snvs_list, filetype = args.filetype, v_scale = v_scale)

    if args.action == "show":

        assert args.hors, "specify HOR-encoded reads"
        hers = pickle.load(open(args.hors, "rb"))
        arrs = [ fillx(her) for her in hers ]
        # units = [ (her, h) for her in hers for h, s, t in fillx(her) if t == "~" ]
        # positions of non-canonical units in read i (index and types)
        nonstd = { i:  [ (ii, t) for ii, (h, s, t) in enumerate(a) if t != "~" ]
            for i, a in enumerate(arrs) } # NOTE: unused??

        assert args.layouts, "no layouts specified, abort."
        layouts = pickle.load(open(args.layouts, "rb"))
        layouts = sorted(layouts, key = lambda x: -len(x))

        # NOTE: precalculate snvs_list
        v_all = var(hers, hor_type = "~",
            err_rate = 0.05, fq_upper_bound = 1.1, comprehensive = True)
        v_major = var(hers, hor_type = "~",
            err_rate = args.err_rate if args.err_rate else 0.05,
            fq_upper_bound = 1.1, comprehensive = False)

        if args.v_scale:
            v_scale = [ float(s) for s in args.v_scale.split(",") ]
        else:
            v_scale = None

        for i in range(len(layouts)): 

            if args.onlys and i not in [ int(s) for s in args.onlys.split(",") ]:
                # if it's not designated as one to be visualized,
                continue

            if len(layouts[i]) > 4:
                print(f"consensus for {i}")
                cons = consensus(layouts[i], hers, arrs, name = f"CS-{i}")
                snvs = var([ hers[li] for li, lk in layouts[i] ])
                snvs_read = var([cons], err_rate = 0.01, comprehensive = True)
                snvs_list = [
                    (v_all, "all"), (v_major, "major"),
                    (snvs, "local"), (snvs_read, "consensus"),
                    (None, "naive")]

                print(f"Total {len(snvs_read)} SNVs on this layout.")
                show_layout(layouts[i], hers, arrs,
                        cons_read = cons, snvs_list = snvs_list, v_scale = v_scale)

