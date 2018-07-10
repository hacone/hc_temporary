## This is for studying SNVs (mismatches) distribution in PacBio reads and short reads, over HOR-encoded regions or each monomer.
## Here is where HOR analysis and SNV analysis come together

# TODO: this is the most dirty code in this project. 
# TODO: rename this script

from collections import Counter
import matplotlib as mpl
mpl.use('SVG')
import matplotlib.pyplot as plt
import numpy as np
import re
import svgwrite
import sys
# for the other data structures
from EncodedRead import *
# for HOR_read
from HOR_segregation import *

def tsne_units(bitvec):

    # TODO; I'll do many types of analysis;
    # distance stats for units neighboring (k=2, 3, 4, ...), intra-read, inter-reads, random globally, in boxplot or in density est.
    # K means clusters and t-SNE embedding (all or each represents one read)

    kmeans = []
    for n_sites in [15, 27, 40]:
        from sklearn.cluster import KMeans
        n_clusters = 3
        kmeans += [ KMeans(n_clusters=n_clusters, random_state=0, verbose=1, n_jobs=-3).fit(bitvec[:,2:(2+n_sites)]) ]
        print(len(kmeans))

    # rest are skipped.
    return 


    for n_sites in [15, 27, 40]:
        from sklearn.manifold import TSNE
        #reduced = TSNE(n_components=2, random_state=0, verbose=1).fit_transform(bitvec[:,2:(2+n_sites)])
        #np.save(f"X_default_units.tsne_{n_sites}snv.npy", reduced)
        reduced = np.load(f"X_default_units.tsne_{n_sites}snv.npy")
        print("shape = " + f"{reduced.shape}")

        for i,j in [(0, 15), (1, 27), (2, 40)]:
            cls = kmeans[i].predict(bitvec[:,2:(2+j)])
            #plt.scatter(reduced[:, 0], reduced[:, 1], c=cls, s=6, alpha=0.5, edgecolors="none", cmap=mycm)
            plt.scatter(reduced[:, 0], reduced[:, 1], c=cls, s=2, alpha=0.4, edgecolors="none")
            plt.savefig(f"units_tsne_{n_sites}snv_{j}.svg")
            plt.close()

# NOTE: this is partly sorted by freq.
SNV_sites = [
        (2, 15, "T"), (2, 135, "A"), (2, 139, "C"), (3, 108, "C"), (3, 59, "G"), (3, 116, "C"), (3, 137, "G"), (4, 92, "T"),
        (4, 122, "C"), (4, 5, "G"), (5, 45, "G"), (5, 109, "G"), (5, 65, "T"), (6, 15, "C"), (6, 137, "C"), (6, 41, "G"),
        (6, 89, "C"), (7, 127, "A"), (8, 105, "A"), (8, 75, "C"), (9, 138, "A"), (9, 148, "G"), (9, 163, "T"), (9, 78, "G"),
        (10, 101, "G"), (10, 43, "A"), (11, 22, "C"), (6, 36, "C"), (10, 137, "G"), (2, 90, "A"), (4, 25, "C"), (11, 68, "G"),
        (2, 154, "G"), (4, 35, "G"), (11, 154, "C"), (5, 3, "G"), (11, 12, "G"), (3, 132, "T"), (11, 39, "C"), (9, 162, "A") ]

# TODO: this should be outside
def consecutive_units(hors):
    """ get consecutive (default) units in a read"""
    arrays = []
    for h in filter(lambda x: x[2] == "~", hors):
        if arrays and h[0] == arrays[-1][1] + 12:
            arrays[-1][1] = h[0]
        else:
            arrays += [ [h[0], h[0]] ]
    #return [ [i,j] for i, j in arrays if i < j ]
    return arrays if any([ i < j for i, j in arrays ]) else []

# TODO: take specifier of units? or can i use whole reps, ignoring some sites conveniently?
def bit_vector_repr(ers, sites, range_of_units = None):
    """ obtain bit vector representation as np.ndarray 
        param: 
            ers = HOR encoded reads
            sites = SNV sites as iterable of SNV (cf. defined as in EncodedRead)
        out:
            [rid, aid, hor_idx, SNV_states*]
    """

    if range_of_units:
        bv = []
        for ri, h in range_of_units:
            bv += [ri, -1, h]
            for ki, p, b in sites:
                bv += [1] if any([ (sp.pos, sp.base) == (p, b) for sp in ers[ri].mons[h+ki].monomer.snvs ]) else [0]
        return np.array(bv).reshape(len(range_of_units), 3+len(sites))

    bv, n_units = [], 0
    for ri, er in enumerate(ers):
        # TODO: consecutive units detection should not be here
        for ai, (k, l) in enumerate(consecutive_units(er.hors)): 
            for s in range(k, l+1, 12):
                #ris += [ri]
                #ais += [ai]
                bv += [ri, ai, s]
                n_units += 1
                for ki, p, b in sites:
                    # NOTE: seems big
                    if [ sp for sp in er.mons[s+ki].monomer.snvs if sp.pos == p and sp.base == b ]:
                        bv += [1]
                    else:
                        bv += [0]
    return np.array(bv).reshape(n_units, 3 + len(sites))

def cluster(hers):
    """ describe current cluster and proceed """

    # Take initial set of default units.
    bag_of_units = [ (ri, h) for ri, er in enumerate(hers) for h, _, t in er.hors if t == "~" ]
    n_units = len(bag_of_units)
    print(f"{n_units} units found in {len(hers)} reads.", flush=True)

    # Detect variants
    def detect_snvs(units):
        n_units = len(units)
        counter = Counter()
        for ri, i in units:
            for j in range(2, 12):
                counter.update([ (j, s.pos, s.base) for s in hers[ri].mons[i+j].monomer.snvs ])

        print("i\tk\tpos\tbase\tfreq\t%freq")
        for i, ((k, p, b), c) in enumerate(counter.most_common(8000)):
            if 0.05 < c/n_units and c/n_units < 0.8:
                print(f"{i}\t{k}\t{p}\t{b}\t{c}\t{100*c/n_units}%")

        return [ (k, p, b) for (k, p, b), c in counter.most_common(8000) if (0.05 < c/n_units) and (c/n_units < 0.8) ]


    snv_sites = detect_snvs(bag_of_units)
    print(f"{len(snv_sites)} SNV sites defined:")

    # TODO: organize below
    brep = bit_vector_repr(hers, snv_sites, range_of_units = bag_of_units)
    print(f"brep.shape = {brep.shape}")

    # define tentative clustering
    from sklearn.cluster import KMeans
    km_cls = KMeans(n_clusters=3, random_state=0, n_jobs=-3).fit(brep[:,3:]).predict(brep[:,3:])
    print("done clustering")

    # t-SNE
    if False:
        from sklearn.manifold import TSNE
        reduced = TSNE(n_components=2, random_state=0, verbose=0).fit_transform(brep[:,3:])
        print("done t-SNE")
        plt.scatter(reduced[:,0], reduced[:,1], c=km_cls, s=2, alpha=0.4, edgecolors="none")
        plt.savefig(f"units_tsne_{len(snv_sites)}snv_{len(bag_of_units)}units.svg")
        plt.close()
        # TODO: PCA
        from sklearn.decomposition import PCA
        reduced = PCA(n_components=2, random_state=0).fit_transform(brep[:,3:])
        print("done PCA")
        plt.scatter(reduced[:,0], reduced[:,1], c=km_cls, s=2, alpha=0.4, edgecolors="none")
        plt.savefig(f"units_pca_{len(snv_sites)}snv_{len(bag_of_units)}units.svg")
        plt.close()

    import pandas as pd
    df_units = pd.DataFrame(brep, columns = ["read_id", "array_id", "unit_idx"] + list(range(len(snv_sites)))).assign( cls = km_cls )

    outfile = "kmeans_cluster.desc.dat"
    print("mean, variance (original; later they should be plotted)")
    pd.DataFrame(df_units.mean()).T.to_csv(outfile, sep="\t", float_format = "%.3f", mode = "a")
    pd.DataFrame(df_units.var()).T.to_csv(outfile, sep="\t", float_format = "%.3f", mode = "a")
    print("count, mean, variance (clusters; later they should be plotted)")
    df_units.groupby(by = "cls").count().to_csv(outfile, sep="\t", float_format = "%.3f", mode = "a")
    df_units.groupby(by = "cls").mean().to_csv(outfile, sep="\t", float_format = "%.3f", mode = "a")
    df_units.groupby(by = "cls").var().to_csv(outfile, sep="\t", float_format = "%.3f", mode = "a")

    # TODO: loop for each cluster
    for cls_i in range(3):
        c0_units = [ (t.read_id, t.unit_idx) for t in df_units.itertuples() if t.cls == cls_i ]
        print("play with the cluster" + f" {cls_i}, which had {len(c0_units)} units.", flush = True)
        snv_sites = detect_snvs(c0_units)
        print(f"{len(snv_sites)} SNV sites defined:")
        # NOTE: rest are just copied to see how it works.
        brep = bit_vector_repr(hers, snv_sites, range_of_units = c0_units)
        km_cls = KMeans(n_clusters=3, random_state=0, n_jobs=-3).fit(brep[:,3:]).predict(brep[:,3:])

        # t-SNE
        from sklearn.manifold import TSNE
        reduced = TSNE(n_components=2, random_state=0, verbose=0).fit_transform(brep[:,3:])
        print("done t-SNE")
        plt.scatter(reduced[:,0], reduced[:,1], c=km_cls, s=2, alpha=0.4, edgecolors="none")
        plt.savefig(f"units_tsne_{len(snv_sites)}snv_{len(c0_units)}units.svg")
        plt.close()

        # TODO: PCA
        from sklearn.decomposition import PCA
        reduced = PCA(n_components=2, random_state=0).fit_transform(brep[:,3:])
        print("done PCA")
        plt.scatter(reduced[:,0], reduced[:,1], c=km_cls, s=2, alpha=0.4, edgecolors="none")
        plt.savefig(f"units_pca_{len(snv_sites)}snv_{len(c0_units)}units.svg")
        plt.close()

        _df_units = pd.DataFrame(brep, columns = ["read_id", "array_id", "unit_idx"] + list(range(len(snv_sites)))).assign( cls = km_cls )
        outfile = "kmeans_cluster.desc.dat"
        print(f"mean for {len(snv_sites)} sites, {len(c0_units)} units.", file = open(outfile, "a"))
        pd.DataFrame(_df_units.mean()).T.to_csv(outfile, sep="\t", float_format = "%.3f", mode = "a")
        print(f"variance for {len(snv_sites)} sites, {len(c0_units)} units.", file = open(outfile, "a"))
        pd.DataFrame(_df_units.var()).T.to_csv(outfile, sep="\t", float_format = "%.3f", mode = "a")
        print(f"count, mean, variance for 3 clusters defined on {len(snv_sites)} sites, {len(c0_units)} units.", file = open(outfile, "a"))
        _df_units.groupby(by = "cls").count().to_csv(outfile, sep="\t", float_format = "%.3f", mode = "a")
        _df_units.groupby(by = "cls").mean().to_csv(outfile, sep="\t", float_format = "%.3f", mode = "a")
        _df_units.groupby(by = "cls").var().to_csv(outfile, sep="\t", float_format = "%.3f", mode = "a")

# TODO: rename this
# TODO: this will perform iterative classifying. and summary picture at each stage.
def load_array_of_default(hers):

    # TODO: to be defined elsewhere # TODO: reasonably, read from tabulated file, instead of hardcoding
    s20 = [ (2, 135, "A"), (6, 15, "C"), (2, 15, "T"), (8, 105, "A"), (9, 138, "A"), (5, 65, "T"), (10, 43, "A"), (8, 75, "C"),
            (6, 41, "G"), (3, 116, "C"), (7, 127, "A"), (2, 139, "C"), (6, 89, "C"), (4, 5, "G"), (4, 122, "C")]
    s10 = [ (10, 101, "G"), (9, 163, "T"), (4, 92, "T"), (5, 45, "G"), (6, 137, "C"), (5, 109, "G"), (9, 78, "G"), (3, 108, "C"),
            (3, 59, "G"), (3, 137, "G"), (9, 148, "G"), (11, 22, "C")]
    s5 = [ (6, 36, "C"), (10, 137, "G"), (2, 90, "A"), (4, 25, "C"), (11, 68, "G"), (2, 154, "G"), (4, 35, "G"), (11, 154, "C"),
           (5, 3, "G"), (11, 12, "G"), (3, 132, "T"), (11, 39, "C"), (9, 162, "A")]

    ers = pickle.load(open(hers, "rb"))
    # pick the longest stretch of default unit
    n_sites = len(s5) + len(s10) + len(s20) # = 40

    # TODO: already extracted.
    #bitv_rep = np.zeros(len(ers) * n_sites).reshape(len(ers), n_sites)
    acc_bitv_rep = []
    idx_readid = []
    idx_arrid = []
    arr_i = 0
    for ers_i in range(len(ers)):
        er = ers[ers_i]
        arrays = []
        for h in filter(lambda x: x[2] == "~", er.hors):
            if arrays and h[0] == arrays[-1][1] + 12:
                arrays[-1][1] = h[0]
            else:
                arrays += [ [h[0], h[0]] ]
        if (not arrays) or all([ i == j for i,j in arrays ]):
            continue

        # TODO: NOTE: this part output the string representation of the SNV-compressed reads, which was very essential. But functionally this should be elsewhere.
        # TODO: string repr should be derived from bit vec rep.

        is_name_said = False
        for i,j in arrays:

            if not is_name_said:
                print(er.name)
                print(f"idx\t" + " ".join([ s[2] for s in s20 ]) + " | " + " ".join([ s[2] for s in s10 ]) + " | " + " ".join([ s[2] for s in s5 ]))
                is_name_said = True

            # TODO: unit size are given through param
            for s in range(i, j+1, 12):

                def bitv_rep(sites):
                    # TODO: logic is the same so maybe merged with str_rep()
                    bv = [ 0 for i in range(len(sites)) ]
                    for (l, (k, p, b)) in enumerate(sites):
                        if [ sp for sp in er.mons[s+k].monomer.snvs if sp.pos == p and sp.base == b ]:
                            bv[l] = 1
                    return bv

                acc_bitv_rep += bitv_rep(s20 + s10 + s5)
                idx_readid += [ers_i]
                idx_arrid += [arr_i]

                def str_rep(sites):
                    str_rep = [ "." for i in range(len(sites)) ]
                    for (l, (k, p, b)) in enumerate(sites):
                        if [ sp for sp in er.mons[s+k].monomer.snvs if sp.pos == p and sp.base == b ]:
                            str_rep[l] = "#"
                    return str_rep
                str_line = " | ".join([ " ".join(str_rep(s)) for s in [s20, s10, s5] ])
                #print(f"{er.mons[s].begin}\t{str_line}") # TODO: back these three print

            arr_i += 1
            #print("")

        #print("", flush=True)
        print(f"{ers_i} / {len(ers)} : {arr_i}", flush = True)
        
        # TODO: implement later
        #for i in range(len(er.hors)-1):
        #    if (er.hors[i][2] er.hors[i+1][2]):
                # make one

        # for i in [ h for h, _, t in er.hors if t == "~" ]:

    bv_data = np.array(acc_bitv_rep, dtype="int").reshape(len(idx_readid), n_sites)
    idx_col = np.append( np.array(idx_readid, dtype="int").reshape(len(idx_readid), 1),
                         np.array(idx_arrid, dtype="int").reshape(len(idx_readid), 1), axis = 1 )
    bv_data_withid = np.append(idx_col, bv_data, axis = 1)

    # NOTE: then you can play with these data.
    np.save("forty_bits_vector.npy", bv_data)
    np.save("forty_bits_vector_withid.npy", bv_data_withid)
    np.savetxt("forty_bits_vector.txt", bv_data, fmt="%d")
    np.savetxt("forty_bits_vector_withid.txt", bv_data_withid, fmt="%d")

# TODO: resolve all TODO. extract the logic where frequent SNVs are detected.
def draw_all_default(hers, us, something_for_short_read = None):
    """
        generate SVG showing mismatch/SNV distribution over default unit (12-mer, for X).
        later this should be abstract out. args are:
        hers = HOR encoded reads
        us = unit size (can be inferred from hers[i].hors?)
    """

    # open long read # TODO: move this part to there.
    ers = pickle.load(open(hers, "rb"))
    n_units = 0
    snv_counters = [ Counter() for i in range(us) ]
    for er in ers:
        for i in [ h for h, _, t in er.hors if t == "~" ]:
            for j in range(us):
                snv_counters[j].update(er.mons[i+j].monomer.snvs)
            n_units += 1

    # NOTE: temporal TODO: move elsewhere # can be oneliner?
    #print("[ ")
    #for j in range(2, us):
    #    for k, v in [ (k, v) for k, v in snv_counters[j].items() if k != "Monomer" and (v / n_units > 0.05) ]:
    #        print(f"({j}, {k.pos}, \"{k.base}\"),\t{v/n_units}")
    #print("]")

    # open short read # TODO: short read pickle file should not be hardcoded
    vs = pickle.load(open("./variant_sites_50_0.01.pickle", "rb"))
    sr_snv = { re.sub("\.[0-9]+", "", k):v for k,v in list(vs.items()) }

    # draw # TODO: ideally output can be specified by args
    dwg = svgwrite.drawing.Drawing("./mms_X_default.svg", profile = "tiny")
    # TODO: legend for this color scheme?
    b2c = dict(A = "#F8766D", C = "#7CAE00", G = "#00BFC4", T = "#C77CFF")
    b2i = dict(A = 0, C = 1, G = 2, T = 3)

    xoff, yoff = 100, 50

    # X axis
    xx= dwg.g(id="x_axis")
    for i in [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
        for j in range(us):
            # TODO: scale factor (1000) should not be hardcoded.
            x =  xoff + (1000 * i / 100)
            xx.add(dwg.text(f"{i} %", insert = (x, yoff-10), font_size = "20"))
            xx.add(dwg.line((x, yoff + (750*j)), (x, yoff + (750*(j+1)) - 15), stroke = "#000000"))
            # for short read part
            x =  xoff + 1000 + (1000 * i / 100)
            xx.add(dwg.text(f"{i} %", insert = (x, yoff-10), font_size = "20"))
            xx.add(dwg.line((x, yoff + (750*j)), (x, yoff + (750*(j+1)) - 15), stroke = "#000000"))
    dwg.add(xx)

    for j in range(us):

        # mismatch set in pacbio long reads
        mms_pac = dwg.g(id=f"mms_pac_{j}")
        for s, c in snv_counters[j].items():
            print(f"{j}\t{s.pos}\t{s.base}\t{c}")
            mms_pac.add(dwg.rect(
                insert=(xoff, yoff + (750 * j) + (s.pos * 4) + b2i[s.base]), size=(1000 * c / n_units, 1), fill = b2c[s.base], fill_opacity = 0.8))
        dwg.add(mms_pac)

        # Y axis
        yx= dwg.g(id=f"y_axis_{j}")
        for i in [0, 50, 100, 150]:
            yx.add(dwg.text(f"{j}:{i}", insert = (40, yoff + (750*j) + (4*i)), font_size = "20"))
            yx.add(dwg.line((xoff-20, yoff + (750*j) + (4*i)), (xoff, yoff + (750*j) + (4*i)), stroke = "#000000"))
        dwg.add(yx)

    xoff = 1100
    # TODO: skipping 0, 1, is so special for X. I need recover these position somehow...
    for j in range(2, us):
        srj = sr_snv[f"horID_17.mon_{j}"]
        m_total = srj["Monomer"]
        mms_sr = dwg.g(id=f"mms_sr_{j}")
        for pb, c in srj.items():
            if pb == "Monomer":
                continue
            mms_sr.add(dwg.rect(
                insert=(xoff, yoff + (750 * j) + (pb[0] * 4) + b2i[pb[1]]), size=(1000 * c / m_total, 1), fill = b2c[pb[1]], fill_opacity = 0.8))
        dwg.add(mms_sr)

        # y axis for short read part
        yx= dwg.g(id=f"y_axis_sr_{j}")
        for i in [0, 50, 100, 150]:
            yx.add(dwg.text(f"{j}:{i}", insert = (xoff-60, yoff + (750*j) + (4*i)), font_size = "20"))
            yx.add(dwg.line((xoff-20, yoff + (750*j) + (4*i)), (xoff, yoff + (750*j) + (4*i)), stroke = "#000000"))
        dwg.add(yx)

    dwg.save()

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='Breakup encoded read based on the set of assigned monomers.')
    parser.add_argument('action', metavar='action', type=str, help='action to perform: plot, tmp-bitv, tsne, cluster, ...')
    parser.add_argument('--hor-reads', dest='hors', help='pickled hor-encoded long reads')
    parser.add_argument('--unit-size', dest='us', help='default unit size (currently not inferred)')
    parser.add_argument('--bitvec', dest='bitvec', help='bit vector of SNV status (with 2 index columns)')
    args = parser.parse_args()

    if args.action == "plot":
        assert args.hors, "give me HOR-encoded reads"
        assert args.us, "give me unit size of the default unit"
        draw_al_default(args.hors, args.us)
    elif args.action == "tmp-bitv":
        # TODO: this is temporal. analysis based on bit vector representation and generating string representation
        assert args.hors, "give me HOR-encoded reads"
        load_array_of_default(args.hors)
    elif args.action == "tsne":
        assert args.bitvec, "bit vector representation .npy is required"
        bitvec = np.load(args.bitvec)
        tsne_units(bitvec)
    elif args.action == "cluster":
        assert args.hors, "need HOR-encoded reads"
        hers = pickle.load(open(args.hors, "rb"))
        cluster(hers)
    else:
        assert None, "invalid action."

