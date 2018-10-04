## tHIS IS FOR STUDYing SNVs (mismatches) distribution in PacBio reads and short reads, over HOR-encoded regions or each monomer.
## Here is where HOR analysis and SNV analysis come together

# TODO: this is the most dirty code in this project. 
# TODO: rename this script?

from collections import Counter
import numpy as np
import re
import svgwrite
import sys

# NOTE: these 3 lines must be in this order.
import matplotlib as mpl
mpl.use('SVG')
import matplotlib.pyplot as plt

from EncodedRead import *
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

def bit_vector_repr(ers, sites, range_of_units):
    """ obtain bit vector representation as np.ndarray 
        missing units (with spacing of 23~24 mons) are represented as 0-masked units
        param: 
            ers = HOR encoded reads (load from pickle)
            sites = SNV sites as iterable of SNV (cf. defined as in EncodedRead.py)
            range_of_units = (possibly a subset of) units which should be encoded into bit-vec repr.
        out:
            [rid, aid, hor_idx, SNV_states*]
    """
    bv, n_units = [], 0
    last_ri, last_s = -1, -1
    for ri, h in range_of_units:
        # NOTE: here I filling up gaps due to 1 or 2 missing monomer assignment.
        if ri == last_ri and (h == last_s + 22 or h == last_s + 23 or h == last_s + 24):
            bv += [ri, -1, last_s + 12] + [ 0 for i in sites ]
            n_units += 1
        elif ri == last_ri and (h == last_s + 33 or h == last_s + 34 or h == last_s + 35):
            bv += [ri, -1, last_s + 12] + [ 0 for i in sites ]
            bv += [ri, -1, last_s + 24] + [ 0 for i in sites ] # 2nd units
            n_units += 2
        bv += [ri, -1, h]
        for ki, p, b, f in sites:
            bv += [1] if any([ (sp.pos, sp.base) == (p, b) for sp in ers[ri].mons[h+ki].monomer.snvs ]) else [-1]
        last_ri, last_s = ri, h
    return np.array(bv).reshape(len(range_of_units) + n_units, 3+len(sites))

def get_e67(pi):
    """
    return a log odds ratio matrix of shape = (4, n_snv_sites)
    this depends on accuracy estimates p, q, and (local) SNV distribution in a random region, p1.
    p, q are currently scalars and can be defined globally. p1 is matrix of shape = (1, n_snv_sites).
    """

    # variables for alignment score calc
    n = pi.shape[1]
    p, q = 0.80, 0.90
    # p, q = 0.7, 0.7

    e_tilde = np.zeros(4*n).reshape(4, n)
    e_tilde[0,:] = (1-p)*(1-p)*pi + q*q*(1-pi) # 00
    e_tilde[1,:] = p*(1-p)*pi + (1-q)*q*(1-pi) # 01
    e_tilde[2,:] = p*(1-p)*pi + (1-q)*q*(1-pi) # 10 = e_tilde[1,:]
    e_tilde[3,:] = p*p*pi + (1-q)*(1-q)*(1-pi) # 11

    e_hat = np.zeros(4*n).reshape(4, n)
    e_hat[0,:] = ((1-p)*pi + q*(1-pi)) * ((1-p)*pi + q*(1-pi)) # 00
    e_hat[1,:] = p*(1-p)*pi*pi + (1-q)*q*(1-pi)*(1-pi) + (p*q + (1-p)*(1-q))*pi*(1-pi) # 01
    e_hat[2,:] = p*(1-p)*pi*pi + (1-q)*q*(1-pi)*(1-pi) + (p*q + (1-p)*(1-q))*pi*(1-pi) # 10 = e_hat[1,:]
    e_hat[3,:] = (p*pi + (1-q)*(1-pi)) * (p*pi + (1-q)*(1-pi)) # 11

    e67 = np.divide(e_tilde, e_hat)
    return e67

def cluster(hers):
    """ describe current cluster and proceed """

    # Take initial set of default units.
    bag_of_units = [ (ri, h) for ri, er in enumerate(hers) for h, _, t in er.hors if t == "~" ]
    n_units = len(bag_of_units)
    print(f"{n_units} units found in {len(hers)} reads.", flush=True)

    # Detect variants
    def detect_snvs(units):
        """ return: [(monomer index in unit, position in monomer, alt. base, relative frequency)] """
        n_units = len(units)
        counter = Counter()
        for ri, i in units:
            for j in range(2, 12): # I'm skipping the first 2 monomers, where assignment is always wrong
                counter.update([ (j, s.pos, s.base) for s in hers[ri].mons[i+j].monomer.snvs ])

        print("id\tmon_id\tpos\talt.base\tfreq\t%freq")
        for i, ((k, p, b), c) in enumerate(counter.most_common(1000)):
            if 0.05 < c/n_units and c/n_units < 0.8:
                print(f"{i}\t{k}\t{p}\t{b}\t{c}\t{100*c/n_units:.2f}%")

        # TODO: redundant calculation.
        return [ (k, p, b, c/n_units) for (k, p, b), c in counter.most_common(1000) if (0.05 < c/n_units) and (c/n_units < 0.8) ]

    def valid_regions(bitvec):
        """ get valid consecutive regions after filling up 22/23/24-mons or 33/34/35-mons gaps """
        regs = []
        b, last_mi, last_ri = 0, bitvec[0,2], bitvec[0,0]
        for i in range(bitvec.shape[0]):
            if (bitvec[i,0] != last_ri) or (bitvec[i,2] > last_mi + 12):
                regs += [(b, i)]
                b = i
            last_mi, last_ri = bitvec[i,2], bitvec[i,0]
        regs += [(b, bitvec.shape[0])]
        return regs

    snv_sites = detect_snvs(bag_of_units)
    print(f"{len(snv_sites)} SNV sites defined.")

    # TODO: organize below
    brep = bit_vector_repr(hers, snv_sites, range_of_units = bag_of_units)

    regs = sorted(valid_regions(brep), key = lambda x: x[0] - x[1])
    print(f"{len(regs)} consecutive regions found")

    # define tentative clustering
    from sklearn.cluster import KMeans
    # clustering without masked units.
    km = KMeans(n_clusters=3, random_state=0, n_jobs=-3).fit( brep[np.any(brep[:,3:] != 0, axis = 1), 3:] )
    km_cls = km.predict(brep[:,3:])
    for i in range(1, brep.shape[0]):
        if brep[i,3] == 0: # if masked
            km_cls[i] = km_cls[i-1]
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
    df_units = pd.DataFrame(brep, columns = ["read_id", "array_id", "unit_idx"] + list(range(len(snv_sites))))
    df_units = df_units.assign( cls = km_cls )
    df_units = df_units.assign( is_hom = df_units.groupby(by = "read_id").cls.transform(lambda x: x.std() == 0 or x.count() < 2) )

    # NOTE: better text representation.
    #print("uid\trid\tcls\tmix\t" + "\t".join([ f"{i}" for i in range(brep.shape[1]-3) ])) # header
    #print("rid\tuidx\tcid\tis_mix\t" + f"{brep.shape[1]-3}_bits" ) # header
    last_rid = 0

    min_units = 12
    long_regs = list(filter(lambda x: x[1] - x[0] >= min_units, regs))
    print(f"calculating alignments for {len(long_regs)} regions with >={min_units} units...")

    def diff_bv(bv1, bv2):
        """ simple correlation of -1/0/1-bit vector """
        assert bv1.shape == bv2.shape, "wrong shapes"
        m = np.multiply(bv1, bv2)
        a = np.sum(np.absolute(m))
        s = (1 + np.sum(m) / (a+1)) / 2
        return s

    def diff_bv2(bv1, bv2):
        """ n11 / (n10+n01+n11) """
        assert bv1.shape == bv2.shape, "wrong shapes"

        m = np.multiply(bv1, bv2) # 1: match, 0: masked, -1: mismatch
        mt = np.multiply((m + 1), m) / 2 # 1: match, 0: masked/mismatch
        n_m = np.sum(np.multiply(m, m) - mt) # n01 or n10
        n11 = np.sum(np.multiply(mt, (bv1 + 1) / 2)) # n11
        return n11 / (1 + n11 + n_m)

    # global pi # NOTE: to be refined
    pi = np.array([ f for k,p,b,f in snv_sites ]).reshape(1, len(snv_sites))
    e67 = get_e67(pi)

    #def diff67(bv1, bv2, e67): # NOTE: I don't need this as e67 is constant, since pi is constant now.
    def diff67(bv1, bv2):
        assert bv1.shape == bv2.shape, "wrong shapes"

        m = np.multiply(bv1, bv2)
        mt = np.multiply((m + 1), m) / 2

        x = np.zeros(len(snv_sites)*4).reshape(4, len(snv_sites))
        x[0,:] = np.sum(np.multiply(mt, (1 - bv1) / 2), axis=0) # n00
        x[1,:] = np.sum((mt - 1) / 2, axis=0) # n01 or n10 # TODO FIXME bug when masked
        # x[2,:] = np.sum((mt - 1) / 2, axis=0) # n01 or n10
        x[3,:] = np.sum(np.multiply(mt, (bv1 + 1) / 2), axis=0) # n11
        return np.sum(np.multiply(x, np.log(e67)))

    # Threshold for reporting alignments
    ms_t = 0.6 # for primary score. correlation (as in diff_bv) or identity (as in diff_bv2)
    ms2_t = 5.0

    targets = dict() # this will accumulates alignment results. # TODO: rename. utilize.

    #for ri in range(0, 5): # for testing
    for ri in range(len(long_regs)):
        r = long_regs[ri]
        targets[ri] = dict()

        #for si in range(ri+1, len(long_regs)): # avoids redundant calc but makes stats from reads difficult to compare.
        target_sis = [ _si for _si in range(len(long_regs)) if _si != ri ]
        for si in target_sis:
            targets[ri][si] = [] # this accumulate results of alignments
            s = long_regs[si]

            l_r = r[1]-r[0]
            l_s = s[1]-s[0]

            # dovetail.
            #for i in range(2, l_s):
            for i in range(2, min([l_s, l_r])):
                rs, re = r[0], r[0]+i
                qs, qe = s[1]-i, s[1]
                # NOTE: use diff_bv or diff_bv2
                ms, ms2 = diff_bv2(brep[rs:re,3:], brep[qs:qe,3:]), diff67(brep[rs:re,3:], brep[qs:qe,3:])
                targets[ri][si] += [ ri, si, rs, re, qs, qe, ms, ms2, ms2/i ]

            # contained/containing
            for i in range(0, abs(l_r - l_s) + 1):
                if l_r > l_s:
                    rs, re = r[0]+i, r[0]+i+l_s
                    qs, qe = s[0], s[1]
                else:
                    rs, re = r[0], r[1]
                    qs, qe = s[1]-i-l_r, s[1]-i
                ms, ms2 = diff_bv2(brep[rs:re,3:], brep[qs:qe,3:]), diff67(brep[rs:re,3:], brep[qs:qe,3:])
                targets[ri][si] += [ ri, si, rs, re, qs, qe, ms, ms2, ms2/l_s ]

            # dovetail
            for i in range(min([l_s, l_r])-1, 1, -1):
                rs, re = r[1]-i, r[1]
                qs, qe = s[0], s[0]+i
                ms, ms2 = diff_bv2(brep[rs:re,3:], brep[qs:qe,3:]), diff67(brep[rs:re,3:], brep[qs:qe,3:])
                targets[ri][si] += [ ri, si, rs, re, qs, qe, ms, ms2, ms2/i ]

            targets[ri][si] = np.array(targets[ri][si]).reshape(l_s+l_r-3, 9)

        print(f"\nRef. region:\t{ri}\t{hers[brep[r[0],0]].name}")
        print(f"read_id\treg_id\treg_b\treg_e\treg_l")
        print(f"{brep[r[0],0]}\t{ri}\t{r[0]}\t{r[1]}\t{r[1]-r[0]}")

        print("\nStructure:")
        for i, t in enumerate(df_units.iloc[r[0]:r[1],:].itertuples()):
            line = f"{r[0]+i}\t{t.read_id}\t{t.unit_idx}\t{t.cls}\t" \
                + (".\t" if t.is_hom else "M\t") \
                + "".join([ ("X" if e == 1 else ("o" if e == 0 else "-")) for e in brep[r[0]+i,3:] ])
            print(line)

        print("\nDistribution of Correlation (Identity) (%):")
        for i, m in enumerate( sorted([ np.max(targets[ri][x], axis=0)[6] for x in target_sis ], key = lambda x: -x)[:30] ):
            print(f"{i+1}\t{100*m:.2f}")

        # for the best 25 target reads
        n_shown = 0
        for si in sorted(target_sis, key=lambda x: -1 * np.max(targets[ri][x], axis=0)[6])[:25]:
        #for si in sorted(list(range(ri+1, ri+11)), key=lambda x: -1 * np.max(targets[ri][x], axis=0)[6])[:50]:

            t = targets[ri][si]

            # for the best 1 alignments configuration
            rank = sorted([ (i, t[i,6], t[i,8]) for i in range(t.shape[0]) ], key = lambda x: -x[1])
            if n_shown < 5 or t[rank[0][0], 6] > ms_t:

                n_shown += 1
                ii = rank[0][0]
                print(f"\n\nTo:\t{si}\t{hers[brep[int(t[ii,4]),0]].name}")

                for i,c,r in rank[:5]:
                    print(f"{ri}\t{si}\t{int(t[i,2])}\t{int(t[i,3])}\t{int(t[i,4])}\t{int(t[i,5])}\t{100*t[i,6]:.2f}\t{t[i,7]:.2f}\t{t[i,8]:.2f}")

                print("\nAlignments:\n")
                rs, re = int(t[ii,2]), int(t[ii,3])
                qs, qe = int(t[ii,4]), int(t[ii,5])
                r0, r1 = long_regs[ri]
                q0, q1 = long_regs[si]
                if r0 < rs:
                    for i, tu in enumerate(df_units.iloc[r0:rs,:].itertuples()):
    if 0 < rs:
                        line = f"{r0+i}\t{tu.read_id}\t{tu.unit_idx}\t{tu.cls}\t" \
                            + (".\t" if tu.is_hom else "M\t") \
                            + "".join([ ("X" if e == 1 else ("o" if e == 0 else "-")) for e in brep[r0+i,3:] ]) \
                            + f"\t*\t*\t*\t*\t*"
                        print(line)
                elif q0 < qs:
                    for i, tu in enumerate(df_units.iloc[q0:qs,:].itertuples()):
                        line = f"*\t*\t*\t*\t*\t" \
                            + "".join([ ("X" if e == 1 else ("o" if e == 0 else "-")) for e in brep[q0+i,3:] ]) \
                            + f"\t{q0+i}\t{tu.read_id}\t{tu.unit_idx}\t{tu.cls}\t" \
                            + ("." if tu.is_hom else "M")
                        print(line)

                def str_pair(e, d):
                    # TODO: shorten this further? we may have dict somewhere?
                    if e == 1:
                        return "#" if d == 1 else ("." if d == -1 else "0")
                    if e == -1:
                        return "'" if d == 1 else (" " if d == -1 else "0")
                    if e == 0:
                        #return "0" if d == 1 else ("0" if d == -1 else "0")
                        return "0"

                for i, tu in enumerate(df_units.iloc[rs:re,:].itertuples()):
                    qu = df_units.iloc[qs+i,:]
                    line = f"{rs+i}\t{tu.read_id}\t{tu.unit_idx}\t{tu.cls}\t" \
                        + (".\t" if tu.is_hom else "M\t") \
                        + "".join([ str_pair(e, d) for e,d in zip(brep[rs+i,3:], brep[qs+i,3:]) ]) \
                        + f"\t{qs+i}\t{qu.read_id}\t{qu.unit_idx}\t{qu.cls}\t" \
                        + ("." if qu.is_hom else "M")
                    print(line)

                if re < r1:
                    for i, tu in enumerate(df_units.iloc[re:r1,:].itertuples()):
                        line = f"{re+i}\t{tu.read_id}\t{tu.unit_idx}\t{tu.cls}\t" \
                            + (".\t" if tu.is_hom else "M\t") \
                            + "".join([ ("X" if e == 1 else ("o" if e == 0 else "-")) for e in brep[re+i,3:] ]) \
                            + f"\t*\t*\t*\t*\t*"
                        print(line)
                elif qe < q1:
                    for i, tu in enumerate(df_units.iloc[qe:q1,:].itertuples()):
                        line = f"*\t*\t*\t*\t*\t" \
                            + "".join([ ("X" if e == 1 else ("o" if e == 0 else "-")) for e in brep[qe+i,3:] ]) \
                            + f"\t{qe+i}\t{tu.read_id}\t{tu.unit_idx}\t{tu.cls}\t" \
                            + ("." if tu.is_hom else "M")
                        print(line)

    return

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

        if False:
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
    parser = argparse.ArgumentParser(description='SNV distribution') # TODO: explain
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

