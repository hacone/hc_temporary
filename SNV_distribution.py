## tHIS IS FOR STUDYing SNVs (mismatches) distribution in PacBio reads and short reads, over HOR-encoded regions or each monomer.
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

def bit_vector_repr(ers, sites, range_of_units = None):
    """ obtain bit vector representation as np.ndarray 
        missing units (with spacing of 23~24 mons) are represented as 0-masked units
        param: 
            ers = HOR encoded reads
            sites = SNV sites as iterable of SNV (cf. defined as in EncodedRead)
        out:
            [rid, aid, hor_idx, SNV_states*]
    """

    if range_of_units:
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

    else: # if range_of_units is not specified 
        bv, n_units = [], 0
        for ri, er in enumerate(ers):
            last_s = -100 # dummy
            # TODO: consecutive units detection should not be here
            for ai, (k, l) in enumerate(consecutive_units(er.hors)): 
                if k == last_s + 23:
                    bv += [ri, ai-1, last_s + 12] + [ 0 for i in sites ]
                    n_units += 1
                elif k == last_s + 24:
                    bv += [ri, ai-1, last_s + 12] + [ 0 for i in sites ]
                    n_units += 1
                for s in range(k, l+1, 12):
                    bv += [ri, ai, s]
                    n_units += 1
                    last_s = s
                    for ki, p, b, f in sites:
                        # NOTE: seems big
                        if [ sp for sp in er.mons[s+ki].monomer.snvs if sp.pos == p and sp.base == b ]:
                            bv += [1]
                        else:
                            bv += [-1]
        return np.array(bv).reshape(n_units, 3 + len(sites))


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

    regs = valid_regions(brep)
    print(f"{len(regs)} consecutive regions found")

    # define tentative clustering
    from sklearn.cluster import KMeans
    km = KMeans(n_clusters=3, random_state=0, n_jobs=-3).fit( brep[np.any(brep[:,3:] != 0, axis = 1), 3:] )
    km_cls = km.predict(brep[:,3:])
    for i in range(1, brep.shape[0]):
        if brep[i,3] == 0: # if masked
            km_cls[i] = km_cls[i-1]
    print("done clustering")

    # variables for alignment score calc
    p, q = 0.80, 0.90
    # p, q = 0.7, 0.7

    # global pi # NOTE: to be refined
    pi = np.array([ f for k,p,b,f in snv_sites ]).reshape(1, len(snv_sites))
    print(f"pi = {pi}")

    e_tilde = np.zeros(4*len(snv_sites)).reshape(4, len(snv_sites))
    e_tilde[0,:] = (1-p)*(1-p)*pi + q*q*(1-pi) # 00
    e_tilde[1,:] = p*(1-p)*pi + (1-q)*q*(1-pi) # 01
    e_tilde[2,:] = p*(1-p)*pi + (1-q)*q*(1-pi) # 10
    e_tilde[3,:] = p*p*pi + (1-q)*(1-q)*(1-pi) # 11

    e_hat = np.zeros(4*len(snv_sites)).reshape(4, len(snv_sites))
    e_hat[0,:] = ((1-p)*pi + q*(1-pi)) * ((1-p)*pi + q*(1-pi)) # 00
    e_hat[1,:] = p*(1-p)*pi*pi + (1-q)*q*(1-pi)*(1-pi) + (p*q + (1-p)*(1-q))*pi*(1-pi) # 01
    e_hat[2,:] = p*(1-p)*pi*pi + (1-q)*q*(1-pi)*(1-pi) + (p*q + (1-p)*(1-q))*pi*(1-pi) # 10
    e_hat[3,:] = (p*pi + (1-q)*(1-pi)) * (p*pi + (1-q)*(1-pi)) # 11

    e67 = np.divide(e_tilde, e_hat)
    print("e67 = ")
    print(e67)

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
    print("rid\tuidx\tcid\tis_mix\t" + f"{brep.shape[1]-3}_bits" ) # header
    last_rid = 0

    # TODO: organize here
    regs = sorted(regs, key = lambda x: x[0] - x[1])

    if False:
        for ii, jj in regs:
            print(f"L={jj-ii}")
            for i, t in enumerate(df_units.iloc[ii:jj,:].itertuples()):
                line = f"{t.read_id}\t{t.unit_idx}\t{t.cls}\t" \
                    + (".\t" if t.is_hom else "M\t") \
                    + "".join([ ("X" if e == 1 else ("o" if e == 0 else "-")) for e in brep[ii+i,3:] ])
                print(line)
            print("")

    # TODO test some alignments here!
    long_regs = list(filter(lambda x: x[1] - x[0] > 11, regs))
    print("calculating alignments for {len(long_regs)} regions...")

    def diff_bv(bv1, bv2):
        """ simple correlation of -1/0/1-bit vector """
        assert bv1.shape == bv2.shape, "wrong shapes"
        m = np.multiply(bv1, bv2)
        a = np.sum(np.absolute(m))
        s = (1 + np.sum(m) / a) / 2
        return s

    def diff67(bv1, bv2):
        assert bv1.shape == bv2.shape, "wrong shapes"

        m = np.multiply(bv1, bv2)
        mt = np.multiply((m + 1), m) / 2

        x = np.zeros(len(snv_sites)*4).reshape(4, len(snv_sites))
        x[0,:] = np.sum(np.multiply(mt, (1 - bv1) / 2), axis=0) # n00
        x[1,:] = np.sum((mt - 1) / 2, axis=0) # n01 or n10
        # x[2,:] = np.sum((mt - 1) / 2, axis=0) # n01 or n10
        x[3,:] = np.sum(np.multiply(mt, (bv1 + 1) / 2), axis=0) # n11

        return np.sum(np.multiply(x, np.log(e67)))

    ms_t = 0.9
    ms2_t = 0.0
    for r in long_regs[:20]:
        print(f"\n--->>> --->>> r = {r}")
        for s in list(filter(lambda x: (x[1]-x[0] <= r[1]-r[0]) and (x != r), long_regs)):
            print(f"\n-->> r = {r}; s = {s}")
            show = False

            l_r = r[1]-r[0]
            l_s = s[1]-s[0]

            # dovetail.
            for i in range(2, l_s):
                ms = diff_bv(brep[r[0]:r[0]+i,3:], brep[s[1]-i:s[1],3:])
                ms2 = diff67(brep[r[0]:r[0]+i,3:], brep[s[1]-i:s[1],3:])
                if (ms > ms_t and ms2/i > ms2_t and i > 2):
                    show = True
                    print(f"d-\t{ (r[0], r[0]+i) } ~ { (s[1]-i, s[1]) }\t{i}\t{ms:.3f}\t{ms2:.2f}\t{ms2/i:.2f}")
                #r[0], r[0]+i
                #s[1]-i, s[1]
            # contained
            for i in range(0, l_r - l_s + 1):
                ms = diff_bv(brep[r[0]+i:r[0]+i+l_s,3:], brep[s[0]:s[1],3:])
                ms2 = diff67(brep[r[0]+i:r[0]+i+l_s,3:], brep[s[0]:s[1],3:])
                if (ms > ms_t and ms2/l_s > ms2_t and l_s > 2):
                    show = True
                    print(f"ct\t{ (r[0]+i, r[0]+i+l_s) } ~ { (s[0], s[1]) }\t{l_s}\t{ms:.3f}\t{ms2:.2f}\t{ms2/l_s:.2f}")
                #(r[0]+i), (r[0]+i)+l_s = r[1]
                #s[0], s[1]
            # dovetail
            #for i in range(l_s-1, 1, -1):
            for i in range(l_s-1, 4, -1):
                ms = diff_bv(brep[r[1]-i:r[1],3:], brep[s[0]:s[0]+i,3:])
                ms2 = diff67(brep[r[1]-i:r[1],3:], brep[s[0]:s[0]+i,3:])
                if (ms > ms_t and ms2/i > ms2_t and i > 3):
                    show = True
                    print(f"d+\t{ (r[1]-i, r[1]) } ~ { (s[0], s[0]+i) }\t{i}\t{ms:.3f}\t{ms2:.2f}\t{ms2/i:.2f}")
                #r[1] - i, r[1]
                #s[0], s[0] + i

            if show:
                print(f"L={r[1]-r[0]}")
                for i, t in enumerate(df_units.iloc[r[0]:r[1],:].itertuples()):
                    line = f"{r[0]+i}\t{t.read_id}\t{t.unit_idx}\t{t.cls}\t" \
                        + (".\t" if t.is_hom else "M\t") \
                        + "".join([ ("X" if e == 1 else ("o" if e == 0 else "-")) for e in brep[r[0]+i,3:] ])
                    print(line)
                print("------------------------------")
                print(f"L={s[1]-s[0]}")
                for i, t in enumerate(df_units.iloc[s[0]:s[1]:].itertuples()):
                    line = f"{s[0]+i}\t{t.read_id}\t{t.unit_idx}\t{t.cls}\t" \
                        + (".\t" if t.is_hom else "M\t") \
                        + "".join([ ("X" if e == 1 else ("o" if e == 0 else "-")) for e in brep[s[0]+i,3:] ])
                    print(line)
                print("------------------------------")

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

def work_assembly(bitvec):
    """ take bitvector with rid, uidx, aid. perform something. """
    # TODO: what to do with missing unit?
    p = 0.85
    q = 0.95
    e00 = log ( 0.5*(1-p)*(1-p) + 0.5*q*q ) / ( 2*0.5*0.5*q*(1-p))
    e11 = log ( 0.5*(1-q)*(1-q) + 0.5*p*p ) / ( 2*0.5*0.5*p*(1-q))
    e01 = log ( 0.5*p*(1-p) + 0.5*q*(1-q) ) / ( 0.5*0.5*(p*q + (1-p)*(1-q)))
    def diff(i,j,k):
        # TODO: write up
        A = bitvec[i:i+k,3:]
        B = bitvec[j:j+k,3:] 
        S = (1-A)*(1-B)*e00 + (1-A)*B*e01 + A*(1-B)*e01 + A*B*e11
        score = sum(S)

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

