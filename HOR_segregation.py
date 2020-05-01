# segregate/encode reads by HOR patterns contianed in it.
from collections import Counter 
from EncodedRead import *
import numpy as np
import os
import pickle
import re
import functools
from more_itertools import windowed

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import svgwrite

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
# import seaborn as sns

# NOTE: here's data structure for HOR encoded reads (cf. EncodedRead.py for other lower-level representation)
# hors := (index, size, symbol)
HOR_Read = namedtuple("HOR_Read", ("name", "mons", "hors", "length", "ori"))

def load_encoded_reads(pickles, n_max_reads = None):
    reads = []
    for picklefile in [ l.strip() for l in open(pickles, "r").readlines() ]:
        rds = pickle.load(open(picklefile, "rb"))
        reads += [ r for r in rds if len(r.mons) > 29 ]
        # print(f"{len(reads)} reads found... " + "loaded " + picklefile, flush = True)
        if n_max_reads and (len(reads) > n_max_reads):
            break
    return reads

def monomers_in_reads(reads, ref = "d0.fasta.fai"):
    """
    counts occurences of each reference monomer in each read.
    returns ndarray of shape (nreads, nmons)
    """

    # TODO: get this from param
    l = np.loadtxt(ref, dtype = "U20", delimiter = "\t", usecols = (0))
    mon_to_id = { i:n for n, i in enumerate(l) }
    nmons = len(mon_to_id)
    occ = np.zeros(len(reads)*nmons).reshape(len(reads), nmons)
    for i, r in enumerate(reads):
        for m in r.mons:
            occ[i, mon_to_id[m.monomer.name]] += 1
    return occ

def cluster_reads(occ, n_clusters = 40):
    """
    perform clustering of reads based on monomer sharing.
    an obtained model is returned and optionally made persistent to specified file by pickling.
    """
    # normalize occ appropriately
    occ_n = occ.copy()
    for i in range(occ.shape[0]):
        occ_n[i] = occ[i] / sum([ occ[i,j] for j in range(occ.shape[1]) ])

    ## NOTE: may i use the subset??
    from sklearn.cluster import KMeans
    kmeans = KMeans(n_clusters=n_clusters, random_state=0, verbose=1, n_jobs=-3).fit(occ_n)
    return kmeans.predict(occ_n)


def segregate_reads(pickles, outdir):
    # pickles = [ pickles_dir + path for path in os.listdir(pickles_dir) if path.endswith('.ori.pickle') ]
    reads = load_encoded_reads(pickles)
    cls = cluster_reads(monomers_in_reads(reads))
    cls_reads = [ [] for i in range(40) ]

    for i, r in enumerate(reads):
        cls_reads[cls[i]] += [r]
    print("encoded read list is segregated, writing out")

    for i, c in enumerate(cls_reads):
        print(f"C-{i+1} has {len(c)} reads")
        pickle.dump(c, open(f"{outdir}/C_{i+1}.pickle", "wb"))

def print_reads(pkl):
    reads = pickle.load(open(pkl, "rb"))
    for r in sorted(reads, key=lambda x: -len(x.mons)):
        print(f"\n{r.name}\t{len(r.mons)}\t{r.length}")
        last_end = 0
        for m in r.mons:
            print(f"{m.monomer.name}\t{m.begin}\t{m.end}\t{m.begin-last_end}\t{len(m.monomer.snvs)}\t{m.ori}")
            last_end = m.end

def type_stats(pkl):
    reads = pickle.load(open(pkl, "rb"))
    print("\t".join([s for s in "JDWMRxAB123450"]))
    for r in sorted(reads, key=lambda x: -len(x.mons)):
        total = len(r.mons)
        # monomer type J, D, W, M, R, or x
        supfam = Counter([ m.monomer.name[2] for m in r.mons ])
        # archetyp A or B
        arch = Counter([ m.monomer.name[1] for m in r.mons ])
        # num 1-5, or 0
        monnum = Counter([ m.monomer.name[3] for m in r.mons ])
        print("\t".join([ f"{supfam[f]}" for f in "JDWMRx" ]) + "\t" +\
              "\t".join([ f"{arch[t]}" for t in "AB" ]) + "\t" +\
              "\t".join([ f"{monnum[n]}" for n in "123450" ] + [f"{total}"]))

# NOTE: currently not visible from the menu
def draw_cluster(pickles, outfile, precomputed = False):
    """ t-SNE / PCA embedding of reads based on monomer sharing to visualize the clusters. """

    n_clusters = 40
    reads = load_encoded_reads(pickles)
    # NOTE: use 20k reads with most monomers 
    occ = monomers_in_reads(sorted(reads,
        key=lambda x: -len(x.mons))[:20000])

    l = np.loadtxt("d0.fasta.fai", dtype = "U20", delimiter = "\t", usecols = (0))
    mon_to_id = { n: i for i, n in enumerate(l) }
    id_to_mon = { i: n.replace("horID_", "").replace(".mon_", "-")
            for i, n in enumerate(l) }
    nmons = len(mon_to_id)

    if not precomputed:
        # Normalize
        occ_n = occ.copy()
        for i in range(occ.shape[0]):
            occ_n[i] = occ[i] / sum([ occ[i,j] for j in range(occ.shape[1]) ])
        tsne_red = TSNE(n_components=2, random_state=0, verbose=1).fit_transform(occ_n)
        np.save(f"{outfile}.tsne.npy", tsne_red)
        print("tsne done")
        pca_red = PCA(n_components=2, random_state=0).fit_transform(occ_n)
        print("pca done")
        np.save(f"{outfile}.pca.npy", pca_red)
    else:
        tsne_red = np.load(f"{outfile}.tsne.npy")
        pca_red = np.load(f"{outfile}.pca.npy")

    # calculate clusters (n_clusters = 40)
    occ_all = monomers_in_reads(sorted(reads,
        key=lambda x: -len(x.mons)))
    cls = cluster_reads(occ_all, n_clusters = n_clusters)[:20000]

    # mycm = plt.get_cmap("tab20b") + plt.get_cmap("tab20c") # TODO: how can i concatenate?
    print("Cluster\t" + "\t".join([
        f"{id_to_mon[i]}" for i in range(len(id_to_mon)) ]),
        file=open(f"{outfile}.tsv", "a"))

    for c in range(n_clusters):
        print(f"C{c}\t" + "\t".join([
            f"{d}" for d in np.sum(occ[[ cls[i] == c for i in range(20000) ], :], axis=0)]),
            file=open(f"{outfile}.tsv", "a"))

    plt.scatter(
            tsne_red[:, 0], tsne_red[:, 1], c=cls,
            s=6, alpha=0.5, edgecolors="none", cmap=plt.get_cmap("tab20b"))
    for c in range(n_clusters):
        xs = [ tsne_red[i, 0] for i in range(20000) if cls[i] == c ]
        ys = [ tsne_red[i, 1] for i in range(20000) if cls[i] == c ]
        plt.text(sum(xs)/len(xs), sum(ys)/len(ys), f"{c}",
            horizontalalignment = "center", verticalalignment = "center")
    plt.savefig(f"{outfile}.tsne.svg")
    plt.close()

    plt.scatter(
            pca_red[:, 0], pca_red[:, 1], c=cls,
            s=6, alpha=0.5, edgecolors="none", cmap=plt.get_cmap("tab20b"))
    for c in range(n_clusters):
        xs = [ pca_red[i, 0] for i in range(20000) if cls[i] == c ]
        ys = [ pca_red[i, 1] for i in range(20000) if cls[i] == c ]
        plt.text(sum(xs)/len(xs), sum(ys)/len(ys), f"{c}",
            horizontalalignment = "center", verticalalignment = "center")
    plt.savefig(f"{outfile}.pca.svg")
    plt.close()

def extract_kmonomers(pkl, k):

    def ren(s):
        s = re.sub("horID_", "", s)
        s = re.sub(".mon_", "-", s)
        return s

    reads = pickle.load(open(pkl, "rb"))

    c = Counter()
    for r in reads:
        is_ori_plus = True if len([ m for m in r.mons if m.ori == '+' ]) / len(r.mons) >= 0.5 else False
        if is_ori_plus:
            for i in range(len(r.mons)-k+1):
                if all([r.mons[j+1].begin - r.mons[j].end < 100 for j in range(i, i+k-1)]):
                    c.update( [ "\t".join([ ren(m.monomer.name) for m in r.mons[i:i+k] ]) ] )
        else:
            rev_mons = list(reversed(r.mons))
            for i in range(len(rev_mons)-k+1):
                if all([rev_mons[j+1].end - rev_mons[j].begin < 100 for j in range(i, i+k-1)]):
                    c.update( [ "\t".join([ ren(m.monomer.name) for m in rev_mons[i:i+k] ]) ] )

    for i, n in c.items():
        print(f"{n}\t{i}")

def draw_HOR(hers_list, path, cmap, with_unit = None):

    # load color map
    cmap = { k:v for k,v in [ l.strip("\n").split('\t') for l in open(cmap, "r").readlines() if (len(l) > 1) and (l[0] != "#") ] }

    import hashlib
    def m2c(n):
        if n in cmap:
            return cmap[n]
        elif n[:3] == "M=M":
            return "#aaaaaa"
        else:
            return f"#{hashlib.md5(n.encode()).hexdigest()[:6]}"

    if with_unit:
        sample_hers = [ (s, [ h for h in sh if with_unit in [ t for i, s, t in h.hors ] ][:80]) for s, sh in hers_list ]

    def leftmost_target(her): 
        if her.ori == "+":
            return min([ her.mons[int(_mix)].begin
                for _mix, s, elm in her.hors if elm == with_unit ])
        else:
            return min([ -1 * her.mons[int(_mix)].end
                for _mix, s, elm in her.hors if elm == with_unit ])

    def rightmost(her): 
        if her.ori == '+':
            return max([ her.mons[int(_mix) + s -1].end for _mix, s, e in her.hors ])
        else:
            return max([ -1 * her.mons[int(_mix) + s -1].begin for _mix, s, e in her.hors ])

    xs, ys = 0.05, 25 # scale
    max_lm = max([0] + [ leftmost_target(her) for s, hers in sample_hers for her in hers ])

    read_seen = 0
    next_reads = [ her for s, hers in sample_hers for her in hers ][read_seen:read_seen + 500]
    max_x = xs * max([0] + [ rightmost(her) + max_lm - leftmost_target(her) for her in next_reads ])
    max_y = ys * len(next_reads)

    dwg = svgwrite.drawing.Drawing(
            path + f"-0.svg",
            size = (max_x + 10, max_y + 40 * ys),
            style="font-size:10;font-family:Arial;")

    chunk = 0
    y_acc = 0 # n of rows

    for sample, hers in sample_hers:

        y_acc += 1
        dwg.add(dwg.text(f"{sample}", insert = (max_lm * xs, (y_acc + 0.6) * ys), style="font-size:15;font-family:Arial;font-weight:bold;"))

        for i, her in enumerate(hers):
            y_acc += 1
            read_seen += 1
            read = dwg.g(id=her.name) 
            lmt = leftmost_target(her)

            for hi, (_mix, _size, elm) in enumerate(her.hors):
                mix, size = int(_mix), int(_size)
                if her.ori == '+':
                    b, e = her.mons[mix].begin, her.mons[mix + size -1].end
                else:
                    b, e = -1 * her.mons[mix].end, -1 * her.mons[mix + size -1].begin
                b += max_lm - lmt
                e += max_lm - lmt
                read.add(dwg.rect(
                    insert=((b + 5) * xs, (y_acc + 0.2) * ys), size=((e-b-10) * xs, 0.7 * ys),
                    fill = m2c(elm), fill_opacity = 0.5))
                #if elm[0:3] == "M=M":
                #    read.add(dwg.text(".", insert = ((b+10) * xs, (y_acc + 0.4) * ys)))
                if elm[0:2] == "M=":
                    read.add(dwg.text(".", insert = ((b + 10) * xs, (y_acc + 0.4) * ys)))
                else:
                    read.add(dwg.text(f"{elm}", insert = ((b + 10) * xs, (y_acc + 0.4) * ys)))

            dwg.add(read)
            print("", end = ".")

            if y_acc > 500:
                dwg.save()
                print(f"done chunk {chunk}")
                chunk += 1
                y_acc = 0
                # new chunk
                next_reads = [ her for s, hers in sample_hers for her in hers ][read_seen:read_seen + 500]
                max_x = xs * max([0] + [ rightmost(her) + max_lm - leftmost_target(her)
                    for her in next_reads ])
                max_y = ys * len(next_reads)
                dwg = svgwrite.drawing.Drawing(path + f"-{chunk}.svg",
                        size = (max_x + 10, max_y + 40 * ys),
                        style="font-size:10;font-family:Arial;")

    dwg.save()


# NOTE: for X only; TODO move to the impl of print-hor
def cenpb_vars(rd, h, summary = False):
    """ takes a unit (specified by encoded read obj and index),
    and returns 9x7 (0,3,4,5,7,10,11) matrix which expresses statuses of variants 
    on CENP-B box recognition motif. The default pattern would be 
    [
    [0,0,0,1,0,0,0,0,0], // 0
    [0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0], // 4*
    [0,0,0,1,0,0,0,0,0], // 5
    [0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,0,0],
    [0,0,0,0,0,0,0,1,1]  // 11*
    ].
    If `summary` is set, just returns tuple like (59, 4),
    that is, the number of intact sites, motifs """

    v = { (k, t.pos) for k in [0, 3, 4, 5, 7, 10, 11] for t in rd.mons[h+k].monomer.snvs }
    os = { 0:0, 3:0, 4:-1, 5:0, 7:0, 10:0, 11:-1 } # offset
    res = np.array([
        1 if (k, p + os[k]) in v else 0
        for k in [0, 3, 4, 5, 7, 10, 11]
        for p in [129, 130, 131, 132, 137, 140, 141, 142, 143]]).reshape(7, 9)

    # variants present in ref, and their back mutations
    vb = { (k, t.pos, t.base) for k in [0, 5, 11] for t in rd.mons[h+k].monomer.snvs }
    res[0,3] = 0 if (0, 132, 'G') in vb else 1
    res[3,3] = 0 if (5, 132, 'G') in vb else 1
    res[6,7] = 0 if (11, 141, 'G') in vb else 1
    res[6,8] = 0 if (11, 142, 'G') in vb else 1

    if summary:
        intact_motifs = np.all(res == 0, axis = 1).sum()
        intact_sites = 63 - np.sum(res)
        return (intact_sites, intact_motifs)
    else:
        return res

# NOTE: CENP-B option is currently only for X
# NOTE: returns candidates tracts of novel HOR variant
def print_HOR_read(r, show_cenpbbxx = False):

    coords = [ 0 ] * len(r.hors)
    hv_str = ["."] * len(r.hors)

    for n, (_idx, _size, elem) in enumerate(r.hors):

        idx, size = int(_idx), int(_size)
        if r.ori == '+':
            b, e = r.mons[idx].begin, r.mons[idx + size - 1].end
            gap = 0 if idx == 0 else r.mons[idx].begin - r.mons[idx-1].end # gap before me.
        else:
            b, e = r.mons[idx].end, r.mons[idx + size - 1].begin
            gap = 0 if idx == 0 else -(r.mons[idx].end - r.mons[idx-1].begin)

        coords[n] = (b, e, gap)

        if elem[:5] == "M=HOR":
            hv_str[n] = "m" if gap < 100 else "M"
        if elem[:3] == "M=M":
            hv_str[n] = "m" if gap < 100 else "M"
        elif elem[:3] == "Rev":
            hv_str[n] = "r" if gap < 100 else "R"
        elif elem[:2] == "M=":
            hv_str[n] = "a" if gap < 100 else "A"
        else:
            hv_str[n] = "h" if gap < 100 else "H"

    # search variant HOR candidates
    elems = [ e for i, s, e in r.hors ]
    hv_new_idx = functools.reduce(
            lambda x, y: x | y,
            [ set(range(mt.start(1)+1, mt.end(1)-1)) for mt
                in re.finditer(r"(?=([hH]m+h))", "".join(hv_str)) ],
            set())
    hv_new_hash = { mt.start(1) + 1 : f"{hash(tuple(elems[ mt.start(1) + 1: mt.end(1) - 1 ])):x}"[-8:] for mt
            in re.finditer(r"(?=([hH]m+h))", "".join(hv_str)) }

    for n, (_idx, _size, elem) in enumerate(r.hors):
        b, e, gap = coords[n]
        idx, size = int(_idx), int(_size)
        nvars = sum([ len(m.monomer.snvs) for m in r.mons[idx:idx+size] ])
        if n in hv_new_hash: 
            print( f"{r.name}\t{b}\t{e}\t{idx}\t{size}\t{elem}\t" +\
                   f"{gap}\t{nvars}\t{100.0*nvars/abs(e-b):.2f}\t" +\
                   "new=" + hv_new_hash[n])
        else:
            print( f"{r.name}\t{b}\t{e}\t{idx}\t{size}\t{elem}\t" +\
                   f"{gap}\t{nvars}\t{100.0*nvars/abs(e-b):.2f}\t" +\
                   ("new" if n in hv_new_idx else "."))

    print("")

    # 5-mer
    kmer = [ tuple(elems[ mt.start(1) : mt.end(1) ]) for mt in re.finditer(r"(?=([hH]hhhh))", "".join(hv_str)) ]
    cand = [ tuple(elems[ mt.start(1) + 1 : mt.end(1) - 1 ]) for mt in re.finditer(r"(?=([hH]m+h))", "".join(hv_str)) ]

    return kmer, cand

def print_HOR(pkl, show_cenpbbxx = False):
    """ taking a pickled HOR encoded reads, outputs HOR structure of the reads. """

    hors = pickle.load(open(pkl, "rb"))
    print("#readname\tbegin\tend\tidx\tsize\telem\tgap\tvars")
    candidates = Counter()
    hor_kmer = Counter()

    for r in sorted(hors, key=lambda x: -len(x.mons)):
        kmer, cand = print_HOR_read(r, show_cenpbbxx)
        candidates.update(cand)
        hor_kmer.update(kmer)

    print("##########")
    print("\n".join([ "#" + f"{hash(p):x}"[-8:] + f"\t{c}\t{len(p)}\t" + "\t".join(p) for p, c in candidates.most_common()  if c >= 2 ]))

    print("##########")
    print("\n".join([ "%" + f"{hash(p):x}"[-8:] + f"\t{c}\t" + "\t".join(p) for p, c in hor_kmer.most_common()  if c >= 2 and c < 100 ]))


def count_HOR_kmer(pkl):
    """ taking a pickled HOR encoded reads, outputs HOR k-mer counter. """
    counter = Counter()
    hors = pickle.load(open(pkl, "rb"))
    for r in sorted(hors, key=lambda x: -len(x.mons)):
        fourmers = list(windowed([elem for _idx, _size, elem in r.hors], 5, fillvalue="."))
        fourmers = [ fm for fm in fourmers if not any([ ("M=" in f) for f in list(fm) ]) ]
        counter.update(fourmers)
    for t, c in counter.most_common(n=500):
        print(f"{c}\t" + "\t".join([ f"{e}" for e in t ]))

def HOR_encoding(pkl, path_merged, path_patterns):

    # NOTE: copied from kmer_analysis.py
    def load_dict(path, sep = "\t"):
        return { k:v for k,v in [l.strip("\n").split(sep)
            for l in open(path, "r").readlines() if (len(l) > 1) and (l[0] != "#")] }

    def load_patterns(path):
        return { tuple([ a for a in lsp[1:] if a]) : lsp[0] 
                for lsp in [l.strip("\n").split("\t") for l in open(path, "r").readlines() if (len(l) > 1) and (l[0] != "#")] }

    merged = load_dict(path_merged)
    patterns = load_patterns(path_patterns)
    pat_str = { p : "#".join(p) for p, c in patterns.items() }
    pat_size = { len(p) for p in patterns.keys() }

    def ren(s, rc = False): # rename
        if s in merged:
            s = merged[s]
        return s + "_RC" if rc else s

    def hor_encode_read(er):

        is_ori_plus = True if len([ m for m in er.mons if m.ori == '+' ]) / len(er.mons) >= 0.5 else False
        if is_ori_plus:
            mons = er.mons
        # list large gap position (gap after me)
            gaps = [ mons[i+1].begin - mons[i].end for i in range(len(mons)-1) ] + [0]
        # renamed monomers in the read
            ren_mons = [ ren(m.monomer.name, True if m.ori == "-" else False) for m in mons ]
        else:
            mons = list(reversed(er.mons))
            gaps = [ mons[i].begin - mons[i+1].end for i in range(len(mons)-1) ] + [0]
            ren_mons = [ ren(m.monomer.name, True if m.ori == "+" else False) for m in mons ]

        # find patterns : NOTE: this can be a bit faster
        found = []
        for ps in pat_size:
            for i in range(len(mons) - ps + 1):
                if all([ gaps[j] < 100 for j in range(i, i + ps - 1) ]):
                    rd_str = "#".join(ren_mons[i:i+ps])
                    found += [ (i, i+ps, c) for p, c in patterns.items() if len(p) == ps and rd_str == pat_str[p] ]

        # find best layout of patterns
        s = [ 0 for i in range(len(mons) + 1) ]
        b = [ 0 for i in range(len(mons) + 1) ]
        t = [ "" for i in range(len(mons)) ]

        for i in range(len(mons)):
            # TODO: this is naive impl
            found = [ f for f in found if not f[1] <= i  ]
            s[i], b[i], t[i] = max( [ (s[f[0]-1] + 2 * (f[1] - f[0]) - 1, f[0]-1, f) for f in found if f[1]-1 == i ] + [(s[i-1] - 1, i-1, f"M={ren_mons[i]}")] )

        # report.
        result = []
        last = len(mons)-1
        while last >= 0:
            if t[last][0] == "M":
                g = gaps[last-1]
                result += [ (last, 1, t[last], g) ]
            else:
                g = gaps[t[last][0]-1]
                size = t[last][1]-t[last][0]
                result += [ (last-size+1, size, t[last][2], g) ] # idx, size, symbol, gap
            last = b[last]

        return HOR_Read(name = er.name, mons = mons, hors = [ h[:-1] for h in sorted(result) ],
                length = er.length, ori = '+' if is_ori_plus else '-')

    reads = pickle.load(open(pkl, "rb"))
    return [ hor_encode_read(r) for r in reads ]

# TODO: just copied from get_mismatches.py
if __name__ == '__main__':

    # TODO: write menu
    import argparse
    parser = argparse.ArgumentParser(description='Breakup encoded read based on the set of assigned monomers.')
    parser.add_argument('action', metavar='action', type=str, help='action to perform: segregate, print(temporary), kmer, encode-hor, prinr-hor, show...')

    # for segregate
    parser.add_argument('--pickles', dest='pickles', help='fofn of pickled encoded reads')
    parser.add_argument('--outdir', dest='outdir', help='output directory to which clustered reads are to be put')

    parser.add_argument('--reads', dest='reads', help='pickled encoded reads')
    parser.add_argument('--hor-reads', dest='hors', help='pickled encoded reads with hor encoding')
    parser.add_argument('-k', dest='k', help='k for k-mer analysis')
    parser.add_argument('--ref-fa', dest='reffile', help='path to reference .fa file')
    parser.add_argument('--bam', dest='bamfile', help='path to BAM format file (where short reads map into monomer ref)')
    parser.add_argument('--merged', dest='merged', help='path to tab-delimited table of monomers to be merged for pattern matching')
    parser.add_argument('--patterns', dest='patterns', help='path to tab-delimited patterns file')
    parser.add_argument('--out', dest='outfile', help='path into which hor encoded reads will be written')

    parser.add_argument('--cenpbbxx', dest='cenpbbxx', action='store_const', const=True, help='add variants data on CENP-box motif (X only)')

    # for draw-hor
    parser.add_argument('--cmap', dest='cmap', help='color map used for HOR drawing')
    parser.add_argument('--with-unit', dest='with_unit', help='reads containing the specified unit will be drawn')

    #parser.add_argument('', dest='', help='')
    args = parser.parse_args()

    #bamfile = "data/SRR1997411/alignments/SRR1997411.join.aligned.sort.bam"
    #reffile = "data/monomers/MigaKH.HigherOrderRptMon.fa"


    if args.action == "segregate":
        assert args.pickles, "pickle fofn is not specified"
        assert args.outdir, "output dir is not specified"
        segregate_reads(args.pickles, args.outdir)

    elif args.action == "draw-cluster":
        assert args.pickles, "pickle fofn is not specified"
        assert args.outfile, "specify output svg filename"
        draw_cluster(args.pickles, args.outfile, precomputed = False)
        #draw_cluster(args.pickles, args.outfile, precomputed = True)

    elif args.action == "print":
        #NOTE: this is temporary
        assert args.reads, "encoded reads pickle is not specified"
        #assert args.reffile, "ref file is missing"
        print_reads(args.reads)
        
    elif args.action == "type-stats":
        assert args.reads, "encoded reads pickle is not specified"
        type_stats(args.reads)

    elif args.action == "kmer":
        assert args.reads, "encoded reads pickle is not specified"
        assert args.k, "k-mer size is not specified"
        extract_kmonomers(args.reads, int(args.k))

    elif args.action == "encode-hor":
        assert args.reads, "encoded reads pickle is not specified"
        assert args.merged, "why not merging monomers?"
        assert args.patterns, "patterns file is not specified"
        assert args.outfile, "specify output path for HOR encoded read"
        hor_reads = HOR_encoding(args.reads, args.merged, args.patterns)
        pickle.dump(hor_reads, open(args.outfile, "wb"))

    elif args.action == "print-hor":
        assert args.hors, "specify path to HOR encoded read"
        print_HOR(args.hors, show_cenpbbxx = args.cenpbbxx)
        # print_HOR(args.hors)

    elif args.action == "count-hor":
        assert args.hors, "specify path to HOR encoded read"
        count_HOR_kmer(args.hors)

    elif args.action == "draw-hor":
        # write the structure of the reads to svg
        assert args.hors, "specify path to HOR encoded read"
        assert args.outdir, "output path is not specified"
        assert args.cmap, "specify path to cmap (even if it's empty)"

        hers_path = [ l.strip("\n").split('\t') for l in open(args.hors, "r").readlines() if (len(l) > 1) and (l[0] != "#") ]
        hers_list = [ (s, sorted(pickle.load(open(h, "rb")), key = lambda x: -len(x.mons))) for s, h in hers_path ]

        # hor_reads = pickle.load(open(args.hors, "rb"))
        if args.with_unit:
            ulist = [ l.strip("\n").split('\t')[0] for l in open(args.with_unit, "r").readlines() if len(l) ]
            print(ulist)
            for u in ulist:
                draw_HOR(hers_list, args.outdir + f"v{u[1:]}", args.cmap, with_unit = u)

        #else:
        #    draw_HOR(sorted(hor_reads, key = lambda x: -len(x.mons)), args.outdir, args.cmap)

    elif args.action == "NOP":
        assert args.bamfile, "bam file is missing"
        #absolute_frequency_distribution(args.bamfile)

