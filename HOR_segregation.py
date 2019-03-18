# segregate/encode reads by HOR patterns contianed in it.

from collections import Counter 
from EncodedRead import *
import numpy as np
import os
import pickle
import re
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import svgwrite

# NOTE: here's data structure for HOR encoded reads (cf. EncodedRead.py for other lower-level representation)
# hors := (index, size, symbol)
HOR_Read = namedtuple("HOR_Read", ("name", "mons", "hors", "length", "ori"))


def load_encoded_reads(pickles, n_max_reads = None):
    reads = []
    for picklefile in [ l.strip() for l in open(pickles, "r").readlines() ]:
        rds = pickle.load(open(picklefile, "rb"))
        reads += [ r for r in rds if len(r.mons) > 29 ]
        print(f"{len(reads)} reads found... " + "loaded " + picklefile, flush = True)
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

# NOTE: currently not visible from the menu
def draw_cluster(pickles, outfile, precomputed = False):
    """ t-SNE / PCA embedding of reads based on monomer sharing to visualize the clusters. """

    reads = load_encoded_reads(pickles)
    # NOTE: use 20k reads with most monomers 
    occ = monomers_in_reads(sorted(reads,
        key=lambda x: -len(x.mons))[:20000])

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
    cls = cluster_reads(occ)

    # NOTE: charm
    import matplotlib
    matplotlib.use('Agg')
    # mycm = plt.get_cmap("tab20b") + plt.get_cmap("tab20c") # TODO: how can i concatenate?
    plt.scatter(
            tsne_red[:, 0], tsne_red[:, 1], c=cls,
            s=6, alpha=0.5, edgecolors="none", cmap=plt.get_cmap("tab20b"))
    plt.savefig(f"{outfile}.tsne.svg")
    plt.scatter(
            pca_red[:, 0], pca_red[:, 1], c=cls,
            s=6, alpha=0.5, edgecolors="none", cmap=plt.get_cmap("tab20b"))
    plt.savefig(f"{outfile}.pca.svg")

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

# TODO; write up
def show_HOR(hors):

    def show_svg_HOR(dwg, ers, hors):
        """
        write out encoded reads with HOR detected to SVG.
        """
        import hashlib
        m2c = lambda n: f"#{hashlib.md5(n.encode()).hexdigest()[:6]}"
        b2c = dict(A = "#F8766D", C = "#7CAE00", G = "#00BFC4", T = "#C77CFF")

        def add_read(dwg, er, hor, offsets, thickness = 8):
            """
            (from show-reads.py)
            er: an encoded read
            hor: detected hor units ( [begin, end, type] for now )
            offsets: tuple as (xoff, yoff)
            """

            read = dwg.g(id=er.name)
            # TODO: clean up here
            gaps = [ er.mons[i+1].begin - er.mons[i].end for i in range(len(er.mons)-1) ] + [0] # gap after i
            gaps_inserted = []

            for i, m in enumerate(er.mons):
                #read.add(dwg.rect(insert=(b, y_offset), size=(e-b, thickness), fill=f"#{name_to_color(name)}"))
                mx, my = offsets[0] + 100 * (i+len(gaps_inserted)), offsets[1]
                read.add(dwg.rect(
                    insert=(mx, my), size=(100*0.9, thickness),
                    fill = m2c(ren(m.monomer.name)), fill_opacity = 0.5)) # backbone
                read.add(dwg.text(f"{ren(m.monomer.name)}", insert = (mx+10, my-4)))

                for snv in m.monomer.snvs: # SNVs
                    #read_shape.add(dwg.line(
                    read.add(dwg.circle(center = (mx + 90 * snv.pos / 180 , my), r = 2,  fill = b2c[snv.base]))

                if gaps[i] > 50:
                    read.add(dwg.text(f"{gaps[i]}", insert = (mx+90, my)))
                    gaps_inserted += [ i for j in range(round((gaps[i]-60)/180)) ]

            for b, e, t in hor: # HOR structure
                gi = len([i for i in gaps_inserted if i<b])
                
                mx0, mx1, my = offsets[0] + 100 * (b+gi), offsets[0] + 100 * (e+gi), offsets[1]
                read.add(dwg.rect(insert = (mx0+10, my+8), size = (mx1-mx0+80, 3), fill = m2c(f"{t}")))
                read.add(dwg.text(f"{t},{b}-{e}:{gi}", insert = (mx0, my+20)))

            return dwg.add(read)

        # TODO: clean up
        x, y = 10, 10
        max_off = max( [ min( [a[0] for a in hors[er.name] if a[2] == "D39"] ) for er in ers ] ) + 80
        for er in ers:
            #add_read(dwg, er, hor, (x, y))
            xoff = min ( [ a[0] for a in hors[er.name] if a[2] == "D39" ] )
            gaps = [ er.mons[i+1].begin - er.mons[i].end for i in range(len(er.mons)-1) ] + [0] # gap after i
            ng = sum([ round(gaps[i]-60)/180 for i in range(xoff) if gaps[i] > 60 ])
            print(f"xoff = {xoff}")
            add_read(dwg, er, hors[er.name], ((x + (max_off - xoff - ng) * 100), y))
            y += 45

    YFP = "D39"
    # hor_show[r.name] = [ (res[0], res[0]+res[1]-1, res[2]) for res in sorted(result) if not res[2][0] == "M" ]

    # show HOR in SVG
    dwg = svgwrite.drawing.Drawing("./hors.svg")
    show_svg_HOR(dwg, ers_show, hor_show)
    dwg.save()

def print_HOR_read(r):
    """ exposed for later use. """
    for _idx, _size, elem in r.hors:

        idx, size = int(_idx), int(_size)

        if r.ori == '+':
            b, e = r.mons[idx].begin, r.mons[idx + size - 1].end
            gap = 0 if idx == 0 else r.mons[idx].begin - r.mons[idx-1].end # gap before me.
        else:
            b, e = r.mons[idx].end, r.mons[idx + size - 1].begin
            gap = 0 if idx == 0 else -(r.mons[idx].end - r.mons[idx-1].begin)

        nvars = sum([ len(m.monomer.snvs) for m in r.mons[idx:idx+size] ])
        print( f"{r.name}\t{b}\t{e}\t{idx}\t{size}\t{elem}\t{gap}\t{nvars}")
    print("")

def print_HOR(pkl):
    """ taking a pickled HOR encoded reads, outputs HOR structure of the reads. """

    hors = pickle.load(open(pkl, "rb"))
    print("readname\tbegin\tend\tidx\tsize\telem\tgap\tvars")
    for r in sorted(hors, key=lambda x: -len(x.mons)):
        print_HOR_read(r)


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
    pat_size = { len(p) for p in patterns.keys() }

    def ren(s): # rename
        s = re.sub("horID_", "", s)
        s = re.sub(".mon_", "-", s)
        return s if not s in merged else merged[s]

    def hor_encode_read(er):

        is_ori_plus = True if len([ m for m in er.mons if m.ori == '+' ]) / len(er.mons) >= 0.5 else False
        if is_ori_plus:
            mons = er.mons
        # list large gap position (gap after me)
            gaps = [ mons[i+1].begin - mons[i].end for i in range(len(mons)-1) ] + [0]
        else:
            mons = list(reversed(er.mons))
            gaps = [ mons[i].begin - mons[i+1].end for i in range(len(mons)-1) ] + [0]
        # renamed monomers in the read
        ren_mons = [ ren(m.monomer.name) for m in mons ]

        # find patterns : NOTE: this can be a bit faster
        found = []
        for ps in pat_size:
            for p, c in patterns.items():
                if len(p) == ps:
                    pat_str = "#".join([f"{s}" for s in p])
                    for i in range(len(mons)-ps+1):
                        if all([ gaps[j] < 100 for j in range(i,i+ps-1) ]) and "#".join([f"{s}" for s in ren_mons[i:i+ps]]) == pat_str:
                            found += [(i, i+ps, c)]

        # find best layout of patterns
        s = [ 0 for i in range(len(mons) + 1) ]
        b = [ 0 for i in range(len(mons) + 1) ]
        t = [ "" for i in range(len(mons)) ]

        for i in range(len(mons)):
            # TODO: this is naive impl
            found = [ f for f in found if not f[1] <= i  ]
            s[i], b[i], t[i] = max( [ (s[f[0]-1] + 2 * (f[1] - f[0]) - 1, f[0]-1, f) for f in found if f[1]-1 == i ] + [(s[i-1] - 1, i-1, f"M={ren(mons[i].monomer.name)}")] )

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
        draw_cluster(args.pickles, args.outfile)

    elif args.action == "print":
        #NOTE: this is temporary
        assert args.reads, "encoded reads pickle is not specified"
        #assert args.reffile, "ref file is missing"
        print_reads(args.reads)

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
        print_HOR(args.hors)

    elif args.action == "show":
        # both stdout and svg; NOTE; or only to svg? cf action print
        if not args.hors:
            assert args.reads, "please specify either HOR encoded reads or monomer encoded reads"
            assert args.merged, "why not merging monomers?"
            assert args.patterns, "patterns file is not specified"
            hor_reads = HOR_encoding(args.reads, args.merged, args.patterns)
        else:
            hor_reads = pickle.load(open(args.hors, "rb"))

        show_HOR(hor_reads)

    elif args.action == "NOP":
        assert args.bamfile, "bam file is missing"
        #absolute_frequency_distribution(args.bamfile)


