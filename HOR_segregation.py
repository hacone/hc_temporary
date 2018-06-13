# segregate reads by HOR patterns contianed in it.

# NOTE: charm
import matplotlib
matplotlib.use('Agg')

from collections import Counter 
from EncodedRead import *
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
import re
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

# TODO: get this from param
l = np.loadtxt("cluster.s14.SRR3189741.fa.fai", dtype = "U20", delimiter = "\t", usecols = (0))
mon_to_id = { i:n for n, i in enumerate(l) }
nmons = len(mon_to_id)


def load_encoded_reads(pickles, n_max_reads = None):
    reads = []
    for picklefile in pickles:
        reads += pickle.load(open(picklefile, "rb"))
        print(f"{len(reads)} reads found... " + "loaded " + picklefile, flush = True)
        if n_max_reads and (len(reads) > n_max_reads):
            break
    return reads

def monomers_in_reads(reads):
    """
    counts occurences of each reference monomer in each read.
    returns ndarray of shape (nreads, nmons)
    """
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

    from sklearn.cluster import KMeans
    kmeans = KMeans(n_clusters=n_clusters, random_state=0, verbose=1, n_jobs=-3).fit(occ_n)
    # pickle.dump(kmeans, open("{n_clusters}_cls.kmeans.pkl", 'wb'))

    return kmeans.predict(occ_n)


def segregate_reads(pickles_dir):
    pickles = [ pickles_dir + path for path in os.listdir(pickles_dir) if path.endswith('.pickle') ]
    #reads = load_encoded_reads(pickles, 20000)
    reads = load_encoded_reads(pickles)
    cls = cluster_reads(monomers_in_reads(reads))
    cls_reads = [ [] for i in range(40) ]

    for i, r in enumerate(reads):
        cls_reads[cls[i]] += [r]
    print("encoded read list is segregated, writing out")

    for i, c in enumerate(cls_reads):
        print(f"C-{i+1} has {len(c)} reads")
        pickle.dump(c, open(f"./encoded_read_clusters/C_{i+1}.pickle", "wb"))

def print_reads(pkl):
    reads = pickle.load(open(pkl, "rb"))
    for r in reads[:200]:
        print(f"{r.name}\t{len(r.mons)}\t{r.length}")
        last_end = 0
        for m in r.mons:
            print(f"{m.monomer.name}\t{m.begin}\t{m.end}\t{m.begin-last_end}\t{len(m.monomer.snvs)}")
            last_end = m.end

# NOTE: currently not visible from the menu
def tsne_reads(occ):
    """ t-SNE embedding of reads based on monomer sharing to visualize the clusters. """

    # reduced = TSNE(n_components=2, random_state=0, verbose=1).fit_transform(all_data)
    # np.save("20k_reads.tsne.npy", reduced)
    reduced = np.load("20k_reads.tsne.npy")
    # mycm = plt.get_cmap("tab20b") + plt.get_cmap("tab20c") # TODO: how can i concatenate?
    mycm = plt.get_cmap("tab20b")
    plt.scatter(reduced[:, 0], reduced[:, 1], c=cls, s=6, alpha=0.5, edgecolors="none", cmap=mycm)
    plt.savefig("TSNE_encoded_reads_20k_40c.svg")

def extract_kmonomers(pkl, k):

    def ren(s):
        s = re.sub("horID_", "", s)
        s = re.sub(".mon_", "-", s)
        return s

    reads = pickle.load(open(pkl, "rb"))

    c = Counter()
    for r in reads:
        for i in range(len(r.mons)-k+1):
            if all([r.mons[j+1].begin - r.mons[j].end < 100 for j in range(i, i+k-1)]):
                c.update( [ "\t".join([ ren(m.monomer.name) for m in r.mons[i:i+k] ]) ] )
    for i, n in c.items():
        print(f"{n}\t{i}")

def HOR_encoding(pkl, path_merged, path_patterns):


    # NOTE: copied from kmer_analysis.py
    def load_dict(path, sep = "\t"):
        return { k:v for k,v in [l.strip("\n").split(sep) for l in open(path, "r").readlines()] }
    merged = load_dict(path_merged)

    def ren(s): # rename
        s = re.sub("horID_", "", s)
        s = re.sub(".mon_", "-", s)
        return s if not s in merged else merged[s]

    def load_patterns(path):
        return { tuple(m[1:]) : m[0] for m in [l.strip("\n").split("\t") for l in open(path, "r").readlines()] }
    patterns = load_patterns(path_patterns)
    pat_size = { len(p) for p in patterns.keys() }

    reads = pickle.load(open(pkl, "rb"))
    for r in reads[:20]:
        # list large gap position
        gaps = [ r.mons[i+1].begin - r.mons[i].end for i in range(len(r.mons)-1) ] + [0]
        # renamed monnomers in the read
        ren_mons = [ ren(m.monomer.name) for m in r.mons ]

        # find patterns
        found = []
        for ps in pat_size:
            for i in [i for i in range(len(r.mons)-ps+1) if all([ gaps[j] < 100 for j in range(i,i+ps-1) ])]:
                found += [ (i, i+ps, c) for p,c in patterns.items()
                        if (len(p) == ps) and "#".join([f"{s}" for s in ren_mons[i:i+ps]]) == "#".join([f"{s}" for s in p]) ]
        # print(found, flush=True)

        # find best layout of patterns
        # s, b = np.zeros(len(r.mons)+1), np.zeros(len(r.mons)+1)
        s = [ 0 for i in range(len(r.mons) + 1) ]
        b = [ 0 for i in range(len(r.mons) + 1) ]
        t = [ "" for i in range(len(r.mons)) ]

        for i in range(len(r.mons)):
            # TODO: this is naive impl
            found = [ f for f in found if not f[1] <= i  ]
            s[i], b[i], t[i] = max( [ (s[f[0]-1] + 2 * (f[1] - f[0]) - 1, f[0]-1, f) for f in found if f[1]-1 == i ] + [(s[i-1] - 1, i-1, f"M={ren(r.mons[i].monomer.name)}")] )

        # report.
        result = []
        last = len(r.mons)-1
        while last >= 0:
            if t[last][0] == "M":
                g = gaps[last-1]
                result += [ (last, 1, t[last], g) ]
            else:
                g = gaps[t[last][0]-1]
                result += [ (last, t[last][1]-t[last][0], t[last][2], g) ]

            last = b[last]

        print(r.name, flush = True)
        print("idx\tsiz\telm\tgap")
        print("\n".join([ "\t".join([ f"{e}" for e in line ]) for line in sorted(result) ]), flush = True)


# This print out the composition of each cluster.
# NOTE: header
#print("Monomer\t" + "\t".join([f"{i+1}" for i in range(40)]))
#for i in range(len(l)):
#    print(l[i] + "\t" + "\t".join([ f"{mon_occ[j,i]}" for j in range(40) ]))

# TODO: just copied from get_mismatches.py
if __name__ == '__main__':

    # TODO: write menu
    import argparse
    parser = argparse.ArgumentParser(description='Breakup encoded read based on the set of assigned monomers.')
    parser.add_argument('action', metavar='action', type=str, help='action to perform: segregate, print(temporary), kmer, encode-hor, ...')
    parser.add_argument('--pickle-dir', dest='pickle_dir', help='directory where pickled encoded reads are placed')
    parser.add_argument('--reads', dest='reads', help='pickled encoded reads')
    parser.add_argument('-k', dest='k', help='k for k-mer analysis')
    parser.add_argument('--ref-fa', dest='reffile', help='path to reference .fa file')
    parser.add_argument('--bam', dest='bamfile', help='path to BAM format file (where short reads map into monomer ref)')
    parser.add_argument('--merged', dest='merged', help='path to tab-delimited table of monomers to be merged for pattern matching')
    parser.add_argument('--patterns', dest='patterns', help='path to tab-delimited patterns file')
    #parser.add_argument('', dest='', help='')
    args = parser.parse_args()

    #bamfile = "data/SRR1997411/alignments/SRR1997411.join.aligned.sort.bam"
    #reffile = "data/monomers/MigaKH.HigherOrderRptMon.fa"


    if args.action == "segregate":
        assert args.pickle_dir, "pickle dir is not specified"
        assert args.reffile, "ref file is missing"
        #segregate_reads(args.pickle_dir)
        segregate_reads("./pacbio/blast/blast_assignment/")

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
        HOR_encoding(args.reads, args.merged, args.patterns)

    elif args.action == "NOP":
        assert args.bamfile, "bam file is missing"
        #absolute_frequency_distribution(args.bamfile)


