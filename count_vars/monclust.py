## Calculate sequence-based grouping for 1868 monomers; these basic data will be used to subset monomers.
from collections import defaultdict
import numpy as np
import pysam
import re

def openFasta(path):
    """ open fasta as simple dict """
    from Bio.SeqIO import FastaIO
    with open(path) as handle:
        return dict(FastaIO.SimpleFastaParser(handle))

def save_pairwise_edit_distance(monfile, outfile):
    """ pairwise edit distance is saved in the specified file """

    monomers_list = list(openFasta(monfile).values())
    dim =  len(monomers_list)
    d = np.zeros((dim, dim), dtype=np.int)

    from Bio import pairwise2
    for i in range(dim):
        print(f"{i}/{dim}")
        for j in range(dim):
            if i == j:
                d[i,j] = 0 # put 0
            elif i > j:
                d[i,j] = d[j,i] # put d[i][j] = d[j][i]
            else:
                score = pairwise2.align.globalms(monomers_list[i], monomers_list[j],
                        0, -1, -1, -1, score_only = True) # calc d[i][j]
                d[i,j] = -score

    np.save(outfile, d)

def save_tsne_embedding():

    from sklearn.manifold import TSNE
    d = np.load("monomers.dist.npy")
    d_embedded = TSNE(n_components=2, metric="precomputed").fit_transform(d)
    np.save("monomers.tsne.npy", d_embedded)

    
def print_close_pairs():

    monomers_dict = openFasta("./MigaKH.HigherOrderRptMon.fa")
    monomers_list = list(monomers_dict.values())
    monomers_name = list(monomers_dict.keys())

    dim =  len(monomers_list)
    d = np.load("monomers.dist.npy")

    skips = set()

    for i in range(dim):

        if i in skips:
            continue
        #print(f"{i} - {monomers_name[i]}")
        print(f"{re.sub(' .*', '', monomers_name[i])}\tcls_{i}_d12")
        similars = { j for j in range(i+1, dim) if (j not in skips) & (d[i,j] <= 12) }
        for j in similars:
            #print(f"\t{j} - d = {d[i,j]}: {monomers_name[i]} - {monomers_name[j]}")
            print(f"{re.sub(' .*', '', monomers_name[j])}\tcls_{i}_d12")
        skips |= similars

    #with open(f"d20.fa", "w") as outfa:
    #    for i in [i for i in range(dim) if (i not in skips)]:
    #        print(f">{monomers_name[i]}", file=outfa)
    #        print(f"{monomers_list[i]}", file=outfa)

    return skips

# MonomerClusters.pi[]
# MonomerClusters.edops[(from, to)]

def cluster(mon, d, sep_t = 10, monomer_freq = defaultdict(lambda: 0)):
    """
    perform clustering of monomers given, among which the distance matrix is precalculated and given as `d`.
    each cluster must be separated at least by `sep_t`. the most frequent (in short reads) monomers are chosen as the representatives.
    caution: results are written into standard output.
    """

    names = [re.sub("[^.]*$", "", n)[:-1] for n in mon.keys()]
    dim = len(names)
    skips = set()

    def find_similars(s):
        "collect elements close to a member of s"
        similars = set()
        for i in s:
            similars |= { j for j in range(dim) if d[i,j] <= sep_t }
        return similars

    def report(cls):
        "write out the formatted cluster information"
        cluster = sorted( [ dict(name=names[j], freq=monomer_freq[names[j]], idx=j) for j in cls ],
                key=lambda x: -(x["freq"]) )
        rep = cluster[0] # representative must be the most frequent one
        for c in cluster:
            print(f"{c['name']}\t{rep['name']}\t{c['freq']}\t{d[c['idx'], rep['idx']]}")

    def report_repseq(cls):
        "write out a Fasta record for representative sequence."
        cluster = sorted( [ dict(name=names[j], freq=monomer_freq[names[j]], idx=j) for j in cls ],
                key=lambda x: -(x["freq"]) )
        rep = cluster[0] # representative must be the most frequent one
        for n, s in mon.items():
            if re.sub("[^.]*$", "", n)[:-1] == rep["name"]:
                print(">" + rep["name"])
                print(s)

    # declare the format
    print("name\tcluster\tfreq\td_to_repr")

    for i in range(dim):
        if i in skips:
            pass
        else:
            curr_clus, p_end = {i}, False
            while not p_end:
                next_clus = find_similars(curr_clus)
                if len(next_clus) == len(curr_clus):
                    skips |= curr_clus
                    p_end = True
                curr_clus = next_clus
            # report current cluster for i
            #report(curr_clus)
            report_repseq(curr_clus)

def plot_embedding():

    import matplotlib.pyplot as plt
    d_embedded = np.load("monomers.tsne.npy")
    print(d_embedded.shape)

    skips = print_close_pairs()
    idx_to_plt = [ i for i in range(1868) if i not in skips ]

    plt.scatter(d_embedded[:,0], d_embedded[:,1], s = 1)
    plt.scatter(d_embedded[idx_to_plt,0], d_embedded[idx_to_plt,1], s = 1, c = "red")
    plt.savefig(f"t-SNE_d30.svg")

#plot_embedding()

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Interact with monomer databases.')
    parser.add_argument('action', metavar='action', type=str, help='action to perform: distmat, cluster, ...')
    parser.add_argument('--mons', dest='monfile', help='path to monomers.fa')
    parser.add_argument('--out', dest='outfile', help='path to output')
    parser.add_argument('--dist', dest='distmatfile', help='path to distance matrix numpy object')
    parser.add_argument('--cluster-size', dest='cluster_size', help='the distance by which clusters must be separated each other')
    parser.add_argument('--bam', dest='bamfile', help='bamfile where short reads are aligned to monomers, used in selecting representatives')
    args = parser.parse_args()

    print(args.action)

    if args.action == "distmat":
        assert args.monfile, "monomers database is not specified. aborting."
        if args.outfile:
            save_pairwise_edit_distance(args.monfile, args.outfile)
        else:
            save_pairwise_edit_distance(args.monfile)

    if args.action == "cluster":
        assert args.monfile, "monomers database is not specified. aborting."
        assert args.distmatfile, "distance matrix is not specified. aborting."
        assert args.cluster_size, "cluster size is not specified. aborting."
        mons = openFasta(args.monfile)
        dm = np.load(args.distmatfile)
        if args.bamfile:
            bam = pysam.AlignmentFile(args.bamfile)
            monomer_freq = { re.sub("[^.]*$", "", refname)[:-1] : len(list(bam.fetch(contig=refname))) for refname in bam.references }
            cluster(mons, dm, int(args.cluster_size), monomer_freq)
        else:
            cluster(mons, dm, int(args.cluster_size))

    else:
        print(f"unknown action. {args.action}")
