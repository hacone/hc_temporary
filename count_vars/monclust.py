## Calculate sequence-based grouping for 1868 monomers; these basic data will be used to subset monomers.
import numpy as np
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

def cluster(mon, distmat, monfreq = defaultdict(1)):
    skips = set()
    for i in range(len(mon.items())) if not i in skips:



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
    else:
        print(f"unknown action. {args.action}")
