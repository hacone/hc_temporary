import numpy as np
import pandas as pd

from scipy.cluster.hierarchy import dendrogram, linkage
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

def tabulate(path, threshold = 5.0):

    t = pd.read_table(path, sep = "\s+", header = None,
            names = ["sample", "mon", "pos", "bas", "freq", "rfreq", "pval", "cfd"])

    majors = [ (row[0], row[1], row[2], row[3], row[5]) \
       for i, row in t.loc[t.rfreq >= threshold].iterrows() ]

    # samples = sorted(list({ e[0] for e in majors }))
    samples = ["Mende", "Esan", "Maasai",
        "CHM13CLR", "Finnish", "Toscani", "Ashkenazi",
        "Gujarati", "Dai", "HG005",
        "PuertoRican", "Peruvian",
        "CHM13Hifi", "NA12878Hifi"]
    id2sample = { i : s for i, s in enumerate(samples) }
    sample2id = { id2sample[i] : i for i in id2sample }

    major_sites = sorted(list({ (e[1], e[2], e[3]) for e in majors }))
    id2site = { i : t for i, t in enumerate((major_sites)) }
    site2id = { id2site[i] : i for i in id2site }

    freq_mat = np.full((14, len(major_sites)), threshold)
    for s, m, i, b, rf in majors:
        freq_mat[sample2id[s], site2id[(m, i, b)]] = rf

    return freq_mat, samples, major_sites

def projection(freq_mat, samples, opath):

    pca = PCA()
    pca.fit(freq_mat[:12, :])
    tm = pca.fit_transform(freq_mat[:12])
    fig, ax = plt.subplots()
    ax.scatter(tm[:, 0], tm[:, 1])
    for i in range(len(samples[:12])):
        ax.annotate(samples[i][:3], (tm[i,0], tm[i,1]))
    plt.savefig(f"{opath}.pca.png")

    red = TSNE(n_components=2, random_state=0).fit_transform(freq_mat[:12])
    fig, ax = plt.subplots()
    ax.scatter(red[:, 0], red[:, 1])
    for i in range(len(samples[:12])):
        ax.annotate(samples[i][:3], (red[i,0], red[i,1]))
    plt.savefig(f"{opath}.tsne.png")


def hclust(freq_mat, samples, opath):

    #for method in ["ward", "single", "average", "weighted", "complete"]:
    met_met = [ (method, metric) 
            for method in ["single", "complete", "average", "weighted"]
            for metric in ["euclidean", "cosine", "correlation"]] +\
            [ ("ward", "euclidean") ]

    for method, metric in met_met:
        Z = linkage(freq_mat, method = method, metric = metric)
        plt.figure(figsize=(15, 8))
        plt.title(f"{method}-{metric}")
        plt.ylabel('Samples')
        plt.xlabel(f"Distance")
        dendrogram(Z, labels = samples, orientation = "left")
        plt.savefig(f"{opath}.dendro.{method}_{metric}.png")

        print(f"{method}_{metric};")

def pretty_print(freq_mat, samples, major_sites, opath):

    id2sample = { i : s for i, s in enumerate(samples) }
    sample2id = { id2sample[i] : i for i in id2sample }

    id2site = { i : t for i, t in enumerate((major_sites)) }
    site2id = { id2site[i] : i for i in id2site }

    out = open(opath, "w")
    print(f"mon\tpos\tbas\t" +\
           "\t".join([f"{s}" for s in samples]) +\
           "\tAVG", file=out)

    for m, i, b in major_sites:
        print(f"{m+1}\t{i}\t{b}\t" + "\t".join([
            f"{freq_mat[sample2id[s], site2id[(m, i, b)]]:.2f}" for s in samples]) +\
            f"\t{np.mean(freq_mat, axis = 0)[site2id[(m,i,b)]]:.2f}", file=out)

for hor in [5, 11, 12, 16]:
    ipath = f"all.{hor}m.{hor}mW.vars"
    opath = f"{hor}m.vars.tab"
    # ("all.all.vars", "all.vars.tab")

    freq_mat, samples, major_sites = tabulate(ipath, threshold = 0.0)
    pretty_print(freq_mat, samples, major_sites, opath)

    continue

    opath += ".p10"
    projection(freq_mat, samples, opath)
    hclust(freq_mat, samples, opath)
    hclust(freq_mat[:12], samples[:12], opath + "_CLR")
    # hclust(ipath, opath)
    # pretty_print(ipath, opath)
