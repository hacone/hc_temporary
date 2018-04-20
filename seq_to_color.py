## map sequences into color space, translating pairwise sequence distances into perceptual difference of colors

## sequence distance can be measured as edit-distance between the two.
## perceptual difference of colors is a topic extensively studied, but here, simply, CIELAB is used

## TODO: organize the usage
## param: mode=

from skimage.color import lab2rgb 
import numpy  as np 
import sys, re

def openFasta(path):
    """ open fasta as simple dict """

    from Bio.SeqIO import FastaIO
    with open(path) as handle:
        return dict(FastaIO.SimpleFastaParser(handle))
        #return dict(list(FastaIO.SimpleFastaParser(handle))[:30])

# TODO: export this
def seq_to_edit_distance(seq_dict):
    """ pairwise edit distance is calculated from dict of type {seqname, seq}"""

    from Bio import pairwise2
    #dim =  len(list(seq_dict.keys()))
    dim =  len(seq_dict.keys())
    d = np.zeros((dim, dim), dtype=np.int)
    s = list(seq_dict.values())

    for i in range(dim):
        print(i)
        sys.stdout.flush()

        for j in range(dim):
            if i == j:
                d[i,j] = 0 # put 0
            elif i > j:
                d[i,j] = d[j,i] # put d[i][j] = d[j][i]
            else:
                score = pairwise2.align.globalms(s[i], s[j], 0, -1, -1, -1, score_only = True) # calc d[i][j]
                d[i,j] = -score
    return d

# TODO: export this
# by default it returns [0,255]*3 RGB
def distMat2colors(dm, lab = False):

    # tsne into 3dim
    from sklearn.manifold import TSNE
    dm_3dim = TSNE(n_components=3, metric="precomputed").fit_transform(dm)
    n = dm_3dim.shape[0]
    #print(dm_3dim)

    # embedding is transformed to fit in CIELAB
    from skimage.color import lab2rgb

    #rgb1 = lab2rgb(dm_3dim.reshape(n, 1, 3))
    #print(rgb1.reshape(n, 3))

    #b, t = dm_3dim.min(0), dm_3dim.max(0)
    b, t = dm_3dim.min(), dm_3dim.max()
    #valid_lab = (dm_3dim - b) * np.array([99.99, 254.99, 254.99]) / (t - b) - np.array([0, 127, 127]) # range of value = 100, +-127, +-127
    valid_lab = (dm_3dim - b) * np.array([100.0, 100, 100]) / (t - b) # percentage ??
    #valid_lab = (np.array([1, 2, 2]) * (dm_3dim - b) / (t - b)) - np.array([0, 1, 1]) # range of value = [0,1] for each

    if lab:
        # return CIELAB for each elem
        # return valid_lab
        print("lab[0-1]")
        print(valid_lab)
        print("3dim")
        print(dm_3dim)
        return valid_lab
    else:
        # return RGB
        rgb2 = lab2rgb(valid_lab.reshape(n, 1, 3)).reshape(n, 3)
        print("rgb2")
        print(rgb2)
        print("rgb2[0-255]")
        print(rgb2 * 256)
        return rgb2 * 256


path = "./data/monomers/d0.fa"
#path = "./data/monomers/MigaKH.HigherOrderRptMon.fa"
seq_dict = openFasta(path)
names, sequences = list(seq_dict.keys()), list(seq_dict.values())

# dm = np.load("./count_vars/monomers.dist.npy")
dm = seq_to_edit_distance(seq_dict)

cols = distMat2colors(dm)
with open(re.sub(".*/", "", f"{path}") + ".colors", "w") as f:
    for i in range(len(names)):
        f.write(f"{names[i]}\t{cols[i,0]}\t{cols[i,1]}\t{cols[i,2]}\n")

cols = distMat2colors(dm, lab=True)
with open(re.sub(".*/", "", f"{path}") + ".colors.lab", "w") as f:
    for i in range(len(names)):
        f.write(f"{names[i]}\t{cols[i,0]}\t{cols[i,1]}\t{cols[i,2]}\n")

