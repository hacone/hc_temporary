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
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

l = np.loadtxt("cluster.s14.SRR3189741.fa.fai", dtype = "U20", delimiter = "\t", usecols = (0))
mon_to_id = { i:n for n, i in enumerate(l) }

blast_assign_dir = "./pacbio/blast/blast_assignment/"
pickles = [ blast_assign_dir + path for path in os.listdir(blast_assign_dir) if path.endswith('.pickle') ]
print(len(pickles))

#picklefile = "./pacbio/blast/blast_assignment/m160427_051149_00116_c100976532550000001823226708101601_s1_p0.pickle"

read_list, occ_list = [], []

for picklefile in pickles:
    print("handling " + picklefile)
    with open(picklefile, "rb") as f:
        p = pickle.load(f)
        for r in p:
            occ = np.zeros(len(mon_to_id))
            n = len(r.mons)
            for m in r.mons:
                #occ[mon_to_id[m.monomer.name]] += 1
                occ[mon_to_id[m.monomer.name]] += 1/n
            occ_list += list(occ)
            read_list += [r.name]

all_data = np.array(occ_list).reshape(len(read_list), len(mon_to_id))

#print(all_data)
print(all_data.shape)

#reduced = PCA(n_components=2).fit_transform(all_data)
#plt.scatter(reduced[:, 0], reduced[:, 1])
#plt.savefig("PCA_encoded_reads.svg")

reduced = TSNE(n_components=2, random_state=0).fit_transform(all_data)
plt.scatter(reduced[:, 0], reduced[:, 1])
plt.savefig("TSNE_encoded_reads.svg")

# TODO t-SNE
