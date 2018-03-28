### Clustering of monomers based on neighboring data in short read

import numpy as np
preds = np.load("preds.npy")
nexts = np.load("nexts.npy")

assert(preds.shape[0] == preds.shape[1] &
       preds.shape[1] == nexts.shape[0] &
       nexts.shape[0] == nexts.shape[1])

ndim = preds.shape[0]

dist = np.zeros(preds.shape)

# convert neighboring frequency data to distance 
for i in range(ndim):
    for j in range(ndim):
        if preds[i,j] + nexts[i,j] == 0:
            dist[i,j] = 2
        else:
            dist[i,j] = 1 / (preds[i,j] + nexts[i,j])

print(dist)

from sklearn.manifold import TSNE
dist_embedded = TSNE(n_components=2, metric="precomputed").fit_transform(dist)

# np.save("monomers.tsne.npy", d_embedded)
# d_embedded = np.load("monomers.tsne.npy")

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
plt.scatter(dist_embedded[:,0], dist_embedded[:,1], s = 1)
plt.savefig("Neighboring.t-SNE.svg")
