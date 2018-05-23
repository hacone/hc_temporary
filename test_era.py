from EncodedRead import *
import pickle

with open("./data/pac_mon_aligned/er.pickle.1", "rb") as f:
    ers = pickle.load(f)

for i in range(10):
    print(f"alignment for {i}")
    print(pairwise_encoded(ers[i], ers[i]))
