# Calculate monomer composition for each read in single pickle file
# Output as tab-delimited table

from collections import Counter 
import numpy as np
import pickle
from EncodedRead import *

if __name__ == '__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Calculate monomer composition for each read in single pickle file')
    parser.add_argument('--reads', dest='reads', help='monomer encoded reads (single pickle file)')
    parser.add_argument('--ref', dest='ref', help='list of monomers to be counted (such as .fai file)')
    args = parser.parse_args()

    assert args.reads, "Need --reads"
    assert args.ref, "Need --ref"

    reads = pickle.load(open(args.reads, "rb"))
    l = np.loadtxt(args.ref, dtype = "U20", delimiter = "\t", usecols = (0))
    mon_to_id = { n : i for i, n in enumerate(l) }

    nreads, nmons = len(reads), len(mon_to_id)
    occ = np.zeros(nreads * nmons).reshape(nreads, nmons)

    # header.
    print("#readname            \t" + "\t".join(l))
    # contents.
    for i, r in enumerate(reads):
        for m in r.mons:
            occ[i, mon_to_id[m.monomer.name]] += 1
        print(f"{r.name:<20}\t" \
                + "\t".join([ f"{int(occ[i, mi]):d}" for mi in range(nmons) ]))
