# calculate statistics on k-mers within encoded reads in each read cluster.
from collections import Counter
import sys

def load_dict(path, sep = "\t"):
    return { k:v for k,v in [l.strip("\n").split(sep) for l in open(path, "r").readlines() if (len(l) > 1) and l[0] != "#" ] }

def load_kmer_stats(path):
    return { tuple(m[1:]) : int(m[0]) for m in [l.strip("\n").split("\t") for l in open(path, "r").readlines()] }

def consolidate_kmer_stats(kmer_stats, merged):
    """
    kmer_stats = dict(k-tuple = freq)
    merged = dict(mon = pi(mon))
    """
    consolidated = Counter()
    for k,v in kmer_stats.items():
        consolidated[ tuple([merged[m] if m in merged else m for m in k]) ] += v
    return consolidated

def emit_kmer_stats(kmer_stats):
    for k,v in kmer_stats.items():
        print(f"{v}\t" + "\t".join(k))

#print(sys.argv)
#sys.exit()

# TODO: parametrize
# d = load_dict("./dict.dat")
d = load_dict(sys.argv[1])

# TODO: parametrize
#ks = load_kmer_stats("./kmer_raw_analysis/K30C7.tsv")
ks = load_kmer_stats(sys.argv[2])

c = consolidate_kmer_stats(ks, d)

emit_kmer_stats(c)
