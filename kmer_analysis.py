# calculate statistics on k-mers within encoded reads in each read cluster.
from collections import Counter

def load_dict(path, sep = "\t"):
    return { k:v for k,v in [l.strip("\n").split(sep) for l in open(path, "r").readlines()] }

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

# TODO: parametrize
d = load_dict("./dict.dat")
# TODO: parametrize
ks = load_kmer_stats("./kmer_raw_analysis/K30C7.tsv")

c = consolidate_kmer_stats(ks, d)
emit_kmer_stats(c)
