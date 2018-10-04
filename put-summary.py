import pickle
from Layout import *

lid = -1

def sum_los(los, desc=""):
    def sum_lo(lo):
        global lid
        lid += 1
        return "\n".join(f"{lid:>3} {i:>4} {ai:>2} {x:>4}" for (i,ai),x in sorted(lo.reads, key=lambda a: a[1]))
    return "\n\n".join(f"> Layout {desc}: {k}\n" + sum_lo(los[k]) for k in range(len(los)) if len(los[k].reads) > 4)

allo = pickle.load(open("./layouts.pickle", "rb"))
print("Layouts for all reads:\n\n" + sum_los(allo, "All"), file=open("allo.dat", "w"))
print("1")

nslo = pickle.load(open("./non-slip-layouts.pickle", "rb"))
print("Layouts for non-slip reads:\n\n" + sum_los(nslo, "Non-slip"), file=open("nslo.dat", "w"))
print("2")

for i, slo_pickle in enumerate(list(os.scandir("./slip-layouts/"))):
    print(slo_pickle)
    slo = pickle.load(open(slo_pickle, "rb"))
    print(sum_los(slo, f"Slip-{i}"), file=open("slo.dat", "a"))
