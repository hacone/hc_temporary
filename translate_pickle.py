# this translate pickle into namedtuple-supported pickle, which I assume for visualization using show-reads.py
from EncodedRead import *
import itertools
import pickle
import re

#refname = "m150214_101911_42172_c100745922550000001823159107071516_s1_p0/102606/0_18451"
#path_aln = "./0508/to_hx/m150214_101911_42172_c100745922550000001823159107071516_s1_p0_102606_0_18451"
#refname = "m150214_101911_42172_c100745922550000001823159107071516_s1_p0/129527/0_13731"
refname = "m150517_050859_42172_c100780482550000001823167008251562_s1_p0/11771/42_13976" ## score of 99
path_aln = "./0508/to_hx/" + re.sub("/", "_", refname)

path_reads = "/home/shubhakar/pacbio_encode_seq_no_embed_lists.pickle"
with open(path_reads, "rb") as f:
    reads = pickle.load(f)

path_snvdict = "/home/shubhakar/all_valid_snvs.pickle"
with open(path_snvdict, "rb") as f:
    snvdict = pickle.load(f)

def parse_mon(mon, read, dp = 0):

    if mon == "GAP":
        return None
    elif type(mon) is str:
        m = mon
    elif type(mon) is list:
        if mon[0] == "GAP":
            m = mon[1]
        else:
            m = mon[0]

    print(m)
    print(mon)
    b, e, name = m.split(":")
    b, e = int(b), int(e)
    
    snv_key = f"{b}#{e}#{name}#{read}"
    if snv_key in snvdict:
        snvs = [ SNV(pos = k, base = v) for k,v in snvdict[snv_key].items() ]
    else:
        snvs = []

    return AssignedMonomer(
        begin = b + dp, end = e + dp,
        monomer = Monomer(
        name = name,
        snvs = list(snvs)))
        #snvs = snvdict[f"{b}#{e}#{name}#{read}"] ))

def parse_read(a, refmons):
    """ return ER from a line of alignment file"""
    asp = a.split()
    qread = asp[1]
    disp = int(asp[5]) - int(asp[7])

    def refpos(i):
        if i >= len(refmons):
            return refmons[-1].begin + (171 * i-len(refmons))
        if i >= 0:
            if refmons[i]:
                return refmons[i].begin
            else:
                return -1000
        else:
            return refmons[0].begin + (171 * i)

    mons = [ parse_mon(m, qread, disp * 171) for m in reads[qread] ]

    def displace(m, i):
        if not m:
            return None
        else:
            return AssignedMonomer(begin = refpos(i), end = refpos(i) + 171, monomer = mons[i].monomer)

    #mons = filter(None, [ displace(mons[i], i) for i in range(len(mons)) ])
    mons = filter(None, [ displace(mons[i], i) for i in range(int(asp[7]), int(asp[8])+1) ])

    return EncodedRead(
        name = qread,
        mons = list(mons),
        length = int(asp[4]))

# set up ref
refmons = [ parse_mon(m, refname, 0) for m in reads[refname] ]

ref = EncodedRead(name = f"REF:{refname}",
        mons = filter(None, [ parse_mon(m, refname, 0) for m in reads[refname] ]),
        length = 0)
    
with open(path_aln, "r") as f:
    ers = [ parse_read(a, refmons) for a in list(f.readlines()) ]
    ers = [ er for er in ers if er.length > 50 ]

with open("pileup0507.pickle", "wb") as f:
    pickle.dump([ref] + ers, f)

g = """

readname -> list [ "pos:pos:name" or ["GAP", "pos:pos:name"] ]

# set up the reference read
mons = [ AssignedMonomer(
    begin = read[i].begin,
    end = read[i].end,
    monomer = Monomer(
        name = read[i].monomername,
        snvs = snvdict[f"{mon.begin}#{mon.end}#{read.name}"] )) for i in aligned_idx ]

ref = EncodedRead(name = f"REF:{read.name}", mons = mons.copied, length = 0)

# set up the aligned reads
mons = [ AssignedMonomer(
    begin = ref.mons[i].begin,
    end = ref.mons[i].end,
    monomer = Monomer(
        name = ,
        snvs = snvdict[f"{mon.begin}#{mon.end}#{read.name}"] )) for i in aligned_idx ]

er = EncodedRead(name = , mons = mons, length = 0)

for 


m150214_101911_42172_c100745922550000001823159107071516_s1_p0_110956_0_11878
m150214_101911_42172_c100745922550000001823159107071516_s1_p0_129527_0_13731
m150214_101911_42172_c100745922550000001823159107071516_s1_p0_144841_0_4437
m150214_101911_42172_c100745922550000001823159107071516_s1_p0_68457_0_10309
m150214_144322_42172_c100745922550000001823159107071517_s1_p0_100867_0_7361
m150214_144322_42172_c100745922550000001823159107071517_s1_p0_35932_0_11092
m150214_144322_42172_c100745922550000001823159107071517_s1_p0_50303_0_20970
m150214_144322_42172_c100745922550000001823159107071517_s1_p0_56059_0_12244
m150323_105907_42172_c100754442550000001823169207081530_s1_p0_106938_5275_20471
"""
