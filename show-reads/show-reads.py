import collections
import itertools
import math
import re

import pysam
import svgwrite
#sam = pysam.AlignmentFile("./minaln_k8_80cells.sam")

def get_read_region(aligned_segment):
    ## get read_start and read_end
    ## (aligned coordinates in an original read including hard clipped part)
    a = aligned_segment
    read_start, read_end = a.query_alignment_start, a.query_alignment_end
    # correction for hard clipping
    if (a.cigartuples[0][0] == 5):
        read_start += a.cigartuples[0][1]
        read_end   += a.cigartuples[0][1]
    return (read_start, read_end)

def hash_color(data):
    ## random colors(3 elem tuple of 0-255) based on hash.
    ## data must be understood as binary
    import hashlib
    hs_hex = hashlib.md5(data).hexdigest()[:6]
    return (int(hs_hex[0:2], 16), int(hs_hex[2:4], 16), int(hs_hex[4:6], 16))

def get_projection():
    d = {}
    with open("./d12.cluster.def") as f:
        for line in f:
#            print( line.split() )
            s = line.split()
            #d[s[0]] = s[1]
            d[s[0]] = s[0] # no projection
    return d


# may this one take a dict
# def name_to_color(name):
def name_to_color(name, color_dict):
    # for raw lab as rgb
    den = 100
    # for rgb
    # den = 256

    # TODO: adhoc failsafe for the name problem
    name = re.sub("[^.]*$", "", name)[:-1]

    r = round(float(color_dict[name][0]) * 255 / den)
    g = round(float(color_dict[name][1]) * 255 / den)
    b = round(float(color_dict[name][2]) * 255 / den)
    return "{:02X}{:02X}{:02X}".format(r, g, b)


def add_read(dwg, name, monomers, y_offset, thickness = 10):
    """
    monomers: [(name, begin, end)]
    """
    read = dwg.g(id=name)
    for n, b, e in monomers:
        read.add(dwg.rect(insert=(b, y_offset), size=(e-b, thickness), fill=f"#{name_to_color(name)}"))
    return dwg.add(read)

## TODO: obviously temporal implementation
def pickle_to_svg(dwg):

    import pickle
    import itertools
    f = open("/glusterfs/yoshimura/monomer_list.pickle", "rb")
    me = pickle.load(f)
    print(f"got {len(me)} reads")

    y_offset = 15
    for r, mons in itertools.islice(me.items(), 200):
        dwg.add(dwg.text(r, insert=(0, y_offset - 12), fill='black', font_size = "10"))
        add_read(dwg, r, [(n, int(b), int(b)+171) for (b, e, n) in mons], y_offset)
        # [(n, b, e) for (b, e, n) in mons]
        y_offset += 30


def show_svg(sam, dwg, offset = 0, color_dict = collections.defaultdict(lambda: "000000"), label = ""):
    """
    visualize monomer-encoded long reads in svg
    input sam file must be sorted by long reads
    add one reads into dwg; 
    """
    import hashlib

    # d = get_projection()

    y = offset
    
    pit = itertools.groupby(sam.fetch(until_eof=True), key=lambda x: x.query_name)

    for readname, alns in itertools.islice(pit, 11):
    #for readname, alns in itertools.islice(gen_aln_for_a_read(sam), 11):

        print(f"{readname}")
        alns = list(itertools.filterfalse(lambda x: x.is_secondary, alns))
        aln_regs = [get_read_region(a) for a in alns]
        alns_n_regs = sorted(zip(alns, aln_regs), key=lambda x: x[1][0])

        #read_shape = dwg.defs.add(dwg.g(id=readname))
        read_shape = dwg.add(dwg.g(id=readname))

        for a, reg in alns_n_regs:
            read_shape.add(dwg.text(label + " : " + readname, insert=(0, y-6), fill='black', font_size = "10"))
            read_shape.add(dwg.rect(insert=(reg[0]/5,y), size=((reg[1]-reg[0])/5, 8),
                fill=f"#{name_to_color(a.reference_name, color_dict)}"))
                #fill=f"#{hashlib.md5(d[a.reference_name].encode()).hexdigest()[:6]}"))

        y += 350


def show_txt(sam):
    # for read in gen_aln_for_a_read(sam):
    for read in itertools.islice(gen_aln_for_a_read(sam), 5):
        # filter secondary and sort by alignment_start in the read
        readname, alns = read[0], read[1]
        alns = list(itertools.filterfalse(lambda x: x.is_secondary, alns))
        aln_regs = [get_read_region(a) for a in alns]
        alns_n_regs = sorted(zip(alns, aln_regs), key=lambda x: x[1][0])
        print(f"## {readname} - {len(alns_n_regs)} monomers")
        print(f"sign\trs\tre\tmid               \tms\tme\tmlen")

        for a, reg in alns_n_regs:
            hc = hash_color(a.reference_name.encode())
            signature = f"\x1b[48;5;{hc[0]}m  \x1b[48;5;{hc[1]}m  \x1b[48;5;{hc[2]}m  \x1b[0m"

            print(a.get_aligned_pairs(matches_only = True, with_seq = False)) 
            #print(a.get_aligned_pairs(matches_only = True, with_seq = True)) # need MD tag

            print(signature + "\t" +\
                  f"{reg[0]}\t{reg[1]}\t" \
                  f"{a.reference_name}\t{a.reference_start}\t{a.reference_end}\t"\
                  f"{sam.lengths[a.reference_id]}") 


if __name__ == "__main__":

    # TODO: write menu
    import argparse
    parser = argparse.ArgumentParser(description='Visualize encoded reads in SAM or [aligned encoded reads file].')
    #parser.add_argument('action', metavar='action', type=str, help='action to perform: distmat, ...')
    parser.add_argument('--sam', dest='samfile', help='path to SAM format file')
    parser.add_argument('--svg', dest='svgfile', help='path to output svg file')
    args = parser.parse_args()

    #dwg = svgwrite.drawing.Drawing("200_reads.svg")

    # if colors are specified as option, set it up
    color_dict = {} 
    color_file = "../MigaKH.HigherOrderRptMon.fa.colors.lab"
    with open(color_file, "r") as colors:
        for l in colors.readlines():
            ll = l.strip("\n").split("\t")
            ## most concise one; ignore all after the last dot
            rname = re.sub("[^.]*$", "", ll[0])[:-1]
            color_dict[rname] = (ll[1], ll[2], ll[3])

    #print(args.action)

    #if args.action == "distmat":
    if True:
        #assert args.samfile, "input sam file is not specified. aborting."
        #assert args.svgfile, "output file is not specified. aborting."

        #dwg = svgwrite.drawing.Drawing(args.svgfile)
        #alignments = pysam.AlignmentFile(args.samfile)
        #show_svg(alignments, 10, dwg)
        #dwg.save()
        
        #sam_8k = pysam.AlignmentFile("./nomd_k8_143cells_primary_hd1khd")
        #sam_8km = pysam.AlignmentFile("./md_k8_s25_m015_143cells_primary_hd1khd")
        #sam_7km = pysam.AlignmentFile("./md_k7_s25_m015_143cells_primary_hd1khd")

        dwg = svgwrite.drawing.Drawing("./testmin/s1m0.01.svg")
        s, m = 1, 0.01
        off = 10
        for k in [8,7]:
            for w in [5,4,3,2]:
                sam = pysam.AlignmentFile(f"../test_minialign/k{k}.w{w}.s1.m001.sam")
                show_svg(sam, dwg, offset = off, color_dict = color_dict, label = f"k={k} w={w} s=1 m=0.01")
                off += 25
        k, w = 6, 5
        sam = pysam.AlignmentFile(f"../test_minialign/k{k}.w{w}.s1.m001.sam")
        show_svg(sam, dwg, offset = off, color_dict = color_dict, label = f"k={k} w={w} s=1 m=0.01")
        dwg.save()


        #show_svg(sam_8km, dwg, 25)
        #show_svg(sam_7km, dwg, 40)
        
        s, m = 20, 0.2
        off = 10
        dwg = svgwrite.drawing.Drawing("./testmin/s20m0.2.svg")
        for k in [8,7]:
            for w in [5,4,3,2]:
                sam = pysam.AlignmentFile(f"../test_minialign/k{k}.w{w}.s20.m02.sam")
                show_svg(sam, dwg, offset = off, color_dict = color_dict, label = f"k={k} w={w} s=20 m=0.2")
                off += 25
        k = 6
        for w in [5,4,3]:
            sam = pysam.AlignmentFile(f"../test_minialign/k{k}.w{w}.s20.m02.sam")
            show_svg(sam, dwg, offset = off, color_dict = color_dict, label = f"k={k} w={w} s=20 m=0.2")
            off += 25
        dwg.save()

        s, m = 50, 0.3
        off = 10
        dwg = svgwrite.drawing.Drawing("./testmin/s50m0.3.svg")
        for k in [8,7,6]:
            for w in [5,4,3,2]:
                sam = pysam.AlignmentFile(f"../test_minialign/k{k}.w{w}.s50.m30.sam")
                show_svg(sam, dwg, offset = off, color_dict = color_dict, label = f"k={k} w={w} s=50 m=0.3")
                off += 25
        k, w = 5, 5
        sam = pysam.AlignmentFile(f"../test_minialign/k{k}.w{w}.s50.m30.sam")
        show_svg(sam, dwg, offset = off, color_dict = color_dict, label = f"k={k} w={w} s=50 m=0.3")
        dwg.save()

    else:
        print(f"unknown action. {args.action}")

    # show_svg(sam_8k,10, dwg)
    #pickle_to_svg(dwg)


