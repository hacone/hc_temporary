import pysam
import itertools

#sam = pysam.AlignmentFile("./minaln_k8_80cells.sam")

sam_8k = pysam.AlignmentFile("./nomd_k8_143cells_primary_hd1khd")
sam_8km = pysam.AlignmentFile("./md_k8_s25_m015_143cells_primary_hd1khd")
sam_7km = pysam.AlignmentFile("./md_k7_s25_m015_143cells_primary_hd1khd")

# TODO: may I use itertools.groupby ??
def gen_aln_for_a_read(sam):
    ## returns a generator which gives (readname, list_of_alignment)
    ls, last = [], ""
    for r in sam.fetch(until_eof=True):
        if (r.query_name != last):
            if (ls):
                yield (last, ls)
            ls, last = [], r.query_name
        ls.append(r)  
    if (ls):
        yield (last, ls)

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
    with open("./d10.cluster.def") as f:
        for line in f:
#            print( line.split() )
            s = line.split()
            d[s[0]] = s[1]
    return d

import svgwrite
dwg = svgwrite.drawing.Drawing()

def show_svg(sam, offset):
    """
    visualize monomer-encoded long reads in svg
    input sam file must be sorted by long reads
    """
    import svgwrite
    import hashlib

    d = get_projection()

    y = offset
    for read in itertools.islice(gen_aln_for_a_read(sam), 11):

        readname, alns = read[0], read[1]
        alns = list(itertools.filterfalse(lambda x: x.is_secondary, alns))
        aln_regs = [get_read_region(a) for a in alns]
        alns_n_regs = sorted(zip(alns, aln_regs), key=lambda x: x[1][0])

        #read_shape = dwg.defs.add(dwg.g(id=readname))
        read_shape = dwg.add(dwg.g(id=readname))

        for a, reg in alns_n_regs:
            read_shape.add(dwg.rect(insert=(reg[0],y), size=(reg[1]-reg[0], 10),
                fill=f"#{hashlib.md5(d[a.reference_name].encode()).hexdigest()[:6]}"))

        y += 55



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

    # if __name__ == "__main__":

show_svg(sam_8k,10)
show_svg(sam_8km,25)
show_svg(sam_7km,40)

dwg.save("d20.svg")
