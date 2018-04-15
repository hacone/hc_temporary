import pysam
import itertools

sam = pysam.AlignmentFile("./minaln_k8_80cells.sam")


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

