# Alignments.py
# evolved, rewritten from fragments in SNV_distribution.py

from collections import Counter
import numpy as np

from EncodedRead import *
from HOR_segregation import *

# 1 - First, Detect variants
def detect_snvs(units):
    """ return: [(monomer index in unit, position in monomer, alt. base, relative frequency)] """
    n_units = len(units)
    counter = Counter()
    for ri, i in units:
        for j in range(2, 12): # I'm skipping the first 2 monomers, where assignment is always wrong
            counter.update([ (j, s.pos, s.base) for s in hers[ri].mons[i+j].monomer.snvs ])
    return [ (k, p, b, c/n_units) for (k, p, b), c in counter.most_common(1000) if (0.05 < c/n_units) and (c/n_units < 0.8) ]

# 2 - Then, construct the bit representation matrix of shape (n_units, n_sites + 3).
def bit_vector_repr(ers, sites, range_of_units):
    """ obtain bit vector representation as np.ndarray 
        missing units (with spacing of 22~24, 33~35 mons) are represented as 0-masked units
        param: 
            ers = HOR encoded reads (load from pickle)
            sites = SNV sites as iterable of SNV (cf. defined as in EncodedRead.py)
            range_of_units = (possibly a subset of) units which should be encoded into bit-vec repr.
        out:
            [is_not_masked, rid, hid, SNV_states*]
    """
    bv, n_units = [], 0
    last_r, last_h = -1, -1

    # rid = index of read in ers; hid = index (in ers[ri].mons) of monomer at beginning of the HOR unit
    for rid, hid in range_of_units:
        # NOTE: here I filling up gaps due to 1 or 2 missing monomer assignment.
        if rid == last_r and (hid == last_h + 22 or hid == last_h + 23 or hid == last_h + 24):
            bv += [0, rid, last_h + 12] + [ 0 for i in sites ]
            n_units += 1
        elif rid == last_r and (hid == last_h + 33 or hid == last_h + 34 or hid == last_h + 35):
            bv += [0, rid, last_h + 12] + [ 0 for i in sites ]
            bv += [0, rid, last_h + 24] + [ 0 for i in sites ] # 2nd units
            n_units += 2
        bv += [1, rid, hid]
        for k, p, b, f in sites: # NOTE: should be one-liner?
            bv += [1] if any([ (sp.pos, sp.base) == (p, b) for sp in ers[rid].mons[hid+k].monomer.snvs ]) else [-1]
        last_r, last_h = rid, hid
    return np.array(bv).reshape(len(range_of_units) + n_units, 3+len(sites))

# 3 - determine which part of the matrix should be aligned each other.
def valid_regions(bitvec):
    """ get valid consecutive regions after filling up 22/23/24-mons or 33/34/35-mons gaps """
    regs = []
    b, last_mi, last_ri = 0, bitvec[0,2], bitvec[0,1]
    for i in range(bitvec.shape[0]):
        if (bitvec[i,1] != last_ri) or (bitvec[i,2] > last_mi + 12):
            regs += [(b, i)]
            b = i
        last_mi, last_ri = bitvec[i,2], bitvec[i,1]
    regs += [(b, bitvec.shape[0])]
    return regs

def align_two(bv1, bv2):
    """
        calculate possible alignments between 2 regions given as 2 windows of a bits matrix.
        return: [( rs, re, qs, qe, score, aln_len )]
    """

    def match_ratio(bv1, bv2):
        """ scoring scheme: n11 / (n10+n01+n11) """
        assert bv1.shape == bv2.shape, "wrong shapes"
        m = np.multiply(bv1, bv2) # 1: match, 0: masked, -1: mismatch
        mt = np.multiply((m + 1), m) / 2 # 1: match, 0: masked/mismatch
        n_m = np.sum(np.multiply(m, m) - mt) # n01 or n10
        n11 = np.sum(np.multiply(mt, (bv1 + 1) / 2)) # n11
        return int(1000 * n11 / (1 + n11 + n_m)) # per-mille !

    # length of each read
    lr, ls = bv1.shape[0], bv2.shape[0]
    # minimum units to be overlaped
    min_ovlp = 2

    # configurations in dovetail.
    confs = [ (0, i, ls-i, ls) for i in range(min_ovlp, min([ls, lr])) ]
    # contained/containing
    confs += [ (i, i+ls, 0, ls) if lr > ls else (0, lr, ls-i-lr, ls-i) for i in range(0, abs(lr-ls) + 1) ]
    # dovetail.
    confs += [ (lr-i, lr, 0, i) for i in range(min([ls, lr])-1, min_ovlp-1, -1) ]
    
    aln_results = [ (rs, re, qs, qe,
       match_ratio(bv1[rs:re,3:], bv2[qs:qe,3:]),
       sum([ 0 if not (bv1[rs+i,0] * bv2[qs+i,0]) else 1 for i in range(re-rs) ])
      ) for rs, re, qs, qe in confs ]

    # filter for ones with 50% identity and sufficient ovlp
    aln_results = list(filter(lambda x: x[4] > 500 and x[5] >= min_ovlp, aln_results))

    return sorted(aln_results, key = lambda x: -x[4])

def print_aln(regs, ri, qi, bv, res, is_top):
    """
        pretty-printing for an alignment between 2 reads.
        regs: regions
        ri, qi: index of ref/query in regs
        bv: bit representation
        res: one alignment result from align_two()
    """

    def str_pair(e, d):
        # TODO: shorten this further? we may have dict somewhere?
        if e == 1:
            return "#" if d == 1 else ("." if d == -1 else "0")
        if e == -1:
            return "," if d == 1 else ("_" if d == -1 else "0")
        if e == 0:
            return "0" if d == 1 else ("0" if d == -1 else "@")

    # aligned part
    rs, re = res[0], res[1]
    qs, qe = res[2], res[3]
    # coord in bv matrix
    r0, r1 = regs[ri]
    q0, q1 = regs[qi]
    # lengths
    rl, ql = r1-r0, q1-q0

    #print(f"\nAlignment: {res[4]/10:.1f} % match in {res[5]} units.\n")
    print(f"\n     *     \tri\tqi\trs\tre\trl\tqs\tqe\tql\tscr\tlen")
    if is_top:
        print(f"BEST_ALIGN\t{ri}\t{qi}\t{rs}\t{re}\t{rl}\t{qs}\t{qe}\t{ql}\t{res[4]/10:.1f}\t{res[5]}\n")
    else:
        print(f"SUBOPT_ALIGN\t{ri}\t{qi}\t{rs}\t{re}\t{rl}\t{qs}\t{qe}\t{ql}\t{res[4]/10:.1f}\t{res[5]}\n")

    _l = f"statuses of {bv.shape[1]-3} SNVs"
    lines = f"  i\t       ids\tmid\t{_l:^40}\tmid\tids      \tj  "

    # dangling part if any
    for i in range(rs):
        lines += f"\n{i:>3}\t{bv[r0+i,1]:>5};{ri:>4}\t{bv[r0+i,2]:>3}\t" \
            + "".join([ ("x" if e == 1 else ("O" if e == 0 else " ")) for e in bv[r0+i,3:] ]) \
            + f"\t***\t*****;****\t***"

    for i in range(qs):
        lines += f"\n***\t*****;****\t***\t" \
            + "".join([ ("x" if e == 1 else ("O" if e == 0 else " ")) for e in bv[q0+i,3:] ]) \
            + f"\t{bv[q0+i,2]:<3}\t{bv[q0+i,1]:<5};{qi:<4}\t{i:<3}"

    # the aligned part
    for i in range(re-rs):
        lines += f"\n{rs+i:>3}\t{bv[r0+rs+i,1]:>5};{ri:>4}\t{bv[r0+rs+i,2]:>3}\t" \
            + "".join([ str_pair(e, d) for e,d in zip(bv[r0+rs+i,3:], bv[q0+qs+i,3:]) ]) \
            + f"\t{bv[q0+qs+i,2]:<3}\t{bv[q0+qs+i,1]:<5};{qi:<4}\t{qs+i:<3}"

    # the remaining dangling part if any
    for i in range(re, rl):
        lines += f"\n{i:>3}\t{bv[r0+i,1]:>5};{ri:>4}\t{bv[r0+i,2]:>3}\t" \
            + "".join([ ("x" if e == 1 else ("O" if e == 0 else " ")) for e in bv[r0+i,3:] ]) \
            + f"\t***\t*****;****\t***"

    for i in range(qe, ql):
        lines += f"\n***\t*****;****\t***\t" \
            + "".join([ ("x" if e == 1 else ("O" if e == 0 else " ")) for e in bv[q0+i,3:] ]) \
            + f"\t{bv[q0+i,2]:<3}\t{bv[q0+i,1]:<5};{qi:<4}\t{i:<3}"

    print(lines)
    # return lines

def gxfe_temp(nodes, edges):
    """ this directly outputs reads ovlp info in GEXF format. """

    s = '<?xml version="1.0" encoding="UTF-8"?>\n' +\
        '<gexf xmlns="http://www.gexf.net/1.2draft"' +\
        ' xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"' +\
        ' xsi:schemaLocation="http://www.gexf.net/1.2draft' +\
        ' http://www.gexf.net/1.2draft/gexf.xsd" version="1.2">\n' +\
        '<meta lastmodifieddate="1990-01-01"><creator>hacone@mlab</creator><description>reads ovlps</description></meta>\n' +\
        '<graph defaultedgetype="directed">\n'

    # TODO: abstract out?
    s += '<attributes class="node">' +\
        '<attribute id="0" title="n_units" type="integer"/>' +\
        '</attributes>\n'

    s += '<attributes class="edge">' +\
        '<attribute id="0" title="identity" type="float"/>' +\
        '<attribute id="1" title="units_ovlp" type="integer"/>' +\
        '<attribute id="2" title="step" type="integer"/>' +\
        '<attribute id="3" title="delta" type="float"/>' +\
        '</attributes>\n'

    s += "<nodes>\n" + nodes + "</nodes>\n"
    s += "<edges>\n" + edges + "</edges>\n"
    s += "</graph>\n</gexf>"

    return s

def align(hers, alns_pers = None):

    bag_of_units = [ (ri, h) for ri, er in enumerate(hers) for h, _, t in er.hors if t == "~" ]
    n_units = len(bag_of_units)
    print(f"{n_units} units found in {len(hers)} reads.")

    snv_sites = detect_snvs(bag_of_units)
    print(f"{len(snv_sites)} SNV sites defined.")

    brep = bit_vector_repr(hers, snv_sites, range_of_units = bag_of_units)
    print(f"{brep.shape[0] - n_units} are filled in construction of bits representation.")

    regs = sorted(valid_regions(brep), key = lambda x: x[0] - x[1])
    print(f"{len(regs)} consecutive regions found")

    # NOTE: I'm not going to change the definition of above variable for now; they should be stable.
    # NOTE: you may change below, though. 

    min_units = 10
    print(f"{ len([ r for r in regs if r[1]-r[0] >= min_units ]) } regions >= {min_units} units")

    import pickle
    if not alns_pers:
        pw_aln = {}
        for i in range(len(regs)):
            if regs[i][1] - regs[i][0] < min_units:
                continue
            pw_aln[i] = { j : align_two(brep[regs[i][0]:regs[i][1],:], brep[regs[j][0]:regs[j][1],:])
                for j in range(len(regs)) if i != j and regs[j][1]-regs[j][0] >= min_units }

            pw_aln[i] = { k : v for k, v in pw_aln[i].items() if v }
        pw_aln = { k : v for k, v in pw_aln.items() if v }

        pickle.dump(pw_aln, open(f"pairwise_alignments.{min_units}u.pickle", "wb"))
        print(f"{ sum([ len(pw_aln[i]) for i in pw_aln.keys() ]) } alignments saved to pairwise_alignments.{min_units}u.pickle")
    else:
        #pw_aln = pickle.load(open(f"pairwise_alignments.{min_units}u.pickle", "rb"))
        pw_aln = pickle.load(open(alns_pers, "rb"))
        print(f"{ sum([ len(pw_aln[i]) for i in pw_aln.keys() ]) } alignments loaded.")

    # print only aln with >4 units
    for i in pw_aln.keys():
        pw_aln[i] = { j : [ e for e in v if e[5] > 4 ] for j, v in pw_aln[i].items() }
        pw_aln[i] = { k : v for k, v in pw_aln[i].items() if v }
    pw_aln = { k : v for k, v in pw_aln.items() if v }

    nodes, edges = "", "" # for gephi

    #for i, v in list(pw_aln.items())[:10]:
    for i, v in list(pw_aln.items()):

        nodes += f'<node id="{i}" label="{i}"><attvalues><attvalue for="0" value="{regs[i][1]-regs[i][0]}"/></attvalues></node>\n'

        if len(v) < sum([ len(alns) for j, alns in v.items() ]):
            continue # only singles
            # pass # include doubles, triples, ...

        print(f"\nSTATS\t{i}\n" + "rid\tn_targets\tn_alns\n" + f"{i}\t{len(v)}\t{sum([ len(alns) for j, alns in v.items() ])}")

        # print out alignments for this reference.
        print("rid\tqid\tn_alns\tstep\tdelta\tIdent/units_ovlp")
        max_ident = None 
        for j, alns in [ (k, v) for k, v in sorted(list(v.items()), key = lambda x: -max([ y[4] for y in x[1]])) ][:10]:
            if not max_ident:
                max_ident = alns[0][4]

            #print(f"\nqid = {j} had {len(alns)} alns: {[s[4] for s in alns][:5]}...")
            step, delta = alns[0][0] - alns[0][2], alns[0][4] - max_ident

            print(f"{i}\t{j}\t{len(alns)}\t{step}\t{delta:.1f}\t" + " ".join([f"{s[4]/10:.1f}/{s[5]}" for s in alns]))

            # NOTE: some printing for graph viz. abstract out this later. # NOTE: only the best ovlp for the target, for now?
            if j in pw_aln.keys():
                #edges += f'<edge id="{i}-{j}" source="{i}" target="{j}" label="{step}"><attvalues>'
                edges += f'<edge id="{i}-{j}" source="{i}" target="{j}" label="{step}"><attvalues>' if step > 0 else f'<edge id="rev-{i}-{j}" source="{j}" target="{i}" label="{step}"><attvalues>'
                edges += f'<attvalue for="0" value="{alns[0][4]/10:.1f}" />' +\
                        f'<attvalue for="1" value="{alns[0][5]}" />' +\
                        f'<attvalue for="2" value="{step}" />' +\
                        f'<attvalue for="3" value="{delta:.1f}" />' +\
                        '</attvalues></edge>\n'

        # printing structure of alignments
        # TODO: consider a good ordering - containing, negative steps, contained (ns-ps), positive steps
        for j, alns in [ (k, v) for k, v in sorted(list(v.items()), key = lambda x: -max([ y[4] for y in x[1]])) ][:10]:
            ri, qi = i, j
            print_aln(regs, ri, qi, brep, alns[0], is_top = True) 
            for aln in alns[1:3]:
                print_aln(regs, ri, qi, brep, aln, is_top = False) 

    #print(gxfe_temp(nodes, edges), file=open("ovlp_graph.gexf", "w"))
    print(gxfe_temp(nodes, edges), file=open("ovlp_graph_step-direct.gexf", "w"))
    print("done.")

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='perform pair-wise alignment among HOR encoded reads on SNV data') # TODO: explain
    parser.add_argument('action', metavar='action', type=str, help='action to perform: align, ...')
    parser.add_argument('--hor-reads', dest='hors', help='pickled hor-encoded long reads')
    parser.add_argument('--alns', dest='alns', help='pickled pickled alignments')
    args = parser.parse_args()

    if args.action == "align":
        assert args.hors, "need HOR-encoded reads"
        hers = pickle.load(open(args.hors, "rb"))
        if args.alns:
            align(hers, args.alns)
        else:
            align(hers)
    else:
        assert None, "invalid action."
