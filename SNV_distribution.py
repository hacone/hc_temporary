## This is for studying SNVs (mismatches) distribution in PacBio reads and short reads, over HOR-encoded regions or each monomer.
## Here is where HOR analysis and SNV analysis come together

# TODO: this is the most dirty code in this project. 
# TODO: rename this script

from collections import Counter
import re
import svgwrite
import sys
# for the other data structures
from EncodedRead import *
# for HOR_read
from HOR_segregation import *

# TODO: rename this
def load_array_of_default(hers):

    # TODO: to be defined elsewhere
    SNV_sites = [ (2, 15, "T"), (2, 135, "A"), (2, 139, "C"), (3, 108, "C"), (3, 59, "G"), (3, 116, "C"), (3, 137, "G"), (4, 92, "T"), (4, 122, "C"), (4, 5, "G"), (5, 45, "G"), (5, 109, "G"), (5, 65, "T"), (6, 15, "C"), (6, 137, "C"), (6, 41, "G"), (6, 89, "C"), (7, 127, "A"), (8, 105, "A"), (8, 75, "C"), (9, 138, "A"), (9, 148, "G"), (9, 163, "T"), (9, 78, "G"), (10, 101, "G"), (10, 43, "A"), (11, 22, "C") ]
    s20 = [(2, 135, "A"), (6, 15, "C"), (2, 15, "T"), (8, 105, "A"), (9, 138, "A"), (5, 65, "T"), (10, 43, "A"), (8, 75, "C"), (6, 41, "G"), (3, 116, "C"), (7, 127, "A"), (2, 139, "C"), (6, 89, "C"), (4, 5, "G"), (4, 122, "C")]
    s10 = [(10, 101, "G"), (9, 163, "T"), (4, 92, "T"), (5, 45, "G"), (6, 137, "C"), (5, 109, "G"), (9, 78, "G"), (3, 108, "C"), (3, 59, "G"), (3, 137, "G"), (9, 148, "G"), (11, 22, "C")]
    s5 = [(6, 36, "C"), (10, 137, "G"), (2, 90, "A"), (4, 25, "C"), (11, 68, "G"), (2, 154, "G"), (4, 35, "G"), (11, 154, "C"), (5, 3, "G"), (11, 12, "G"), (3, 132, "T"), (11, 39, "C"), (9, 162, "A")]

    ers = pickle.load(open(hers, "rb"))
    # pick the longest stretch of default unit
    for er in ers:

        arrays = []

        for h in filter(lambda x: x[2] == "~", er.hors):
            if arrays and h[0] == arrays[-1][1] + 12:
                arrays[-1][1] = h[0]
            else:
                arrays += [ [h[0], h[0]] ]

        if (not arrays) or all([ i == j for i,j in arrays ]):
            continue

        str_rep = [[], [], []]
        is_name_said = False

        # TODO: NOTE: this part output the string representation of the SNV-compressed reads, which was very essential. But functionally this should be elsewhere.
        for i,j in arrays:

            if not is_name_said:
                print(er.name)
                print(f"idx\t" + " ".join([ s[2] for s in s20 ]) + " | " + " ".join([ s[2] for s in s10 ]) + " | " + " ".join([ s[2] for s in s5 ]))
                is_name_said = True

            # TODO: unit size are given through param
            for s in range(i, j+1, 12):

                str_rep[0] = [ "." for i in range(len(s20)) ]
                for (l, (k, p, b)) in enumerate(s20):
                    if [ sp for sp in er.mons[s+k].monomer.snvs if sp.pos == p and sp.base == b ]:
                        str_rep[0][l] = "#"

                str_rep[1] = [ "." for i in range(len(s10)) ]
                for (l, (k, p, b)) in enumerate(s10):
                    if [ sp for sp in er.mons[s+k].monomer.snvs if sp.pos == p and sp.base == b ]:
                        str_rep[1][l] = "#"

                str_rep[2] = [ "." for i in range(len(s5)) ]
                for (l, (k, p, b)) in enumerate(s5):
                    if [ sp for sp in er.mons[s+k].monomer.snvs if sp.pos == p and sp.base == b ]:
                        str_rep[2][l] = "#"

                print(f"{er.mons[s].begin}\t" + " ".join(str_rep[0]) + " | " + " ".join(str_rep[1]) + " | " + " ".join(str_rep[2]))
            print("")
        
        # TODO: implement later
        #for i in range(len(er.hors)-1):
        #    if (er.hors[i][2] er.hors[i+1][2]):
                # make one

        # for i in [ h for h, _, t in er.hors if t == "~" ]:


# TODO: resolve all TODO.
def draw_all_default(hers, us, something_for_short_read = None):
    """
        generate SVG showing mismatch/SNV distribution over default unit (12-mer, for X).
        later this should be abstract out. args are:
        hers = HOR encoded reads
        us = unit size (can be inferred from hers[i].hors?)
    """

    # open long read 
    ers = pickle.load(open(hers, "rb"))
    n_units = 0
    snv_counters = [ Counter() for i in range(us) ]
    for er in ers:
        for i in [ h for h, _, t in er.hors if t == "~" ]:
            for j in range(us):
                snv_counters[j].update(er.mons[i+j].monomer.snvs)
            n_units += 1

    # NOTE: temporal TODO: move elsewhere # can be oneliner?
    #print("[ ")
    #for j in range(2, us):
    #    for k, v in [ (k, v) for k, v in snv_counters[j].items() if k != "Monomer" and (v / n_units > 0.05) ]:
    #        print(f"({j}, {k.pos}, \"{k.base}\"),\t{v/n_units}")
    #print("]")

    # open short read # TODO: short read pickle file should not be hardcoded
    vs = pickle.load(open("./variant_sites_50_0.01.pickle", "rb"))
    sr_snv = { re.sub("\.[0-9]+", "", k):v for k,v in list(vs.items()) }

    # draw # TODO: ideally output can be specified by args
    dwg = svgwrite.drawing.Drawing("./mms_X_default.svg", profile = "tiny")
    # TODO: legend for this color scheme?
    b2c = dict(A = "#F8766D", C = "#7CAE00", G = "#00BFC4", T = "#C77CFF")
    b2i = dict(A = 0, C = 1, G = 2, T = 3)

    xoff, yoff = 100, 50

    # X axis
    xx= dwg.g(id="x_axis")
    for i in [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
        for j in range(us):
            # TODO: scale factor (1000) should not be hardcoded.
            x =  xoff + (1000 * i / 100)
            xx.add(dwg.text(f"{i} %", insert = (x, yoff-10), font_size = "20"))
            xx.add(dwg.line((x, yoff + (750*j)), (x, yoff + (750*(j+1)) - 15), stroke = "#000000"))
            # for short read part
            x =  xoff + 1000 + (1000 * i / 100)
            xx.add(dwg.text(f"{i} %", insert = (x, yoff-10), font_size = "20"))
            xx.add(dwg.line((x, yoff + (750*j)), (x, yoff + (750*(j+1)) - 15), stroke = "#000000"))
    dwg.add(xx)

    for j in range(us):

        # mismatch set in pacbio long reads
        mms_pac = dwg.g(id=f"mms_pac_{j}")
        for s, c in snv_counters[j].items():
            print(f"{j}\t{s.pos}\t{s.base}\t{c}")
            mms_pac.add(dwg.rect(
                insert=(xoff, yoff + (750 * j) + (s.pos * 4) + b2i[s.base]), size=(1000 * c / n_units, 1), fill = b2c[s.base], fill_opacity = 0.8))
        dwg.add(mms_pac)

        # Y axis
        yx= dwg.g(id=f"y_axis_{j}")
        for i in [0, 50, 100, 150]:
            yx.add(dwg.text(f"{j}:{i}", insert = (40, yoff + (750*j) + (4*i)), font_size = "20"))
            yx.add(dwg.line((xoff-20, yoff + (750*j) + (4*i)), (xoff, yoff + (750*j) + (4*i)), stroke = "#000000"))
        dwg.add(yx)

    xoff = 1100
    # TODO: skipping 0, 1, is so special for X. I need recover these position somehow...
    for j in range(2, us):
        srj = sr_snv[f"horID_17.mon_{j}"]
        m_total = srj["Monomer"]
        mms_sr = dwg.g(id=f"mms_sr_{j}")
        for pb, c in srj.items():
            if pb == "Monomer":
                continue
            mms_sr.add(dwg.rect(
                insert=(xoff, yoff + (750 * j) + (pb[0] * 4) + b2i[pb[1]]), size=(1000 * c / m_total, 1), fill = b2c[pb[1]], fill_opacity = 0.8))
        dwg.add(mms_sr)

        # y axis for short read part
        yx= dwg.g(id=f"y_axis_sr_{j}")
        for i in [0, 50, 100, 150]:
            yx.add(dwg.text(f"{j}:{i}", insert = (xoff-60, yoff + (750*j) + (4*i)), font_size = "20"))
            yx.add(dwg.line((xoff-20, yoff + (750*j) + (4*i)), (xoff, yoff + (750*j) + (4*i)), stroke = "#000000"))
        dwg.add(yx)

    dwg.save()

#draw_all_default()

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description='Breakup encoded read based on the set of assigned monomers.')
    parser.add_argument('action', metavar='action', type=str, help='action to perform: plot, ...')
    parser.add_argument('--hor-reads', dest='hors', help='pickled hor-encoded long reads')
    parser.add_argument('--unit-size', dest='us', help='default unit size (currently not inferred)')
    args = parser.parse_args()

    if args.action == "plot":
        assert args.hors, "give me HOR-encoded reads"
        assert args.us, "give me unit size of the default unit"
        draw_al_default(args.hors, args.us)
    elif args.action == "":
        load_array_of_default()
        sys.exit()

