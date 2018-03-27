### Estimate monomers occurrence from SAM; using constrained gradient ascent.

from pysam import *

samfile = AlignmentFile("./SRR3189743.join.aligned.sort.sam", "r")
