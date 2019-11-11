#!/bin/bash

# plot-length-distribution.sh target.fa.fai output.svg 50

FAI=$1
OUT=${2:-${1%%.fa.fai}.png}
RES=${3:-10}

gnuplot <<EOL
set terminal png
set output "${OUT}"
filter(x,y)=int(x/y)*y
plot "${FAI}" u (filter(\$2,${RES})):(1) smooth frequency with boxes
EOL
