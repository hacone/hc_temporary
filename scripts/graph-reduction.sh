#!/bin/bash

PYTHON=../venv/bin/python3

for s in 80 75 70 65 60; do

$PYTHON ../Edges2Layout.py transitivity \
  --hor-reads ../2019Feb/C10.chrX.hor.pickle \
  --edges ../align0417/edges_best_score.dat \
  --params ${s},30,-1 > .log-${s}

  mkdir -p figures-${s}
  mv *.png figures-${s}
  mv *.svg figures-${s}
  echo "done ${s}"

done
