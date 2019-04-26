#!/bin/bash

PYTHON=../venv/bin/python3

if [[ 1 == 0 ]]; then
for s in 75 70 65 60; do
#for s in 80; do

$PYTHON ../Edges2Layout.py transitivity \
  --hor-reads ../2019Feb/C10.chrX.hor.pickle \
  --edges ../align0417/edges_best_score.dat \
  -o CHM13_CLR_s${s}p30.cons.pickle --prefix CHM13_CLR_s${s}p30 \
  --params ${s},30,-1 > .log-${s}

  #mkdir -p figures-${s}
  #mv *.png figures-${s}
  #mv *.svg figures-${s}
  rm *.png *.svg
  echo "done ${s}"

done
fi


# see alignments between layouts !!
PICKLES="CHM13_CLR_s60p30.cons.pickle,"
PICKLES=${PICKLES}"CHM13_CLR_s65p30.cons.pickle,"
PICKLES=${PICKLES}"CHM13_CLR_s70p30.cons.pickle,"
PICKLES=${PICKLES}"CHM13_CLR_s75p30.cons.pickle,"
PICKLES=${PICKLES}"CHM13_CLR_s80p30.cons.pickle"

for i in {0..23}; do
#for i in 0; do

$PYTHON ../UnitAnalysis.py layout \
  --consensi ${PICKLES} \
  --range 8,0 --park ${i} \
  --params 0.6,0.3,4 --err-rate 0.001 > cl.${i}.fwd.log &

$PYTHON ../UnitAnalysis.py layout \
  --consensi ${PICKLES} \
  --range 8,0 --park ${i} \
  --params 0.6,0.3,4 --err-rate 0.001 --reverse > cl.${i}.rev.log &

sleep 300

done
