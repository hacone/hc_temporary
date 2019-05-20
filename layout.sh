#!/bin/bash

source venv/bin/activate

for s in 60 65 70 75 80 85; do

  python3 Edges2Layout.py layout190512 \
    --hor-reads 2019Feb/C10.chrX.hor.pickle --edges align0417/edges_best_score.dat \
    --params ${s},30,0 --noc --tri > layout-log.${s}-noc-tri

  python3 Edges2Layout.py layout190512 \
    --hor-reads 2019Feb/C10.chrX.hor.pickle --edges align0417/edges_best_score.dat \
    --params ${s},30,0 --noc > layout-log.${s}-noc

done
