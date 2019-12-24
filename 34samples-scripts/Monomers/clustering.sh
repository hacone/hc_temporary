#!/bin/bash

OUT=58mon.cluo

clustalo -i Hum58AlpMon.fa \
  --outfmt clu --resno --output-order tree-order --wrap 200 \
  -o ${OUT} \
  --guidetree-out ${OUT}.gt \
  --clustering-out ${OUT}.cls \
  --threads 12
