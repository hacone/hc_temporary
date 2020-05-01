#!/bin/bash

HCROOT=/work2/hacone/2018/human_centromeres/
PYTHON=$HCROOT/venv/bin/python3
FOFN=./hor.pickles.fofn

while read line; do
        READS=$line
        echo "$PYTHON $HCROOT/HOR_segregation.py count-hor --hor-reads $READS > $READS.5mers"
        $PYTHON $HCROOT/HOR_segregation.py count-hor --hor-reads $READS > $READS.5mers
done < <(cat $FOFN)
