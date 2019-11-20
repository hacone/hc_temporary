#!/bin/bash

HCROOT=/work2/hacone/2018/human_centromeres/
PYTHON=$HCROOT/venv/bin/python3

FOFN=$1
HORDEF=$2
HORDEF_BN=$( basename ${HORDEF%%.def} )

while read line; do
        READS=$line
        $PYTHON $HCROOT/HOR_segregation.py encode-hor \
                --reads $READS --patterns ${HORDEF} --merged /dev/null \
                --out $FOFN.${HORDEF_BN}.tmp.hor.enc
        # write into stdout
        $PYTHON $HCROOT/HOR_segregation.py print-hor \
                --hor-reads $FOFN.${HORDEF_BN}.tmp.hor.enc
        rm $FOFN.${HORDEF_BN}.tmp.hor.enc
done < <(cat $FOFN)
