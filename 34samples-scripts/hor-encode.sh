#!/bin/bash

HCROOT=/work2/hacone/2018/human_centromeres/
PYTHON=$HCROOT/venv/bin/python3

FOFN=$1

# TODO: parametrize HOR definition

while read line; do
        READS=$line
        $PYTHON $HCROOT/HOR_segregation.py encode-hor \
                --reads $READS --patterns ./HOR-pentamers.def --merged /dev/null \
                --out $FOFN.tmp.hor.enc

        # write into stdout
        $PYTHON $HCROOT/HOR_segregation.py print-hor \
                --hor-reads $FOFN.tmp.hor.enc
        rm $FOFN.tmp.hor.enc
done < <(cat $FOFN)
