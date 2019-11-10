#!/bin/bash

for t in $(ls *.txt); do
        wc -w $t
        N=1
        for i in $(cat $t); do
                echo -e ${N}"\t"$(pwd)/$(find ncbitmp/sra/ | grep $i.fastq.stats) 
                N=$(( $N + 1 ))
        done
done > log-proc
