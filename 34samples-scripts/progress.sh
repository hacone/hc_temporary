#!/bin/bash

echo -e "\nno. of CHUNKS."
for dir in $(ls filtered/ ); do echo -e "${dir}\t"$( ls filtered/$dir/split/ | wc -l ); done

#echo -e "\nfinished ENCODING"
#echo -e "> For s14-X17"
#for s in $(ls *-s14-encode.log); do echo $s ; grep Done $s | wc -l ; done
#echo -e "> For 14 mons"
#for s in $(ls *-14m-encode.log); do echo $s ; grep Done $s | wc -l ; done

for sample in $( ls /glusterfs/hacone/blast-tmp/ ); do

        echo "Alignment for $sample : "$( ls -lah /glusterfs/hacone/blast-tmp/$sample/*/*.sort.read.sam.gz | wc -l )
        ls -lahS /glusterfs/hacone/blast-tmp/$sample/*/*.sort.read.sam.gz

        echo "Encoding for $sample : "$( ls -lah /glusterfs/hacone/blast-tmp/$sample/*/*.pickle | wc -l )
        ls -lahS /glusterfs/hacone/blast-tmp/$sample/*/*.pickle
done
