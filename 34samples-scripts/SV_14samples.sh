#!/bin/bash

#ls blast-tmp | transpose > all.16m.summary

HOR=${1:-16m}

# header.

( echo "*"; cut -f2 filtered_alignments/Ashkenazi.58mons.${HOR}.her.summary ) | transpose > .tmp
( echo "*"; cut -f1 filtered_alignments/Ashkenazi.58mons.${HOR}.her.summary ) | transpose >> .tmp

#for sample in $( ls blast-tmp ); do
for sample in Mende Esan Maasai; do
        (echo ${sample}; cut -f3 filtered_alignments/${sample}.58mons.${HOR}.her.summary) \
        | transpose >> .tmp.count
        (echo ${sample}; cut -f8 filtered_alignments/${sample}.58mons.${HOR}.her.summary) \
        | transpose >> .tmp.mpm
done
echo "" >> .tmp
echo "" >> .tmp.mpm

for sample in CHM13CLR CHM13Hifi NA12878Hifi Toscani Finnish Ashkenazi; do
        (echo ${sample}; cut -f3 filtered_alignments/${sample}.58mons.${HOR}.her.summary) \
        | transpose >> .tmp
        (echo ${sample}; cut -f8 filtered_alignments/${sample}.58mons.${HOR}.her.summary) \
        | transpose >> .tmp.mpm
done
echo "" >> .tmp
echo "" >> .tmp.mpm

for sample in Gujarati Dai HG005 ; do
        (echo ${sample}; cut -f3 filtered_alignments/${sample}.58mons.${HOR}.her.summary) \
        | transpose >> .tmp
        (echo ${sample}; cut -f8 filtered_alignments/${sample}.58mons.${HOR}.her.summary) \
        | transpose >> .tmp.mpm
done
echo "" >> .tmp
echo "" >> .tmp.mpm

for sample in PuertoRican Peruvian; do
        (echo ${sample}; cut -f3 filtered_alignments/${sample}.58mons.${HOR}.her.summary) \
        | transpose >> .tmp
        (echo ${sample}; cut -f8 filtered_alignments/${sample}.58mons.${HOR}.her.summary) \
        | transpose >> .tmp.mpm
done

cut -f2 --complement .tmp | column -t -L > all.${HOR}.count.summary
cut -f2 --complement .tmp.mpm | column -t -L > all.${HOR}.density.summary

rm .tmp .tmp.mpm

#  African
# Mende Esan Maasai
#  European
# CHM13CLR CHM13Hifi NA12878Hifi Toscani Finnish Ashkenazi
#  Asian
# Gujarati Dai HG005
#  American
# PuertoRican Peruvian

# and B cells...!
