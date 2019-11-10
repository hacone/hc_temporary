#!/bin/bash

BBMAP=/home/hacone/tools/bbmap/current

for FQ in $(cat readtypes.dat | grep Raw | cut -f1 | sed -e "s#//#/#g" | head -n 5); do
        OUT=${FQ%%.fastq.gz}.subreads.fastq
        LOG=${FQ%%.fastq.gz}.subreads.raw.log
        echo "java -ea -Xmx2g -cp $BBMAP pacbio.RemoveAdapters2 in=$FQ out=$OUT &" | tee -a chop.log
        java -ea -Xmx2g -cp $BBMAP pacbio.RemoveAdapters2 in=$FQ out=$OUT 2>&1 > $LOG &
done

for FQ in $(cat readtypes.dat | grep -v Raw | grep -v CCS | cut -f1 | sed -e "s#//#/#g" | head -n 5); do
        OUT=${FQ%%.fastq.gz}.subreads.fastq
        LOG=${FQ%%.fastq.gz}.subreads.sr.log
        echo "java -ea -Xmx2g -cp $BBMAP pacbio.RemoveAdapters2 in=$FQ out=$OUT &" | tee -a chop.log
        java -ea -Xmx2g -cp $BBMAP pacbio.RemoveAdapters2 in=$FQ out=$OUT 2>&1 > $LOG &
done
