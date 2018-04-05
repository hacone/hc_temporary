#!/bin/bash

SMRTLINK=/bio/package/pacbio/smrtlink/install/smrtlink-release_5.0.1.9585
SMRTTOOLS=$SMRTLINK/bundles/smrttools/install/smrttools-release_5.0.1.9578

MOVIE=m160526_170733_42274_c101014002550000001823231410211647_s1_p0

# $SMRTTOOLS/smrtcmds/bin/bax2bam ${MOVIE}.1.bax.h5 ${MOVIE}.2.bax.h5 ${MOVIE}.3.bax.h5 \
# && rm ${MOVIE}.scraps.bam ${MOVIE}.scraps.bam.pbi

#$SMRTTOOLS/smrtcmds/bin/bamtools convert -format fasta \
#	-in ${MOVIE}.subreads.bam -out ${MOVIE}.subreads.fasta

$SMRTTOOLS/smrtcmds/bin/bamtools filter -length ">1000" -tag "rq:>0.85" -in ${MOVIE}.subreads.bam \
| $SMRTTOOLS/smrtcmds/bin/bamtools convert -format fastq -out ${MOVIE}.filtered.subreads.fastq

# this seems strange, but you need fastq for squeaker
# https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/fasta-to-fastq/fasta_to_fastq.pl

# perl fasta_to_fastq.pl ../data/monomers/d0.fa > ../data/monomers/d0.fq
