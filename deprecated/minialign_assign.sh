#!/bin/bash

## monomer assignment to centromeric pacbio reads by minialign (check better params)
## the resulted sam will be parsed by pysam and then saved in a good format (id + var, alnstart-end, ...

#MINIALIGN=$(pwd)/minialign
# this one seeks unlimited number of seeds, i hope
MINIALIGN=$(pwd)/minialign_centromere/minialign

#REF=$(pwd)/data/monomers/d0.fa
REF=$(pwd)/cluster.s14.SRR3189741.fa

READS_GZ=$(pwd)/pacbio/Centro_Fasta/m160426_225135_00116_c100976532550000001823226708101600_s1_p0.filtered.subreads.fasta.gz

#DIR_FA_GZ=$(pwd)/pacbio/Centro_Fasta/
DIR_FA_GZ=$(pwd)/pacbio/Centro_Fasta_143Cells

# wc <(zcat pacbio/Centro_Fasta/*.fasta.gz)

#./minialign/minialign -xpacbio.clr -t4 -k15 ${REF} <(zcat ${READS_GZ}) > res.sam
#time ./minialign/minialign -xpacbio.clr -t8 -k8 ${REF} <(zcat ${DIR_FA_GZ}/*.fasta.gz) > nomd_k8_143cells.sam
#time ./minialign/minialign -TMD -xpacbio.clr -t8 -k8 ${REF} <(zcat ${DIR_FA_GZ}/*.fasta.gz) > md_k8_143cells.sam
#time ./minialign/minialign -TMD -xpacbio.clr -t8 -k8 -s25 -m0.15 ${REF} <(zcat ${DIR_FA_GZ}/*.fasta.gz) > md_k8_s25_m015_143cells.sam
#time ./minialign/minialign -TMD -xpacbio.clr -t8 -k7 -s25 -m0.15 ${REF} <(zcat ${DIR_FA_GZ}/*.fasta.gz) > md_k7_s25_m015_143cells.sam
#time ./minialign/minialign -TMD -xpacbio.clr -t8 -k6 -s25 -m0.15 ${REF} <(zcat ${DIR_FA_GZ}/*.fasta.gz) > md_k6_s25_m015_143cells.sam
#time ./minialign/minialign -TMD -xpacbio.clr -t8 -k7 -s10 -m0.10 ${REF} <(zcat ${DIR_FA_GZ}/*.fasta.gz) > ./data/pac_mon_aligned/md_k7_s10_m010_143cells.sam

time ${MINIALIGN} -TMD -xpacbio.clr -t8 -k8 -w4 -s50 -m0.30 ${REF} <(zcat ${DIR_FA_GZ}/*.fasta.gz) > ./data/pac_mon_aligned/Cent_k8_w4_s50_m30.sam

#for K in 8 7 6 5 4
#mkdir test_minialign
#time ./minialign/minialign -TMD -xpacbio.clr -t8 -k8 -s1 -m0.01 ${REF} <(zcat ${READS_GZ}) > 
