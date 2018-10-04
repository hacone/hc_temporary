#!/bin/bash

WORK_DIR=/work2/hacone/human_centromeres
MINIALIGN=$WORK_DIR/minialign

REF=$WORK_DIR/data/monomers/d0.fa
READS_GZ=$WORK_DIR/pacbio/Centro_Fasta/m160426_225135_00116_c100976532550000001823226708101600_s1_p0.filtered.subreads.fasta.gz

# DIR_FA_GZ=$(pwd)/pacbio/Centro_Fasta_143Cells

mkdir -p $WORK_DIR/test_minialign/

S=50
M=01

# ./minialign_centromere/minialign

#for K in 8 7; do
#for W in 5 4 3; do
#	time $WORK_DIR/minialign_centromere/minialign -TMD -TNM -xpacbio.clr -t8 -k${K} -w${W} -s${S} -m0.${M} ${REF} <(zcat ${READS_GZ}) \
#		> $WORK_DIR/test_minialign/Cent.k${K}.w${W}.s${S}.m${M}.sam 2> $WORK_DIR/test_minialign/minialign.log
#done; done

K=6
W=5
#time $WORK_DIR/minialign_centromere/minialign -TMD -TNM -xpacbio.clr -t8 -k${K} -w${W} -s${S} -m0.${M} ${REF} <(zcat ${READS_GZ}) \
# > $WORK_DIR/test_minialign/Cent.k${K}.w${W}.s${S}.m${M}.sam 2> $WORK_DIR/test_minialign/minialign.log

#time ./minialign/minialign -TMD -TNM -xpacbio.clr -t8 -k8 -s1 -m0.01 ${REF} <(zcat ${READS_GZ}) > 

OUTFILE=./minialign_coverage_errors_correct.dat

for f in $( ls ./test_minialign/*.sam ); do
	python3 show-reads/EncodedRead.py encode --sam $f --out ${f%%.sam}.pickle
	echo "pickled for $f"

	echo $f >> $OUTFILE
	echo "0.01" >> $OUTFILE
	python3 show-reads/EncodedRead.py correct --read ${f%%.sam}.pickle --vars ./variant_sites_50_0.01.pickle --out er01.pickle \
	| grep "TOTAL" >> $OUTFILE

	echo $f >> $OUTFILE
	echo "0.02" >> $OUTFILE
	python3 show-reads/EncodedRead.py correct --read ${f%%.sam}.pickle --vars ./variant_sites_50_0.02.pickle --out er02.pickle \
	| grep "TOTAL" >> $OUTFILE
done

