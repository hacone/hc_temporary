#!/bin/bash

## This performs monomer-assignment by blast. Output should be organized as sam/bam format file that can be directly parsed by python.
## As it considers monomers to be query, parsing part should be somewhat different. This specification might be changed later.

# Testing against the following version available in hx.
#
# makeblastdb: 2.4.0+
#  Package: blast 2.4.0, build Jun  1 2016 13:39:30
# blastn: 2.4.0+
#  Package: blast 2.4.0, build Jun  1 2016 13:39:30

export READ_DIR=/work2/hacone/human_centromeres/pacbio/Centro_Fasta_143Cells
prep_dbs() {
	zcat ${READ_DIR}/*.fasta.gz > pbreads_centromeres_143cells.fasta &
	makeblastdb -dbtype nucl -in pbreads_centromeres_143cells.fasta -parse_seqids
}

export MONOMER_DB=/work2/hacone/human_centromeres/cluster.s14.SRR3189741.fa
align() {
	MON_NAME=$1
	samtools faidx ${MONOMER_DB} ${MON_NAME} > queries_2/${MON_NAME}.fa
	blastn -db pbreads_centromeres_143cells.fasta -query queries_2/${MON_NAME}.fa -out blast_output_2/${MON_NAME} \
		-word_size 7 -qcov_hsp_perc 60 \
		-outfmt "17 SQ" -parse_deflines
}; export -f align

prep_dbs
samtools faidx ${MONOMER_DB}

mkdir -p queries_2/
mkdir -p blast_output_2/
cut -f1 ${MONOMER_DB}.fai | xargs -P 12 -I % bash -c "align %"
