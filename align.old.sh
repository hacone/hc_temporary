#!/bin/bash

### Align short reads to monomer database.
###
### requires:
### 	fasta_formatter in FASTX Toolkit
###     bwa by Dr. Heng Li

if [[ 0 == 1 ]]; then # skip
echo "skipped"

function multiply() {
	seqfile=$1; mult=$2
	while read line; do
		if [[ ${line:0:1} == ">" ]]; then
			echo $line
		else
			seq $mult | xargs -I{} -P1 -n1 echo $line
		fi
	done < $seqfile
}

multiply data/monomers/MigaKH.HigherOrderRptMon.fa 1 | fasta_formatter -w 60 > data/monomers/single_mon.fa
#multiply data/monomers/MigaKH.HigherOrderRptMon.fa 2 | fasta_formatter -w 60 > data/monomers/double_mon.fa
#multiply data/monomers/MigaKH.HigherOrderRptMon.fa 3 | fasta_formatter -w 60 > data/monomers/triple_mon.fa

bwa index data/monomers/single_mon.fa
#bwa index data/monomers/double_mon.fa
#bwa index data/monomers/triple_mon.fa

fi # skip until here

ACC=$1

bwa mem data/monomers/single_mon.fa ./${ACC}.join.fq > ${ACC}.join.single.sam 2> ${ACC}.join.single.log
echo "done join single"
bwa mem data/monomers/single_mon.fa ./${ACC}.un1.fq > ${ACC}.un1.single.sam 2> ${ACC}.un1.single.log
echo "done un1 single"
bwa mem data/monomers/single_mon.fa ./${ACC}.un2.fq > ${ACC}.un2.single.sam 2> ${ACC}.un2.single.log
echo "done un2 single; all done"

