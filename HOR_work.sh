#!/bin/bash
export CLS=$1
export ACT=$2

## NOTE: suggested workflow;
## 1. you'll create an appropriate mental translating table of monomers in ./mons_dicts/ for each cluster,
##    by checking how kmer table will change accordingly (command: kmer)
## 2. you'll create table of possible HOR patterns in ./HOR_patterns/ (be careful for collision)
## 3. If you are satisfied, then ready, aim, encode! (command: encode)
## 3a. check the result by `grep -C 10 "symbol" ./encode.hor.cls${CLS}.dat | less -N` and/or command: summary; you may jump to 2 again.
## 4. Have you got good HOR structures?

# NOTE: commands implemented are: kmer, raw-kmer, show-def, encode, summary, gap-summary

# TODO: check input
# NOTE: tsv in kmer_raw_analysis was created using kmer_analysis.sh

# synopsis: ./HOR_work.sh 12 kmer 15 | less -N
# available K = 2, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30
if [[ $ACT == "kmer" ]]; then
	export LC_ALL=C
	K=$3
	python kmer_analysis.py mons_dicts/C${CLS}.dat kmer_raw_analysis/K${K}C${CLS}.tsv | sort -k1,1nr
fi

# synopsis: ./HOR_work.sh 12 raw-kmer 15 | less -N
if [[ $ACT == "raw-kmer" ]]; then
	export LC_ALL=C
	K=$3
	cat kmer_raw_analysis/K${K}C${CLS}.tsv | sort -k1,1nr
fi

if [[ $ACT == "show-def" ]]; then
	echo "Identifications made:"
	cat ./mons_dicts/C${CLS}.dat

	echo -e "\nPatterns defined:"
	cat ./HOR_patterns/C${CLS}.dat
fi

# synopsis: ./HOR_work.sh 12 encode
if [[ $ACT == "encode" ]]; then
	mkdir -p HOR_encoded_reads/
	python HOR_segregation.py encode-hor \
		--reads encoded_read_clusters/C_${CLS}.pickle \
		--merged mons_dicts/C${CLS}.dat \
		--patterns HOR_patterns/C${CLS}.dat \
		--out HOR_encoded_reads/C_${CLS}.hor.pickle

fi

if [[ $ACT == "print-hor" ]]; then
	python HOR_segregation.py print-hor \
		--hor-reads HOR_encoded_reads/C_${CLS}.hor.pickle \
		> encode.hor.cls${CLS}.dat
fi

# This output raw read sequence which is assigned as a specified HOR unit
if [[ $ACT == "show-hor" ]]; then
#	NOTE: delegated to print-hor
#	python HOR_segregation.py show \
#		--hor-reads HOR_encoded_reads/C_12.hor.pickle > log
		#--merged mons_dicts/C${CLS}.dat \
		#--patterns HOR_patterns/C${CLS}.dat \
		#--out HOR_encoded_reads/C_${CLS}.hor.pickle

	FORMATTER=~/local/bin/fasta_formatter
	PAT=$3
	REF_PBREADS=pacbio/blast/pbreads_centromeres_143cells.fasta
	HOR_TAB=./encode.hor.cls${CLS}.dat

	if [[ $PAT == "default" ]]; then
		while read line; do
			# TODO: what to do with invalid regions?
			set $line;
			if [[ $2 -lt $3 ]]; then
				samtools faidx ${REF_PBREADS} ${1}:$(($2 - 50))-$(($3 + 50));
			else
				samtools faidx ${REF_PBREADS} ${1}:$(($3 - 50))-$(($2 + 50)) \
				| ${FORMATTER} -o .tmp
				cat <( head -n1 .tmp | sed -e "/^>/s/$/_RC/" ) \
				    <( gawk 'NR>1' .tmp | sed -e "/^[ACGTNacgtn]/y/acgtACGT/tgcaTGCA/" | rev ) \
			        | ${FORMATTER} -w 60
				rm .tmp
			fi
		done < <(grep -e "~" $HOR_TAB | grep -v "\]" | grep -v "\[")
	else
		while read line; do
			set $line;
			if [[ $2 -lt $3 ]]; then
				samtools faidx ${REF_PBREADS} ${1}:$(($2 - 50))-$(($3 + 50)) \
				| ${FORMATTER} -w 60 
			else
				# TODO: revcomp version seems to produce slightly wrong seq?
				samtools faidx ${REF_PBREADS} ${1}:$(($3 - 50))-$(($2 + 50)) \
				| ~/local/bin/fasta_formatter -o .tmp
				cat <( head -n1 .tmp | sed -e "/^>/s/$/_RC/" ) \
				    <( gawk 'NR>1' .tmp | sed -e "/^[ACGTNacgtn]/y/acgtACGT/tgcaTGCA/" | rev ) \
			        | ${FORMATTER} -w 60
				rm .tmp
			fi
		done < <(grep -e $PAT $HOR_TAB)
	fi
fi


# synopsis: ./HOR_work.sh 12 summary
if [[ $ACT == "summary" ]]; then

	nreads=$( grep s1_p0 encode.hor.cls${CLS}.dat | wc -l )
	echo -e "n_reads\t${nreads}"

	tot_pat=0

	echo -e "unit\tsize\tocc\tmons"
	while read line; do
		set $line
		us=$( echo $line | gawk '{ print NF - 1 }' )
		pat=$( echo  $line | cut -d' ' -f1 )
		npat=$( grep -F "$pat" encode.hor.cls${CLS}.dat | wc -l)
		tot_pat=$(( $tot_pat + ( $us * $npat ) ))
		echo -e "${pat}\t${us}\t${npat}\t$(( $us * $npat ))"

	done < <( cat HOR_patterns/C${CLS}.dat | grep -v ^$ | grep -v \# )

	#for p in $( cut -f1 HOR_patterns/C${CLS}.dat | grep -v ^$ | grep -v \# ); do
	#	npat=$( grep -F "$p" encode.hor.cls${CLS}.dat | wc -l)
	#	echo -e "${p}\t${npat}"
	#	tot_pat=$(( $tot_pat + $npat ))
	#done 

	nunenc=$( grep "M=" encode.hor.cls${CLS}.dat | wc -l )
	echo -e "n_enc'ed\t*\t${tot_pat} ("$( echo "scale=2; 100*${tot_pat}/(${nunenc}+${tot_pat})" | bc -l )" %)"
	echo -e "n_unenc'ed\t*\t${nunenc} ("$( echo "scale=2; 100*${nunenc}/(${nunenc}+${tot_pat})" | bc -l )" %)"

fi

if [[ $ACT == "gap-summary" ]]; then

	# TODO; this works only for cls12
	for t in {0..9}; do
		echo -e "\n$t"
		cat encode.hor.cls${CLS}.dat \
		| sed -ne "/s1_p0/h; /\~${t}\]/{ N; /\[$(( $t + 2 ))\~/{ x;p;x;p } }"
		#| sed -ne "/s1_p0/h; /\~${t}\]/{ N; /\[$(( $t + 1 ))\~/{ x;p;x;p } }"
	done
fi
