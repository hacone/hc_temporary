#!/bin/bash
CLS=$1
ACT=$2

## NOTE: suggested workflow;
## 1. you'll create an appropriate mental translating table of monomers in ./mons_dicts/ for each cluster,
##    by checking how kmer table will change accordingly (command: kmer)
## 2. you'll create table of possible HOR patterns in ./HOR_patterns/ (be careful for collision)
## 3. If you are satisfied, then ready, aim, encode! (command: encode)
## 3a. check the result by `grep -C 10 "symbol" ./encode.hor.cls${CLS}.dat | less -N` and/or command: summary; you may jump to 2 again.
## 4. Have you got good HOR structures?

# TODO: check input

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

# synopsis: ./HOR_work.sh 12 encode
if [[ $ACT == "encode" ]]; then
	python HOR_segregation.py encode-hor \
		--reads encoded_read_clusters/C_${CLS}.pickle \
		--merged mons_dicts/C${CLS}.dat \
		--patterns HOR_patterns/C${CLS}.dat \
		> encode.hor.cls${CLS}.dat
fi

# synopsis: ./HOR_work.sh 12 summary
if [[ $ACT == "summary" ]]; then
	nreads=$( grep s1_p0 encode.hor.cls${CLS}.dat | wc -l )
	echo -e "n_reads\t${nreads}"

	tot_pat=0
	for p in $( cut -f1 HOR_patterns/C${CLS}.dat | grep -v ^$ | grep -v \# ); do
		npat=$( grep $p encode.hor.cls${CLS}.dat | wc -l)
		echo -e "${p}\t${npat}"
		tot_pat=$(( $tot_pat + $npat ))
	done 

	nunenc=$( grep "M=" encode.hor.cls${CLS}.dat | wc -l )
	echo -e "n_unenc\t${nunenc} ($(( $nunenc / $nunenc + $tot_pat )) %)"
fi
