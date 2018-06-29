#!/bin/bash

# extract kmer information from each cluster.
# this is to define alternative monomers merged in each HOR, and to define major deviation from the pattern

func() {
c=$1

for k in 2 3 4 5 6 8 10 12 15 20 25 30; do
# NOTE: this .py code will have been changed when you see this script. please infer my intension...
echo "python HOR_segregation.kmer.py kmer -k ${k} --reads encoded_read_clusters/C_${c}.pickle > ./kmer_raw_analysis/K${k}C${c}.tsv"
python HOR_segregation.kmer.py kmer -k ${k} --reads encoded_read_clusters/C_${c}.pickle > ./kmer_raw_analysis/K${k}C${c}.tsv

done
}; export -f func

for c in {1..40}; do echo $c; done | xargs -n1 -P10 -I% bash -c "func %"
