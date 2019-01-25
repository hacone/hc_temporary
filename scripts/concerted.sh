
# Check if concerted evolution is observed.
python3 UnitAnalysis.py --hor-reads pickles/C_15.hor.pickle --skips 0,1 print-snv-evolution | LC_ALL=C sort -k1,1n > 0125.log.diff
for i in {1..20}; do gawk 'NR>2&&$1=='$i'{a+=$2;v+=$3;n+=1} END{print '$i', a/n, v/n, n}' 0125.log.diff ; done
