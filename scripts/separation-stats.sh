#!/bin/bash

# in data.dat colomns are:
# 5: best score, 6: 2nd best score, 7: prom, 8: eov, 9: vars-used, 10: round

echo -e "T\tR\tB\tSB\tRatio\tETPN"

for T in 60 65 70 75 80 85 90 95 99 99.2 99.5 99.7 99.9; do
for R in {0..12}; do # TODO; max round?

R=$(( $R - 1 ))

cat data.dat | gawk '$10=='$R > .data.dat.$R

B=$( cat .data.dat.$R | gawk '$5>='$T | wc -l)
SB=$( cat .data.dat.$R | gawk '$6>='$T | wc -l)

if [[ $B -gt 0 ]]; then
        Ratio=$(echo "scale=3; 100*$SB/$B" | bc)
        echo -e "$T\t$R\t$B\t$SB\t$Ratio\t$(($B - $SB))"
fi

done

echo "" # new line

done


