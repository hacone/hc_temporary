#/bin/bash

TAB=$1
DEF=$2

# synopsis: ./hor-summary.sh <TAB> <DEF>

#nreads=$( grep -v readname $TAB | grep -v ^$ | LC_ALL=C sort | uniq | wc -l )
nmons=$( grep -v readname $TAB | gawk '{ s+=$5 } END{ print s }' )

echo -e "unit\tsize\tocc\tmons"
echo -e "MONS\t1\t${nmons}\t${nmons}"

tot_pat=0

while read line; do

	set $line
	us=$( echo $line | gawk '{ print NF - 1 }' )
	pat=$( echo  $line | cut -d' ' -f1 )
	npat=$( grep -F "$pat" $TAB | wc -l)
	tot_pat=$(( $tot_pat + ( $us * $npat ) ))
	echo -e "${pat}\t${us}\t${npat}\t$(( $us * $npat ))"

done < <( cat ${DEF} | grep -v ^$ | grep -v \# )
