#/bin/bash

TAB=$1

# synopsis: ./hor-summary.sh

#nreads=$( grep -v readname $TAB | grep -v ^$ | LC_ALL=C sort | uniq | wc -l )
nreads=$( grep -v readname $TAB | gawk '{ s+=$5 } END{ print s }' )
echo $nreads

exit

tot_pat=0

echo -e "n_reads\t${nreads}"
echo -e "unit\tsize\tocc\tmons"

while read line; do

	set $line
	us=$( echo $line | gawk '{ print NF - 1 }' )
	pat=$( echo  $line | cut -d' ' -f1 )
	npat=$( grep -F "$pat" $TAB | wc -l)
	tot_pat=$(( $tot_pat + ( $us * $npat ) ))
	echo -e "${pat}\t${us}\t${npat}\t$(( $us * $npat ))"

done < <( cat HOR-pentamers.def | grep -v ^$ | grep -v \# )
