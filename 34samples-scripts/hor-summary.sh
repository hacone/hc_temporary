#/bin/bash

TAB=$1
DEF=$2

# synopsis: ./hor-summary.sh <TAB> <DEF>

nReads=$( cut -f1 $TAB | LC_ALL=C sort | uniq | grep -v ^$ | grep -v ^# | wc -l )
nMons=$( grep -v "#" $TAB | gawk '{ s+=$5 } END{ print s }' )

echo -e "Unit\tSize\tCount\tCount_OK\tC/kM\tM/kM\tC/kM*\tM/kM*"
echo -e "nMons\t1\t${nMons}\t*\t*\t*\t*\t*"
echo -e "nReads\t*\t${nReads}\t*\t*\t*\t*\t*"

while read line; do

	set $line

	US=$( echo $line | gawk '{ print NF - 1 }' )
	PAT=$( echo  $line | cut -d' ' -f1 )

	Count=$( grep -F "$PAT" $TAB | grep -v "#" | wc -l )
	Count_ok=$( grep -F "$PAT" $TAB | grep -v "#" | grep OK | wc -l )

  NORM_FCT=$nMons
  mNORM_FCT=$(( $nMons - ( $US - 1 ) * $nReads ))

  CKM=$( echo "scale=4; 1000.0 * $Count / $NORM_FCT" | bc -l )
  MKM=$( echo "scale=4; 1000.0 * $US * $Count / $NORM_FCT" | bc -l )
  CKMM=$( echo "scale=4; 1000.0 * $Count / $mNORM_FCT" | bc -l )
  MKMM=$( echo "scale=4; 1000.0 * $US * $Count / $mNORM_FCT" | bc -l )

	echo -e "${PAT}\t${US}\t${Count}\t${Count_ok}\t${CKM}\t${MKM}\t${CKMM}\t${MKMM}"

done < <( cat ${DEF} | grep -v ^$ | grep -v \# )

