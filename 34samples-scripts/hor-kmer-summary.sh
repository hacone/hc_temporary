#!/bin/bash

echo -e "TYPE\tSAMPLE\tCAN\tSPOR\tpSPOR\tTAND\tpTAND\n"
for SAMPLE in $( cat ./HOR_SVs_All/samples.lst ); do
  TAB=HOR_SVs_All/pickles/${SAMPLE}.*16m*.3mers
  CAN=$( cat $TAB | grep -e "[^*]16m.*[^*]16m.*[^*]16m" | grep -v Rev | cut -f1 )
  SPOR=$( cat $TAB | grep -e "[^*]16m.*13m9-13.*[^*]16m" | grep -v Rev | cut -f1 )
  TAND=$( cat $TAB | grep -e "13m9-13.*13m9-13.*13m9-13" | grep -v Rev | cut -f1 )
  pSPOR=$( echo "scale=4; 100 * ${SPOR:-0} / $CAN" | bc -l )
  pTAND=$( echo "scale=4; 100 * ${TAND:-0} / $CAN" | bc -l )
  echo -e "16m-13m\t${SAMPLE}\t${CAN}\t${SPOR}\t${pSPOR}\t${TAND}\t${pTAND}"
done

echo -e "\nTYPE\tSAMPLE\tCAN\tSPOR\tpSPOR\tTAND\tpTAND\n"
for SAMPLE in $( cat ./HOR_SVs_All/samples.lst ); do
  TAB=HOR_SVs_All/pickles/${SAMPLE}.*5m*.3mers
  CAN=$( cat $TAB | grep -e "[^*]5m.*[^*]5m.*[^*]5m" | grep -v Rev | cut -f1 )
  SPOR=$( cat $TAB | grep -e "[^*]5m.*6m1.*[^*]5m" | grep -v Rev | cut -f1 )
  TAND=$( cat $TAB | grep -e "6m1.*6m1.*6m1" | grep -v Rev | cut -f1 )
  pSPOR=$( echo "scale=4; 100 * ${SPOR:-0} / $CAN" | bc -l )
  pTAND=$( echo "scale=4; 100 * ${TAND:-0} / $CAN" | bc -l )
  echo -e "5m-6m\t${SAMPLE}\t${CAN}\t${SPOR}\t${pSPOR}\t${TAND}\t${pTAND}"
done
