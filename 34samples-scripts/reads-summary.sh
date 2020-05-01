#!/bin/bash

echo -e "***\t"$( head -n1 ./filtered/CHM13CLR/*.stats | sed -e "s/  */\t/g" | cut -f4,5,7,8,13 )
for SAMPLE in $( cat ./HOR_SVs_All/samples.lst | grep -v -e "B[0-9][0-9][0-9]" -e "AK1" ); do
  echo -e ${SAMPLE}"\t"$( tail -n1 ./filtered/*${SAMPLE}*/*.stats | sed -e "s/  */\t/g"  | cut -f4,5,7,8,13 )
done
