#!/bin/bash

# pick essential layouts

LAYOUT=log-layouts

sed -e "/^$/d" $LAYOUT \
| LC_ALL=C sort -k1,1n -k2,2n -k5,5nr -k4,4nr \
| gawk 'BEGIN{i="*";ai="*"} ($1!=i)||($2!=ai){print;i=$1;ai=$2;f=$5} $5<f{print $0"\tF"}' \
| sort -k6,6n > .forward-ext.log
echo $( cat .forward-ext.log | sed -e "/[FR]/d" | cut -f6 | uniq | wc -l )" forward extending layouts"

sed -e "/^$/d" $LAYOUT \
| LC_ALL=C sort -k1,1n -k2,2n -k4,4nr -k5,5nr \
| gawk 'BEGIN{i="*";ai="*"} ($1!=i)||($2!=ai){print;i=$1;ai=$2;r=$4} $4<r{print $0"\tR"}' \
| sort -k6,6n > .reverse-ext.log
echo $( cat .reverse-ext.log | sed -e "/[FR]/d" | cut -f6 | uniq | wc -l )" reverse extending layouts"

rm .essential-layouts 
for li in $(cat .reverse-ext.log .forward-ext.log | sed -e "/[FR]/d" | cut -f6 | sort -n | uniq ); do
  cat .reverse-ext.log .forward-ext.log \
  | gawk '$6=='$li | LC_ALL=C sort -k 6,6n -k3,3n | uniq >> .essential-layouts
done

exit

rm .forward-ext.log .reverse-ext.log .essential-layouts

sed -e "/^$/d" log-layouts | LC_ALL=C sort -k1,1n -k2,2n -k5,5nr -k4,4nr | gawk 'BEGIN{i="*";ai="*"} ($1!=i)||($2!=ai){print;i=$1;ai=$2}' | sort -k6,6n | cut -f6 | uniq | wc
sed -e "/^$/d" log-layouts | sort -k6,6n | cut -f6 | uniq | wc
sed -e "/^$/d" log-layouts | sort -k6,6n | less
