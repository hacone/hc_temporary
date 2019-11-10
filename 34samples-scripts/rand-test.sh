#!/bin/bash

for i in {1..1000}; do
  RAND=$( cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 12 | head -n 1 )
  echo -e $i"\t"$RAND
done
