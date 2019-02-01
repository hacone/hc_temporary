#!/bin/bash

NVARS=25

A=PCA-9813-units-${NVARS}-vars.png

for B in $( ls PCA-read-${NVARS}-vars-*.png ); do
  composite -blend 40 $A $B .blend.$B
done

LAST=.blend.PCA-read-${NVARS}-vars-0199.png

convert -layers optimize -loop 0 -delay 40 .blend.PCA-read-${NVARS}-vars-* -delay 460 $LAST animation-${NVARS}-vars.gif
rm .blend.*
