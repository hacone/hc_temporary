#!/bin/bash

A=PCA-5000-units-50-vars.png

for B in $( ls PCA-read-50-vars-0*.png ); do
  composite -blend 40 $A $B .blend.$B
done

LAST=.blend.PCA-read-50-vars-0211.png

convert -layers optimize -loop 0 -delay 40 .blend.PCA-read-50-vars-0* -delay 460 $LAST animation.gif
rm .blend.*
