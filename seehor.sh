#!/bin/bash

# TODO: accept params on hors

#sed -e "s/horID_3.mon_2/horID_17.mon_0*/; s/horID_9.mon_10/horID_17.mon_0+/; s/horID_3.mon_3/horID_17.mon_1*/; s/horID_31.mon_3/horID_17.mon_1+/" $1
sed -e "s/3-2/17-0*/g; s/9-10/17-0+/g; s/3-3/17-1*/g; s/31-3/17-1+/g" $1
