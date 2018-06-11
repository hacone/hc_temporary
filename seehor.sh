#!/bin/bash

# TODO: accept params on hors

sed -e "s/horID_3.mon_2/horID_17.mon_0*/; s/horID_9.mon_10/horID_17.mon_0+/; s/horID_3.mon_3/horID_17.mon_1*/; s/horID_31.mon_3/horID_17.mon_1+/" $1
