#!/bin/bash

#sleep 3600
./align_monomers.sh hor-encode filtered/PRJNA438669_Mende/split/ /glusterfs/hacone/blast-tmp/Mende/ 2>&1 > Mende-hor-encode.log &

sleep 3600
./align_monomers.sh hor-encode filtered/PRJNA480858_Peruvian/split/ /glusterfs/hacone/blast-tmp/Peruvian/ 2>&1 > Peruvian-hor-encode.log &

sleep 3600
./align_monomers.sh hor-encode filtered/PRJNA530212_Toscani/split/ /glusterfs/hacone/blast-tmp/Toscani/ 2>&1 > Toscani-hor-encode.log &

sleep 3600
./align_monomers.sh hor-encode filtered/PRJNA530214_Gujarati/split/ /glusterfs/hacone/blast-tmp/Gujarati/ 2>&1 > Gujarati-hor-encode.log &

sleep 3600
./align_monomers.sh hor-encode filtered/PRJNA530216_Esan/split/ /glusterfs/hacone/blast-tmp/Esan/ 2>&1 > Esan-hor-encode.log &

sleep 3600
./align_monomers.sh hor-encode filtered/PRJNA480712_Finnish/split/ /glusterfs/hacone/blast-tmp/Finnish/ 2>&1 > Finnish-hor-encode.log &

sleep 3600
./align_monomers.sh hor-encode filtered/PRJNA530217_Maasai/split/ /glusterfs/hacone/blast-tmp/Maasai/ 2>&1 > Maasai-hor-encode.log &

sleep 3600
./align_monomers.sh hor-encode filtered/PRJNA483067_PuertoRican/split/ /glusterfs/hacone/blast-tmp/PuertoRican/ 2>&1 > PuertoRican-hor-encode.log &

sleep 3600
./align_monomers.sh hor-encode filtered/PRJNA530776_CHM13Hifi/split/ /glusterfs/hacone/blast-tmp/CHM13Hifi/ 2>&1 > CHM13Hifi-hor-encode.log &

sleep 3600
./align_monomers.sh hor-encode filtered/PRJNA540705_NA12878Hifi/split/ /glusterfs/hacone/blast-tmp/NA12878Hifi/ 2>&1 > NA12878Hifi-hor-encode.log &

sleep 3600
./align_monomers.sh hor-encode filtered/PRJNA547614_Dai/split/ /glusterfs/hacone/blast-tmp/Dai/ 2>&1 > Dai-hor-encode.log &

sleep 3600
./align_monomers.sh hor-encode filtered/PRJNA558394_Ashkenazi/split/ /glusterfs/hacone/blast-tmp/Ashkenazi/ 2>&1 > Ashkenazi-hor-encode.log &

sleep 3600
./align_monomers.sh hor-encode filtered/SRX4739017_HG005/split/ /glusterfs/hacone/blast-tmp/HG005/ 2>&1 > HG005-hor-encode.log &
