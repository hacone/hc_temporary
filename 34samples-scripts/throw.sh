#!/bin/bash


if [[ 1 == 0 ]]; then

sleep 7200
./align_monomers.sh all filtered/PRJNA438669_Mende/split/ /glusterfs/hacone/blast-tmp/Mende/ 2>&1 > Mende.log &

sleep 7200
./align_monomers.sh all filtered/PRJNA480858_Peruvian/split/ /glusterfs/hacone/blast-tmp/Peruvian/ 2>&1 > Peruvian.log &

sleep 7200
./align_monomers.sh all filtered/PRJNA530212_Toscani/split/ /glusterfs/hacone/blast-tmp/Toscani/ 2>&1 > Toscani.log &


sleep 14400
./align_monomers.sh all filtered/PRJNA530214_Gujarati/split/ /glusterfs/hacone/blast-tmp/Gujarati/ 2>&1 > Gujarati.log &

sleep 14400
./align_monomers.sh all filtered/PRJNA530216_Esan/split/ /glusterfs/hacone/blast-tmp/Esan/ 2>&1 > Esan.log &

sleep 14400
./align_monomers.sh all filtered/PRJNA480712_Finnish/split/ /glusterfs/hacone/blast-tmp/Finnish/ 2>&1 > Finnish.log &

sleep 14400
./align_monomers.sh all filtered/PRJNA530217_Maasai/split/ /glusterfs/hacone/blast-tmp/Maasai/ 2>&1 > Maasai.log &

sleep 14400
./align_monomers.sh all filtered/PRJNA483067_PuertoRican/split/ /glusterfs/hacone/blast-tmp/PuertoRican/ 2>&1 > PuertoRican.log &

sleep 14400
./align_monomers.sh all filtered/PRJNA530776_CHM13Hifi/split/ /glusterfs/hacone/blast-tmp/CHM13Hifi/ 2>&1 > CHM13Hifi.log &

sleep 14400
./align_monomers.sh all filtered/PRJNA540705_NA12878Hifi/split/ /glusterfs/hacone/blast-tmp/NA12878Hifi/ 2>&1 > NA12878Hifi.log &
fi

# NOTE: 14mons alignment for the last three, start 11/04
# sleep 14400
#./align_monomers.sh all filtered/PRJNA547614_Dai/split/ /glusterfs/hacone/blast-tmp/Dai/ 2>&1 > Dai-14m.log &

#sleep 14400
./align_monomers.sh all filtered/PRJNA558394_Ashkenazi/split/ /glusterfs/hacone/blast-tmp/Ashkenazi/ 2>&1 > Ashkenazi-14m.log &
sleep 3600
./align_monomers.sh all filtered/SRX4739017_HG005/split/ /glusterfs/hacone/blast-tmp/HG005/ 2>&1 > HG005-14m.log &

exit
