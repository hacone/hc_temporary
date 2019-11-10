#!/bin/bash

# export HCBIN=/work2/hacone/2018/human_centromeres

function stats() {

        HCBIN=/work2/hacone/2018/human_centromeres

        RP=$1
        #echo "RP=${RP}"
        DP=${RP%/*}
        #echo "DP=${DP}"
        FA=$( ls $DP | grep fasta$ )
        #echo "FA="${FA}

        $HCBIN/venv/bin/python3 \
        $HCBIN/monomer-composition.py  \
        --reads ${RP} --ref /data/hacone/10mons.lst \
        > ${RP%%.pickle}.compos

        samtools faidx ${DP}/${FA}

        wc -l ${DP}/${FA}.fai > ${RP%%.pickle}.stats
        wc -l ${RP%%.pickle}.compos >> ${RP%%.pickle}.stats
        gawk '$7>1&&$8>1&&$9>1&&$10>1&&$11>1' ${RP%%.pickle}.compos \
                | sort -k7,7nr | wc -l >> ${RP%%.pickle}.stats

        echo "done for RP=${RP}"
        # TODO: actual extraction...

# gawk '$7>1&&$8>1&&$9>1&&$10>1&&$11>1' log.tmp | sort -k7,7nr | less

}; export -f stats;

function print() {
        HCBIN=/work2/hacone/2018/human_centromeres
        RP=$1
        $HCBIN/venv/bin/python3 $HCBIN/HOR_segregation.py print \
                --reads ${RP} > ${RP%%.pickle}.monomers
}; export -f print;

for sample in $( ls /glusterfs/hacone/blast-tmp/ ); do
        echo "start; sample = $sample" 
        find /glusterfs/hacone/blast-tmp/${sample} \
                | grep Hum14AlpMon | grep pickle \
                | xargs -P 12 -I% bash -c "print %"
                #| xargs -P 12 -I% bash -c "stats %"
        echo "done; sample = $sample" 
done
