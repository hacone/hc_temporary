#!/bin/bash

export SEQKIT=/home/hacone/local/bin/seqkit
export HCBIN=/work2/hacone/2018/human_centromeres
export PYTHON=$HCBIN/venv/bin/python3

function stats() {

        RP=$1 # path to mon-enc-reads.pickle
        DP=${RP%/*}
        FA=$( ls $DP | grep fasta$ )

        if [[ 1 == 0 ]]; then
                echo ""
        fi

        $PYTHON $HCBIN/monomer-composition.py  \
        --reads ${RP} --ref /data/hacone/10mons.lst \
        > ${RP%%.pickle}.compos

        samtools faidx ${DP}/${FA}

        wc -l ${DP}/${FA}.fai > ${RP%%.pickle}.stats
        wc -l ${RP%%.pickle}.compos >> ${RP%%.pickle}.stats
        gawk '$7>1&&$8>1&&$9>1&&$10>1&&$11>1' ${RP%%.pickle}.compos \
                | sort -k7,7nr | wc -l >> ${RP%%.pickle}.stats

        echo "done for RP=${RP}"

        rm ${RP%%.pickle}.typeW.fa 2> /dev/null
        gawk '$7>1&&$8>1&&$9>1&&$10>1&&$11>1&&NR>1{ print $1 }' ${RP%%.pickle}.compos \
        | xargs -I% samtools faidx ${DP}/${FA} % >> ${RP%%.pickle}.typeW.fa
        samtools faidx ${RP%%.pickle}.typeW.fa

}; export -f stats;

function print() {
        RP=$1
        $PYTHON $HCBIN/HOR_segregation.py print --reads ${RP} > ${RP%%.pickle}.monomers
}; export -f print;

# for a given mer.pickle, calculate hor-encoding & extract sequences out of reads.
function hor-seq-extract() {

        RP=$1 # path to mon-enc-reads.pickle
        DP=${RP%/*}
        TAB=${RP%%.pickle}.hor.txt

        # HOR encoding.
        $PYTHON $HCBIN/HOR_segregation.py encode-hor \
                --reads ${RP} --patterns ./HOR-pentamers.def --merged /dev/null \
                --out ${RP%%.pickle}.hor.pickle

        echo -n "e"

        $PYTHON $HCBIN/HOR_segregation.py print-hor \
                --hor-reads ${RP%%.pickle}.hor.pickle > ${TAB}

        echo -n "p"

        FA=${DP}/$( ls $DP | grep fasta$ ) # TODO: make sure this wont collide
        if [[ ! -e ${FA}.fai ]]; then
                samtools faidx ${FA}
        fi

        FORMATTER="$SEQKIT seq"
        for HOR in 11mW 5mW2 16mW 12mW2; do
                OUT=${RP%%.pickle}.${HOR}.fa
                rm ${OUT} 2> /dev/null
                while read line; do
                        set $line; R=$1; B=$2; E=$3
                        if [[ $B -lt $E ]]; then
                                RB=$(( $B - 5 ))
                                RE=$(( $E + 5 ))
                                if [[ $RB -gt 0 ]]; then
                                        samtools faidx ${FA} $R:$RB-$RE | $FORMATTER >> ${OUT}
                                fi
                        else
                                # TODO: offset should be 0 in newer version
                                RL=$( grep $R ${FA}.fai | cut -f2 )
                                RB=$(( $RL - (10000 - $E) - 5 ))
                                RE=$(( $RL - (10000 - $B) + 5 ))
                                if [[ $RB -gt 0 ]]; then
                                        samtools faidx ${FA} $R:$RB-$RE \
                                        | $FORMATTER -t dna --reverse --complement 2> /dev/null \
                                        | sed -e "/^>/s/$/_RC/" >> ${OUT}
                                fi
                        fi
                done < <(grep -e "$HOR" $TAB)
                if [[ -e $OUT ]]; then samtools faidx ${OUT} ; fi
        done

        echo -n "X"

}; export -f hor-seq-extract;

for sample in $( ls /glusterfs/hacone/blast-tmp/ ); do
        echo "start; sample = $sample" 
        find /glusterfs/hacone/blast-tmp/${sample} \
                | grep Hum14AlpMon | grep pickle | grep -v hor \
                | xargs -P 12 -I% bash -c "hor-seq-extract %"

                # | xargs -P 12 -I% bash -c "stats %"
                # | xargs -P 12 -I% bash -c "print %"
        echo -e "\ndone; sample = $sample" 
done

