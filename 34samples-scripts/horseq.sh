#!/bin/bash

export SEQKIT=/home/hacone/local/bin/seqkit
export FORMATTER="${SEQKIT} seq -t dna"

#  Remove \n => sed ':a;N;$!ba;s/\n//g'

# iter-consed in.fa out.fa seqname
function iter_consed() {

        IN_FA=$1
        IN_FA_BN=$( basename $IN_FA )
        OUT_FA=${2:-${IN_FA_BN%%.fa}}
        SEQNAME=${3:-${IN_FA_BN%%.fa}_consed}

        $SEQKIT sort -l $IN_FA | $FORMATTER -m 1 -w 0 -s > .tmp.all

        NSEQ=$( cat .tmp.all | wc -l )
        if [[ $NSEQ -gt 6250 ]]; then
                MINI=$( echo "($NSEQ / 2) - 2500" | bc )
        else
                MINI=$( echo "$NSEQ / 10" | bc )
        fi
        MAXI=$( echo "$NSEQ - $MINI" | bc )
        HALF=$( echo "($MINI + $MAXI) / 2" | bc )

        echo "nseq, mini, half, maxi = $NSEQ, $MINI, $HALF, $MAXI"

        gawk "$MINI<NR&&NR<$MAXI" .tmp.all > .tmp.80p

        rm ${OUT_FA}.single.fa 2> /dev/null
        rm ${OUT_FA}.double.fa 2> /dev/null

        for OF in {0..19}; do
                gawk "NR==($OF+$HALF)" .tmp.80p > .tmp.cs
                gawk "NR!=($OF+$HALF)" .tmp.80p >> .tmp.cs

                echo ">${SEQNAME}_s_${OF}" >> ${OUT_FA}.single.fa
                consed .tmp.cs | grep "^[acgt]" | tr "acgt" "ACGT" \
                        | sed ':a;N;$!ba;s/\n//g' >> ${OUT_FA}.single.fa
                echo -n $OF"."

                tail -n 1 ${OUT_FA}.single.fa > .tmp.cs2
                gawk "NR!=($OF+$HALF)" .tmp.80p >> .tmp.cs2

                echo ">${SEQNAME}_d_${OF}" >> ${OUT_FA}.double.fa
                consed .tmp.cs2 | grep "^[acgt]" | tr "acgt" "ACGT" \
                        | sed ':a;N;$!ba;s/\n//g' >> ${OUT_FA}.double.fa
                echo -n $OF"*"
        done

        # stats 20 single-round seqs
        $FORMATTER -m 1 -w 0 -s ${OUT_FA}.single.fa > .tmp.cs
        consed -V -t.04 .tmp.cs > ${IN_FA_BN%%.fa}.consed.log

        # stats 20 double-round seqs
        $FORMATTER -m 1 -w 0 -s ${OUT_FA}.double.fa > .tmp.cs
        consed -V -t.04 .tmp.cs >> ${IN_FA_BN%%.fa}.consed.log

        # wrap @ 60
        $FORMATTER -w 60 ${OUT_FA}.single.fa > .tmp.fa
        mv .tmp.fa ${OUT_FA}.single.fa

        $FORMATTER -w 60 ${OUT_FA}.double.fa > .tmp.fa
        mv .tmp.fa ${OUT_FA}.double.fa

        rm .tmp.{all,80p,cs,cs2,fa} 2> /dev/null
        echo -e "\ndone for ${IN_FA_BN}"

}; export -f iter_consed;

mkdir -p HOR-consensi
for s in $( ls /glusterfs/hacone/blast-tmp/ | grep -v Ashkenazi ); do
        for h in 11mW 5mW2 16mW 12mW2; do
                echo "iter_consed ./HOR-subseqs/$s.$h.fa $s.$h.cons ${s}_${h}"
                iter_consed ./HOR-subseqs/$s.$h.fa $s.$h.cons ${s}_${h}
                mv *.consed.log HOR-consensi/
                mv ${s}.${h}.double.fa HOR-consensi/
                mv ${s}.${h}.single.fa HOR-consensi/
        done
        echo "done for $s"
done
echo "all done"
