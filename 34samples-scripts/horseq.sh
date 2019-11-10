#!/bin/bash

export SEQKIT=/home/hacone/local/bin/seqkit

# This output raw read sequence which is assigned as a specified HOR unit
function extract-horseq() {

        FORMATTER="$SEQKIT seq -t dna"

        TAB=$1; READS=$2; HOR=$3

        echo "$( grep -e "$HOR" $TAB | wc -l ) $HOR seqs found." >&2
        while read line; do
                set $line; R=$1; B=$2; E=$3
                if [[ $B -lt $E ]]; then
                        RB=$(( $B - 10 ))
                        RE=$(( $E + 10 ))
                        samtools faidx ${READS} $R:$RB-$RE | $FORMATTER
                else
                        # TODO: revcomp version seems to produce slightly wrong seq?
                        RL=$( grep $R $READS.fai | cut -f2 )
                        RB=$(( $RL - (10000 - $E) - 10 ))
                        RE=$(( $RL - (10000 - $B) + 10 ))
                        samtools faidx ${READS} $R:$RB-$RE \
                        | $FORMATTER --reverse --complement | sed -e "/^>/s/$/_RC/"
                fi
        done < <(grep -e "$HOR" $TAB)
}

## 
# iter-consed in.fa out.fa seqname
function iter-consed() {
        #  Remove \n => sed ':a;N;$!ba;s/\n//g'

        IN_FA=$1
        OUT_FA=${2:-$IN_FA.cons.fa}
        SEQNAME=${3:-cons}

        # echo "seed: $(head -n 1 $IN_FA | sed -e 's/>//')"
        # $FORMATTER -m 1 -w 0 -s $IN_FA > .seqs.$IN_FA

        $FORMATTER -m 1 -w 0 -s $IN_FA > .tmp.seqs
        rm .cons.$IN_FA

        for N in {0..100}; do

        gawk 'NR>'$N .tmp.seqs > .seqs.$IN_FA
        head -n $N .tmp.seqs >> .seqs.$IN_FA
        consed .seqs.$IN_FA | grep "^[acgt]" | tr "acgt" "ACGT" \
                | sed ':a;N;$!ba;s/\n//g' >> .cons.$IN_FA

        cat .cons.$IN_FA .seqs.$IN_FA > .cons.seqs.$IN_FA

        echo ">$SEQNAME" > .cons2.$IN_FA
        consed .cons.seqs.$IN_FA | grep "^[acgt]" | tr "acgt" "ACGT" \
                | sed ':a;N;$!ba;s/\n//g' >> .cons2.$IN_FA

        $FORMATTER .cons2.$IN_FA > $OUT_FA
        rm .seqs.$IN_FA .cons.$IN_FA .cons.seqs.$IN_FA .cons2.$IN_FA

        done

}

TAB=$1; READS=$2; HOR=$3
FORMATTER="seqkit seq -t dna"
if [[ $TAB == "-h" ]]; then
        echo "horseq.sh TABLE READS HOR"; exit
fi
extract-horseq $TAB $READS $HOR > $HOR.fa
echo "extracted $HOR"

#iter-consed $HOR.fa $HOR.cons.fa
#echo "got consensus"
