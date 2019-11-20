#!/bin/bash

export SEQKIT=/home/hacone/local/bin/seqkit
export HCBIN=/work2/hacone/2018/human_centromeres
export PYTHON=$HCBIN/venv/bin/python3

export RESULT_DIR=/data/hacone/filtered_alignments/
mkdir -p $RESULT_DIR


# assuming 58 monomers alignments, calc monomer composition,
# which is then used for extract chromosome-identigied reads out of .sam.gz
function stats() {

        RP=$1 # path to mon-enc-reads.pickle
        SAMPLE=$2
        DP=${RP%/*}
        FA=$( ls $DP | grep fasta$ )

        if [[ 1 == 0 ]]; then
                echo ""
        fi

        $PYTHON $HCBIN/monomer-composition.py  \
        --reads ${RP} --ref /data/hacone/Monomers/58mons.lst \
        > ${RP%%.pickle}.compos

        if [[ ! -e ${DP}/${FA}.fai ]]; then
                samtools faidx ${DP}/${FA}
        fi

        #wc -l ${DP}/${FA}.fai > ${RP%%.pickle}.stats
        #wc -l ${RP%%.pickle}.compos >> ${RP%%.pickle}.stats
        rm ${RP%%.pickle}.stats 2> /dev/null

        # A: 5-mer(chr11), B: 11-mer(chr1), C: 12-mer(chrX), D: 16-mer(chr17)
        cat ${RP%%.pickle}.compos \
        | gawk '
                BEGIN{
                        OFS="\t";
                        print "#readname", "5m", "11m", "12m", "16m", "D", "J", "M", "R", "W"
                }
                NR>1{
                        A=$2+$3+$4+$5+$6;
                        B=$7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17;
                        C=$18+$19+$20+$21+$22+$23+$24+$25+$26+$27+$28+$29;
                        D=$30+$31+$32+$33+$34+$35+$36+$37+$38+$39+$40+$41+$42+$43+$44+$45;
                        E=$46+$47; F=$48+$49; H=$51+$52; I=$53+$54+$55+$56+$57;
                        print $1, A, B, C, D, E, F, $50, H, I
        }' > ${RP%%.pickle}.compos.sum

        RPDIR=$( dirname $RP )
        SAMGZ=$RPDIR/Hum58AlpMon.fa.sort.read.sam.gz
        while read line; do
                set $line
                if [[ $2 -gt 9 ]]; then
                        # NOTE: paths are hard-coded...
                        zcat $RPDIR/Hum58AlpMon.fa.sort.read.sam.gz \
                        | gawk '$1=="'${1}'"' >> $RPDIR/58mons.5m.sam
                fi
                if [[ $3 -gt 9 ]]; then
                        zcat $RPDIR/Hum58AlpMon.fa.sort.read.sam.gz \
                        | gawk '$1=="'${1}'"' >> $RPDIR/58mons.11m.sam
                fi
                if [[ $4 -gt 9 ]]; then
                        zcat $RPDIR/Hum58AlpMon.fa.sort.read.sam.gz \
                        | gawk '$1=="'${1}'"' >> $RPDIR/58mons.12m.sam
                fi
                if [[ $5 -gt 9 ]]; then
                        zcat $RPDIR/Hum58AlpMon.fa.sort.read.sam.gz \
                        | gawk '$1=="'${1}'"' >> $RPDIR/58mons.16m.sam
                fi
        done < <( grep -v readname ${RP%%.pickle}.compos.sum )

        echo "done for RP=${RP}"

}; export -f stats;

for sample in $( ls /glusterfs/hacone/blast-tmp/ ); do
        echo -e "\nstart; sample = $sample" 

        find /glusterfs/hacone/blast-tmp/${sample} \
                | grep Hum58AlpMon | grep pickle | grep -v hor \
                | xargs -P 12 -I% bash -c "stats % ${sample}"

        echo -e "done stats: $sample" 

        find /glusterfs/hacone/blast-tmp/${sample}  | grep Hum14AlpMon \
                | grep .compos$ | xargs cat > $RESULT_DIR/${sample}.10mons.compos

        find /glusterfs/hacone/blast-tmp/${sample}  | grep Hum58AlpMon \
                | grep .compos$ | xargs cat > $RESULT_DIR/${sample}.58mons.compos

        find /glusterfs/hacone/blast-tmp/${sample} \
                | grep .compos.sum$ | xargs cat > $RESULT_DIR/${sample}.58mons.compos.sum

        echo -e "got compos: $sample" 

        SAM=$( find /glusterfs/hacone/blast-tmp/${sample} | grep Hum58AlpMon | grep .sam.gz | head -n 1 )
        for hor in 5m 11m 12m 16m; do
                samtools view -H $SAM > $RESULT_DIR/${sample}.58mons.${hor}.sam
                find /glusterfs/hacone/blast-tmp/${sample} \
                | grep 58mons.${hor}.sam | xargs cat >> $RESULT_DIR/${sample}.58mons.${hor}.sam
                gzip $RESULT_DIR/${sample}.58mons.${hor}.sam &
                echo -e "compressing: $sample - $hor" 
        done
        echo -e "all done; sample = $sample" 
done

