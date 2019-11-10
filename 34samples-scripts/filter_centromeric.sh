# export SQUEAKR_DIR=../squeakr
# export FASTQ_DIR=Fastq

## filterin-centromeres.sh <sample-name>.fq.fofn
# <sample-name>.fq.fofn : paths to fq or fq.gz files for the sample, separated by the newlines.
# this would create centromeric/<sample-name>/ and FASTA files for centromeric reads.

export HCBIN=/work2/hacone/2018/human_centromeres/
export SQBIN=/work2/hacone/2018/human_centromeres/read_squeakr
export FILT=$SQBIN/squeakr-count

# shared parameters for squeakr
export K=6 ; export S=20
export MONS_CQF=$HCBIN/resource/cqf-refs/Hum14AlpMon.K6S20.ser

#make_ref_cqf() {
  # NOTE: executed on ax02
	#${SQUEAKR_DIR}/squeakr-count -f -k $K -s $S -t 1 -o ./ ../resource/monomers/Hum14AlpMon.fq
  # then it was moved to $MONS_CQF
#}

filter_centromeric() {

	READS_FQ=$1 ; OUTDIR=$2
  if [[ $READS_FQ =~ \.gz$ ]]; then
          TF="-g"
  else
          TF="-f"
  fi

  # TODO: handle .zip file

  OUT=${READS_FQ%%.gz}
  OUT=${OUT%%.fastq}; OUT=${OUT%%.fq}
  OUT=${OUT%%.fasta}; OUT=${OUT%%.fa}
  OUT=${OUT}.centro.fa

	## This is the call for read-squeakr
  if [[ $3 == "dist" ]]; then
          # verbose mode to check distribution
          echo "here-vb; OUT=$OUT; READS_FQ=$READS_FQ"
          ${FILT} $TF -k $K -s $S -r ${MONS_CQF} -t 1 -u 200 -v 1 -o $OUTDIR ${READS_FQ} > /dev/null
  else
          # TODO: check here whether the result is already exsited
          echo "here2; OUT=$OUT; READS_FQ=$READS_FQ"
          ${FILT} $TF -k $K -s $S -r ${MONS_CQF} -t 1 -u 200 -o $OUTDIR ${READS_FQ} > $OUT \
                  && gzip $OUT \
                  && mv $OUT.gz $OUTDIR
  fi

	echo "filter done for $READS_FQ"

}; export -f filter_centromeric

#make_ref_cqf
#echo "Reference CQF generated"

FQ_FOFN=$1
SAMPLE=${FQ_FOFN%%.fq.fofn}
OUTDIR=$(pwd)/filtered/$SAMPLE
mkdir -p $OUTDIR

if [[ $2 == "dist" ]]; then
        head -n 10 $FQ_FOFN | xargs -P 10 -I % bash -c "filter_centromeric % $OUTDIR dist"
        cat $OUTDIR/*.refip > .tmp.${SAMPLE}.10runs.refip
        cut -f3,4 .tmp.${SAMPLE}.10runs.refip | gawk '$1>=1000' \
                | LC_ALL=C sort -k2,2nr > $OUTDIR/${SAMPLE}.10runs.refip
        rm .tmp.${SAMPLE}.10runs.refip $OUTDIR/*.refip
else
        if [[ ! -e $OUTDIR/${SAMPLE}.stats ]]; then
          echo "run filter & stats..."
          cat $FQ_FOFN | xargs -P 12 -I % bash -c "filter_centromeric % $OUTDIR"
          zcat $OUTDIR/*.gz | seqkit stats -a > $OUTDIR/${SAMPLE}.stats
          zcat $OUTDIR/*.gz | seqkit seq -m 1000 | seqkit stats -a >> $OUTDIR/${SAMPLE}.stats
        fi

        if [[ ! -e $OUTDIR/split ]]; then
          ## split into 10Mb chunks in $OUTDIR/split
          AVG_RL=$( gawk 'NR==4{print $7}' $OUTDIR/${SAMPLE}.stats | sed -e "s/,//g" )
          NSEQ=$( echo "scale=0; 10 * 1000 * 1000 / $AVG_RL" | bc -l )
          echo $AVG_RL" bps avg. => splitting into "$NSEQ" seqs / file"

          CDIR=$(pwd)
          cd $OUTDIR
          seqkit seq -m 1000 *.fa.gz | seqkit split -s $NSEQ
          mv stdin.split split
          for f in $(ls split/*.fasta); do
                nf=$(echo $f | sed -e "s/stdin/$SAMPLE/")
                mv $f $nf
          done
          find split/ | grep .fasta | xargs -P 12 gzip
          cd $CDIR
        fi
fi

echo "done"

