export SQUEAKR_DIR=../squeakr
export READ_SQ_DIR=../read_squeakr
export FASTQ_DIR=Fastq

#export TARGET_DIR=Centro_Fasta_143Cells
export TARGET_DIR=stats_for_figures

# shared parameters for squeakr
export K=6
export S=20

#export MONS_CQF=./d0.fq.K$K.S$S.ser
export MONS_CQF=../resource/cqf-refs/d0.fq.K6.S20.ser

make_ref_cqf() {
  # NOTE: executed on ax02
	${SQUEAKR_DIR}/squeakr-count -f -k $K -s $S -t 1 -o ./ ../data/monomers/d0.fq \
	&& mv d0.fq.ser ${MONS_CQF}
}

filter_centromeric() {

	READS_FQ=$1
  IS_G=$2

	# TODO: check here whether the result is already exsited
	## This is the call for read-squeakr

  if [[ $IS_G == "true" ]]; then
    # -- for fq.gz
    ${READ_SQ_DIR}/squeakr-count \
          -g -k $K -s $S -r ${MONS_CQF} -t 1 -v 1 -o . ${FASTQ_DIR}/${READS_FQ} \
          > ${TARGET_DIR}/${READS_FQ%%.fastq.gz}.fasta \
    && gzip ${TARGET_DIR}/${READS_FQ%%.fastq.gz}.fasta

  else
	  # -- for plain fq
	  ${READ_SQ_DIR}/squeakr-count \
          -f -k $K -s $S -r ${MONS_CQF} -t 1 -o . ${FASTQ_DIR}/${READS_FQ} \
          > ${TARGET_DIR}/${READS_FQ%%.fastq}.fasta \
	  && gzip ${TARGET_DIR}/${READS_FQ%%.fastq}.fasta
  fi

	echo "done for $READS_FQ"

}; export -f filter_centromeric

#make_ref_cqf
#echo "Reference CQF generated"

mkdir -p $TARGET_DIR
ls -S $FASTQ_DIR | head -n 3 | grep .fastq.gz \
| xargs -P 3 -I % bash -c "filter_centromeric % true"

# To be verbose;
# ${READ_SQ_DIR}/squeakr-count \
#   -f -k $K -s $S -r ${MONS_CQF} -t 1 -o . -v 1 ${FASTQ_DIR}/${READS_FQ}

echo "all done"
