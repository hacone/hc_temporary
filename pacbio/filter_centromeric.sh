# export READS_FQ=../pacbio/m160526_170733_42274_c101014002550000001823231410211647_s1_p0.filtered.subreads.fastq

# export READS_FQ=../pacbio/Fastq/m160426_225135_00116_c100976532550000001823226708101600_s1_p0.filtered.subreads.fastq.gz
export SQUEAKR_DIR=../squeakr
export READ_SQ_DIR=../read_squeakr

export FASTQ_DIR=Fastq
export TARGET_DIR=Centro_Fasta_143Cells

# shared parameters for squeakr
export K=6
export S=20

export MONS_CQF=./d0.fq.K$K.S$S.ser

make_ref_cqf() {

  # NOTE: executed on ax02
	${SQUEAKR_DIR}/squeakr-count -f -k $K -s $S -t 1 -o ./ ../data/monomers/d0.fq \
	&& mv d0.fq.ser ${MONS_CQF}
}


filter_centromeric() {

	READS_FQ=$1
	# TODO: check here whether the result is already exsited
	## This is the call for read-squeakr

	# -- for plain fq
	# ${READ_SQ_DIR}/squeakr-count -f -k $K -s $S -r ${MONS_CQF} -t 1 -o . ${READS_FQ} > ...

	# -- for fq.gz
	${READ_SQ_DIR}/squeakr-count \
        -g -k $K -s $S -r ${MONS_CQF} -t 1 -o . ${FASTQ_DIR}/${READS_FQ} \
        > ${TARGET_DIR}/${READS_FQ%%.fastq.gz}.fasta \
	&& gzip ${TARGET_DIR}/${READS_FQ%%.fastq.gz}.fasta

	echo "done for $READS_FQ"
}; export -f filter_centromeric

make_ref_cqf

mkdir -p $TARGET_DIR
ls $FASTQ_DIR | grep .fastq.gz \
| xargs -P 8 -I % bash -c "filter_centromeric %"

echo "all done"
