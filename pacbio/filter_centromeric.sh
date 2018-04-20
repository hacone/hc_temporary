# export READS_FQ=../pacbio/m160526_170733_42274_c101014002550000001823231410211647_s1_p0.filtered.subreads.fastq

# export READS_FQ=../pacbio/Fastq/m160426_225135_00116_c100976532550000001823226708101600_s1_p0.filtered.subreads.fastq.gz
export SQUEAKR_DIR=../squeakr
export READ_SQ_DIR=../read_squeakr

# shared parameters for squeakr
export K=6
export S=20

export MONS_CQF=./d0.fq.K$K.S$S.ser

make_ref_cqf() {

	${SQUEAKR_DIR}/squeakr-count -f -k $K -s $S -t 1 -o ./ ../data/monomers/d0.fq \
	&& mv d0.fq.ser ${MONS_CQF}
	# echo "got d0.fq.K$K.S$S.ser"
}


filter_centromeric() {

	READS_FQ=$1

	# TODO: check here whether the result is already exsited

	## This is the call for read-squeakr
	# -- for plain fq
	# ${READ_SQ_DIR}/squeakr-count -f -k $K -s $S -r ${MONS_CQF} -t 1 -o . ${READS_FQ}
	# -- for fq.gz
	${READ_SQ_DIR}/squeakr-count -g -k $K -s $S -r ${MONS_CQF} -t 1 -o . Fastq/${READS_FQ} > Centro_Fasta_143Cells/${READS_FQ%%.fastq.gz}.fasta \
	&& gzip Centro_Fasta_143Cells/${READS_FQ%%.fastq.gz}.fasta
	# there might be better way

	echo "done for $READS_FQ"

}; export -f filter_centromeric

mkdir -p Centro_Fasta_143Cells

ls Fastq/ | grep .fastq.gz \
| xargs -P 8 -I % bash -c "filter_centromeric %"

echo "all done"
