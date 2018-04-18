#!/bin/bash

# (hdf.tgz -> movie.bax) -> subreads.bam -> filtered.subreads.fq
mkdir -p HDFTGZ BAM Fastq

export SMRTLINK=/bio/package/pacbio/smrtlink/install/smrtlink-release_5.0.1.9585
export SMRTTOOLS=$SMRTLINK/bundles/smrttools/install/smrttools-release_5.0.1.9578

# TODO: Note: contains adhoc code.
down_dec_bam() {

	## Download hdf5.tgz
	url="http://sra-download.ncbi.nlm.nih.gov/srapub_files/${1}_${1}_hdf5.tgz"
	#echo "Downloading $url"
	#wget -c -P HDFTGZ $url

	## 22:27

	MOVIE=$( tar tzf HDFTGZ/${1}_${1}_hdf5.tgz | head -1 ); MOVIE=${MOVIE%%\.*}
	echo "Found Movie: $MOVIE"
	echo "Generating ${MOVIE}.subreads.bam"
	if [[ -e BAM/${MOVIE}.subreads.bam ]]; then
		echo "BAM/$MOVIE.subreads.bam is in"
	else
	## Extract from tar.gz and extract subreads
		# cd HDFTGZ && tar xf ${1}_${1}_hdf5.tgz && rm ${MOVIE}.bas.h5 ${MOVIE}.metadata.xml \
		cd HDFTGZ/ && $SMRTTOOLS/smrtcmds/bin/bax2bam ${MOVIE}.1.bax.h5 ${MOVIE}.2.bax.h5 ${MOVIE}.3.bax.h5 \
		&& mv ${MOVIE}.subreads.{bam,bam.pbi} ../BAM/ \
		&& rm ${MOVIE}.scraps.bam ${MOVIE}.scraps.bam.pbi ${MOVIE}.1.bax.h5 ${MOVIE}.2.bax.h5 ${MOVIE}.3.bax.h5
	fi

	## Cleanup hdf.tgz elsewhere, maybe

}; export -f down_dec_bam

bam2fastq() {
	MOVIE=$1
	echo "processing $MOVIE"

	$SMRTTOOLS/smrtcmds/bin/bamtools filter -length ">1000" -tag "rq:>0.75" -in BAM/${MOVIE}.subreads.bam \
	| $SMRTTOOLS/smrtcmds/bin/bamtools convert -format fastq -out Fastq/${MOVIE}.filtered.subreads.fastq \
	&& gzip Fastq/${MOVIE}.filtered.subreads.fastq

	gzip Fastq/${MOVIE}.filtered.subreads.fastq
}; export -f bam2fastq


# 3/30
# ACCESSION_LIST=SRX2010180.Runs.txt
# 4/12
ACCESSION_LIST=SRX2010823.Runs.txt

#cat $ACCESSION_LIST | xargs -P 10 -I % bash -c "down_dec_bam %"

# Movie id list
ls BAM/ | sed -e "s/.subreads.*//" | sort | uniq \
| xargs -P 10 -I % bash -c "bam2fastq %"

# Then you call filter_centromeric.sh

exit

# Other usage;
# $SMRTTOOLS/smrtcmds/bin/bamtools convert -format fasta -in BAM/${MOVIE}.subreads.bam -out Fastq/${MOVIE}.subreads.fasta

# this seems strange, but you need fastq for squeaker
# https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/fasta-to-fastq/fasta_to_fastq.pl
# perl fasta_to_fastq.pl ../data/monomers/d0.fa > ../data/monomers/d0.fq
