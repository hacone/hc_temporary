#!/bin/bash

# (hdf.tgz -> movie.bax) -> subreads.bam -> filtered.subreads.fq

mkdir -p HDFTGZ BAM Fastq

down_dec_bam() {

	## Download hdf5.tgz
	url="http://sra-download.ncbi.nlm.nih.gov/srapub_files/${1}_${1}_hdf5.tgz"
	echo "downloading $url"
	wget -c -P HDFTGZ $url

	## Extract from tar.gz
	MOVIE=$( tar tzf HDFTGZ/${1}_${1}_hdf5.tgz | head -1 ); MOVIE=${MOVIE%%\.*}
	tar xf HDFTGZ/${1}_${1}_hdf5.tgz && rm HDFTGZ/${MOVIE}.bas.h5 HDFTGZ/${MOVIE}.metadata.xml 

	## Extract subreads
	SMRTLINK=/bio/package/pacbio/smrtlink/install/smrtlink-release_5.0.1.9585
	SMRTTOOLS=$SMRTLINK/bundles/smrttools/install/smrttools-release_5.0.1.9578

	cd HDFTGZ/ && $SMRTTOOLS/smrtcmds/bin/bax2bam ${MOVIE}.1.bax.h5 ${MOVIE}.2.bax.h5 ${MOVIE}.3.bax.h5 \
	&& mv ${MOVE}.subreads.bam BAM/ \
	&& rm ${MOVIE}.scraps.bam ${MOVIE}.scraps.bam.pbi ${MOVIE}.1.bax.h5 ${MOVIE}.2.bax.h5 ${MOVIE}.3.bax.h5

	## Cleanup hdf.tgz elsewhere

}; export -f down_dec_bam

ACCESSION_LIST=SRX2010180.Runs.txt
cat $ACCESSION_LIST | xargs -P 10 -I % bash -c "down_dec_bam %"

exit

# http://sra-download.ncbi.nlm.nih.gov/srapub_files/SRR4015715_SRR4015715_hdf5.tgz


#$SMRTTOOLS/smrtcmds/bin/bamtools convert -format fasta \
#	-in ${MOVIE}.subreads.bam -out ${MOVIE}.subreads.fasta

$SMRTTOOLS/smrtcmds/bin/bamtools filter -length ">1000" -tag "rq:>0.85" -in ${MOVIE}.subreads.bam \
| $SMRTTOOLS/smrtcmds/bin/bamtools convert -format fastq -out ${MOVIE}.filtered.subreads.fastq

# this seems strange, but you need fastq for squeaker
# https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/fasta-to-fastq/fasta_to_fastq.pl

# perl fasta_to_fastq.pl ../data/monomers/d0.fa > ../data/monomers/d0.fq
#!/bin/bash

## Memo : hdf5.tar.gz to fastq with adapter-trimming

tar xf 
