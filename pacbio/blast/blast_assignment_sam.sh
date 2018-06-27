#!/bin/bash

export TRIAL=3
export READ_DIR=/work2/hacone/human_centromeres/pacbio/Centro_Fasta_143Cells
export MONOMER_DB=/work2/hacone/human_centromeres/cluster.s14.SRR3189741.fa

## This performs monomer-assignment by blast. Output should be organized as sam/bam format file that can be directly parsed by python.
## As it considers monomers to be query, parsing part should be somewhat different. This specification might be changed later.

# Testing against the following version available in hx.
#
# makeblastdb: 2.4.0+
#  Package: blast 2.4.0, build Jun  1 2016 13:39:30
# blastn: 2.4.0+
#  Package: blast 2.4.0, build Jun  1 2016 13:39:30

prep_dbs() {
	zcat ${READ_DIR}/*.fasta.gz > pbreads_centromeres_143cells.fasta &
	makeblastdb -dbtype nucl -in pbreads_centromeres_143cells.fasta -parse_seqids
}

align() {
	MON_NAME=$1
	samtools faidx ${MONOMER_DB} ${MON_NAME} > queries_${TRIAL}/${MON_NAME}.fa
	blastn -db pbreads_centromeres_143cells.fasta -query queries_${TRIAL}/${MON_NAME}.fa -out blast_output_${TRIAL}/${MON_NAME} \
		-max_target_seqs 1000000 -word_size 7 -qcov_hsp_perc 60 \
		-outfmt "17 SQ" -parse_deflines
}; export -f align

split_by_movie() {
# TODO: This should be coupled with align() later.
# TODO: Use of /tmp should be applied for other parts as well...?
	MON_NAME=$1
	TARGET_DIR=./tmp/blast_assignment/by_monomer/${MON_NAME}
	mkdir -p $TARGET_DIR
	head -n 3 blast_output_${TRIAL}/${MON_NAME} > ${TARGET_DIR}/header
	cat blast_output_${TRIAL}/${MON_NAME} | gawk 'NR>3 { dest=$1; gsub(/\/.*/, "", dest); print >> "'${TARGET_DIR}/'"dest }'

	# TODO this is too dangerous
	# ls $TARGET_DIR/ | grep -v header >> .tmp.movie.lst
}; export -f split_by_movie

collect_by_movie() {

	TMP_DIR=$(pwd)/tmp/blast_assignment/
	mkdir -p ${TMP_DIR}/by_movie/

	RESULT_DIR=./blast_assignment/
	mkdir -p ${RESULT_DIR}


# Generate sam for each movie id
	#while read movie; do
	movie=$1

		echo "processing for $movie"
		cat ${TMP_DIR}/header > ${TMP_DIR}/by_movie/${movie}.sam
		find ${TMP_DIR}/by_monomer/ | grep $movie | sort -R | xargs cat >> ${TMP_DIR}/by_movie/${movie}.sam
		echo "collected records for $movie"

		TARGET=${TMP_DIR}/by_movie/${movie}
# Convert merged.sam to bam with MD tag. # TODO: parametrize reference monomer seq
		samtools view -@ 12 -b -F 4 ${TARGET}.sam > ${TARGET}.bam && \
		samtools calmd -@ 12 -b ${TARGET}.bam ${MONOMER_DB} > ${TARGET}.md.bam 2> /dev/null && \
		rm ${TARGET}.bam
		echo "md tags calculated for ${TARGET}.sam as .md.bam"

# Convert merged sam to bam sorted in reference, then index it.
		samtools sort -@ 12 -o ${TARGET}.sort.bam ${TARGET}.md.bam && \
		samtools index ${TARGET}.sort.bam
		echo "got bam(.bai) sorted by ref (monomers)."

# Convert merged sam to bam sorted in read ID, then get sam to be parsed.
		samtools sort -@ 12 -n -o ${TARGET}.sort.read.bam ${TARGET}.md.bam && \
		samtools view -h ${TARGET}.sort.read.bam > ${TARGET}.sort.read.sam
		gzip ${TARGET}.sort.read.sam
		echo "got sam sorted by read. ready to be parsed. deleting temp files..."

# delete temp files and move the results
		rm ${TARGET}.sort.read.bam ${TARGET}.md.bam
		mv ${TARGET}.sort.bam ${RESULT_DIR}
		mv ${TARGET}.sort.bam.bai ${RESULT_DIR}
		mv ${TARGET}.sort.read.sam.gz ${RESULT_DIR}
		echo "done for $movie."
		# NOTE: in the end, you'll have 1 sam (read-sorted gzipped) and 2 bams (ref-sorted and read-sorted).
	#done < <(ls ${TMP_DIR}/by_monomer/*/ | grep s1_p0 | sed -e "s/.*\///" | sort | uniq) 
	#done < <( sort .tmp.movie.lst | uniq )

	#rm .tmp.movie.lst
	#rm -r ${TMP_DIR}
}; export -f collect_by_movie

merge() { # TODO: deprecated !!
# This merges alignments for monomers into one file, retaining acceptable header.
#if [ 1 == 0 ]; then
#fi
	for f in $( ls blast_output_${TRIAL} | sort ); do
		echo "processing $f"
		head -n 3 blast_output_${TRIAL}/$f >> .headers.tmp
		gawk 'NR>3' blast_output_${TRIAL}/$f >> .records.tmp
	done
	cat <(grep @HD .headers.tmp | head -n1) \
		<(grep @SQ .headers.tmp) \
		<(grep @PG .headers.tmp | head -n1) .records.tmp > merged.sam
	rm .headers.tmp .records.tmp

# Convert merged.sam to bam with MD tag. # TODO: parametrize reference monomer seq
	samtools view -@ 12 -b -F 4 merged.sam > merged.bam && \
	samtools calmd -@ 12 -b merged.bam ${MONOMER_DB} > merged.md.bam && \
	rm merged.bam
	echo "md tags calculated for merged.bam, which is gone already."

# Convert merged sam to bam sorted in reference, then index it.
	samtools sort -@ 12 -o merged.sort.bam merged.md.bam && \
	samtools index merged.sort.bam
	echo "got bam(.bai) sorted by ref (monomers)."

# Convert merged sam to bam sorted in read ID, then get sam to be parsed.
	samtools sort -@ 12 -n -o merged.sort.read.bam merged.md.bam && \
	samtools view -h merged.sort.read.bam > merged.sort.read.sam
	echo "got sam sorted by read. ready to be parsed. deleting temp files..."

# delete temp files
	rm merged.sort.read.bam
	rm merged.md.bam
	echo "done."
}

encode() {
	SAM=$1
	source ../../venv/bin/activate
	echo "-------------------- python3 ../../EncodedRead.py encode_dp --sam ${SAM} --out ${SAM%%p0.*}p0.ori.pickle"
	python3 ../../EncodedRead.py encode_dp --sam ${SAM} --out ${SAM%%p0.*}p0.ori.pickle
}; export -f encode

#prep_dbs
#samtools faidx ${MONOMER_DB}

mkdir -p queries_${TRIAL}/
mkdir -p blast_output_${TRIAL}/

#cut -f1 ${MONOMER_DB}.fai | xargs -P 12 -I % bash -c "align %"
#echo "done alignment"
#merge

#cut -f1 ${MONOMER_DB}.fai | xargs -P 12 -I % bash -c "split_by_movie %"
#echo "done split"

# Construct nice acceptable header.
TMP_DIR=$(pwd)/tmp/blast_assignment

if [[ 1==0 ]]; then
	find ${TMP_DIR}/by_monomer/ | grep header | xargs cat > ${TMP_DIR}/.headers.tmp
	cat <(grep @HD ${TMP_DIR}/.headers.tmp | head -n1) \
		<(grep @SQ ${TMP_DIR}/.headers.tmp | sort | uniq) \
		<(grep @PG ${TMP_DIR}/.headers.tmp | head -n1) > ${TMP_DIR}/header
fi

#ls ${TMP_DIR}/by_monomer/*/ | grep s1_p0 | sed -e "s/.*\///" | sort | uniq | xargs -P 24 -I % bash -c "collect_by_movie %"
#echo "done collection"

find ./blast_assignment | grep ".sam.gz" | sort | uniq | xargs -P 16 -I % bash -c "encode %"
#find ./blast_assignment | grep ".sam.gz" | sort | uniq | head -n 4 | xargs -P 2 -I % bash -c "encode %"

