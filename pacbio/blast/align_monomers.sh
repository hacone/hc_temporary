#!/bin/bash

export HC=/work2/hacone/2018/human_centromeres

# General Settings
#export MONOMER_DB=${HC}/cluster.s14.SRR3189741.fa
export MONOMER_DB=${HC}/resource/monomers/d0.fasta
export READ_DIR=${HC}/pacbio/Centro_Fasta_143Cells
export TMP_DIR=${HC}/pacbio/blast/tmp/

# Required for `sorted`
export SAMTOOLS=samtools

# Required for `encode`
export VENV_ACT=${HC}/venv/bin/activate
export ENCODEPY=${HC}/EncodedRead.py
export MONLIST=${HC}/resource/monomers/s14-for-X.txt

# in case it can't reach the correct one
export PYTHON3=${HC}/venv/bin/python3

## This performs monomer-assignment by blast. Output should be organized as sam/bam format file that can be directly parsed by python.
## As it considers monomers to be query, parsing part should be somewhat different. This specification might be changed later.

# Testing against the following version available in hx.
#
# makeblastdb: 2.4.0+
#  Package: blast 2.4.0, build Jun  1 2016 13:39:30
# blastn: 2.4.0+
#  Package: blast 2.4.0, build Jun  1 2016 13:39:30

align() {

  
  ## This should take the full path of the reads in .fasta.gz

  READ=$1
  LOCAL_TMP=${TMP_DIR}/$( basename $READ )
  REF=${LOCAL_TMP}/reads.fasta

  # TODO: remove this later
  if [[ -e ${LOCAL_TMP}/$( basename $MONOMER_DB ).sam ]]; then
    echo "exist"
    ls -lah ${LOCAL_TMP}/$( basename $MONOMER_DB ).sam
    exit
  else
    echo "not exist; go"
  fi

  ## Prepare directory and database
  mkdir -p ${LOCAL_TMP}
	zcat ${READ} > ${REF}
	makeblastdb -dbtype nucl -in ${REF} -parse_seqids

  echo "Aa. Created reference for ${READ}"

  ## Align all the monomers 
  blastn -db ${REF} -query ${MONOMER_DB} -out ${LOCAL_TMP}/$( basename $MONOMER_DB ).sam \
		-max_target_seqs 1000000 -word_size 7 -qcov_hsp_perc 60 \
		-outfmt "17 SQ" -parse_deflines

  echo "Ab. Aligned for ${READ}"

  # ${SAMTOOLS} faidx ${MONOMER_DB}
  # (temporarily abandoned) ${SAMTOOLS} faidx ${MONOMER_DB} ${MON_NAME} > queries_${TRIAL}/${MON_NAME}.fa

}; export -f align

sorted() {

# Generate sorted, indexed, binarized, and compressed version of the specified SAM file.
# In the end, you'll have 1 sam (read-sorted gzipped) and 2 bams (ref-sorted and read-sorted).

  # NOTE: These should comply with the definition above.
  READ=$1
  LOCAL_TMP=${TMP_DIR}/$( basename $READ )
  # full path, and without the extension
  SAM=${LOCAL_TMP}/$( basename $MONOMER_DB ).sam
  SA=${SAM%%.sam}

    # TODO: temporary
		gzip ${SA}.sort.read.sam
		echo "Bd. Done."
    exit

  if [ ! -e $SAM ]; then
          echo "$SAM not exists"
          exit
  fi

# Convert merged.sam to temporary bam with MD tag.
		${SAMTOOLS} view -@ 12 -b -F 4 ${SAM} > ${SA}.bam && \
		${SAMTOOLS} calmd -@ 12 -b ${SA}.bam ${MONOMER_DB} > ${SA}.md.bam 2> /dev/null && \
		rm ${SA}.bam
		echo "Ba. MD tags calculated for ${SAM} as *.md.bam"

# Convert merged sam to bam sorted in reference, then index it.
		${SAMTOOLS} sort -@ 12 -o ${SA}.sort.bam ${SA}.md.bam && \
		${SAMTOOLS} index ${SA}.sort.bam
		echo "Bb. Got bam(.bai) sorted by ref (monomers)."

# Convert merged sam to bam sorted in read ID, then get sam to be parsed.
		${SAMTOOLS} sort -@ 12 -n -o ${SA}.sort.read.bam ${SA}.md.bam && \
		${SAMTOOLS} view -h ${SA}.sort.read.bam > ${SA}.sort.read.sam
		gzip ${SA}.sort.read.sam
		echo "Bc. Got sam sorted by reads. ready to be parsed. deleting temp files..."

# delete temp files
		rm ${SA}.sort.read.bam ${SA}.md.bam
		echo "Bd. Done."

}; export -f sorted

encode() {

  # This should take the full path of the sam.gz

	SAM=$1
	source ${VENV_ACT}

  # echo debugging...
	echo "python3 ${ENCODEPY} encode_dp --sam ${SAM} --mons ${MONLIST} --out ${SAM%%.sam.gz}.pickle"
  ${PYTHON3} ${ENCODEPY} encode_dp --sam ${SAM} --mons ${MONLIST} --out ${SAM%%.sam.gz}.pickle

  echo "Done encoding."

}; export -f encode

# NOTE: HERE IS THE ENTRY POINT OF THE SCRIPT!


#find $READ_DIR | grep .fasta.gz$ | xargs -P 24 -I % bash -c "align %"

#find $READ_DIR | grep .fasta.gz$ | xargs -P 12 -I % bash -c "sorted %"

find $TMP_DIR | grep .sam.gz$ | xargs -P 16 -I % bash -c "encode %"

