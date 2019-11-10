#!/bin/bash

# General Settings

export HC=/work2/hacone/2018/human_centromeres

# TODO: parametrize this
#export MONOMER_DB=/data/hacone/s14-X17.fa
#export MONOMER_DB=${HC}/cluster.s14.SRR3189741.fa
export MONOMER_DB=${HC}/resource/monomers/Hum14AlpMon.fa

# export READ_DIR=${HC}/pacbio/Centro_Fasta_143Cells
# export READ_DIR=${HC}/pacbio/20190917/reads/
# export TMP_DIR=/glusterfs/hacone/blast-tmp/
export TMP_DIR=$3

# Required for `sorted`
export SAMTOOLS=samtools

# Required for `encode`
export VENV_ACT=${HC}/venv/bin/activate
export ENCODEPY=${HC}/EncodedRead.py

#export MONLIST=${HC}/resource/monomers/s14-for-X.txt
#export MONLIST=${HC}/resource/monomers/s14-for-32.txt
#export MONLIST=${HC}/resource/monomers/s14-for-chr17.txt

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

export WORD=8
export QCOV=60
export PCID=60
export GOPEN=1
export GEXT=1
export PARAM=w$WORD.p$QCOV.i$PCID.go$GOPEN.ge$GEXT

align() {

  ## This should take the full path of the reads in .fasta.gz
  READ=$1
  LOCAL_TMP=${TMP_DIR}/$( basename $READ )
  LOCAL_TMP=${LOCAL_TMP%%.fasta.gz}
  REF=${LOCAL_TMP}/$( basename $READ )
  REF=${REF%%.gz}

  if [[ -e ${LOCAL_TMP}/$( basename $MONOMER_DB ).sam ]]; then
    echo ".sam exist. go for sorting; $READ"
    ls -lah ${LOCAL_TMP}/$( basename $MONOMER_DB ).sam
    exit
  fi

  if [[ -e ${LOCAL_TMP}/$( basename $MONOMER_DB ).sort.read.sam.gz ]]; then
    echo ".sam.gz exist (already done upto sorting); $READ"
    ls -lah ${LOCAL_TMP}/$( basename $MONOMER_DB ).sort.read.sam.gz
    exit
  fi

  ## Prepare directory and database
  mkdir -p ${LOCAL_TMP}
  if [[ ! -e $REF ]]; then
          zcat ${READ} > ${REF}
          makeblastdb -dbtype nucl -in ${REF} -parse_seqids
          echo "Aa. Created reference for ${READ}"
  else
          echo "Aa. Found reference for ${READ}"
  fi


  ## Align all the monomers 
  blastn -db ${REF} -query ${MONOMER_DB} -out ${LOCAL_TMP}/$( basename $MONOMER_DB ).sam \
		-max_target_seqs 1000000 -word_size 8 -qcov_hsp_perc 60 \
		-outfmt "17 SQ" -parse_deflines

  #blastn -db ${REF} -query ${MONOMER_DB} -out ${LOCAL_TMP}/$( basename $MONOMER_DB ).$PARAM.sam \
	#	-max_target_seqs 1000000 -word_size $WORD -qcov_hsp_perc $QCOV -perc_identity $PCID \
  #  -gapopen $GOPEN -gapextend $GEXT \
	#	-outfmt "17 SQ" -parse_deflines
	#-max_target_seqs 1000000 -word_size 8 -qcov_hsp_perc 80 -perc_identity 80 \
	#-max_target_seqs 1000000 -word_size 7 -qcov_hsp_perc 60 \

  echo "Ab. Aligned for ${READ}"

}; export -f align

sorted() {

# Generate sorted, indexed, binarized, and compressed version of the specified SAM file.
# In the end, you'll have 1 sam (read-sorted gzipped) and 2 bams (ref-sorted and read-sorted).

  # NOTE: These should comply with the definition above.
  READ=$1
  LOCAL_TMP=${TMP_DIR}/$( basename $READ )
  LOCAL_TMP=${LOCAL_TMP%%.fasta.gz}

  # full path, and without the extension
  #SAM=${LOCAL_TMP}/$( basename $MONOMER_DB ).$PARAM.sam
  SAM=${LOCAL_TMP}/$( basename $MONOMER_DB ).sam
  SA=${SAM%%.sam}

  if [ ! -e $SAM ]; then
          echo "$SAM not exists"
          exit
  fi

# Convert merged.sam to temporary bam with MD tag.
		${SAMTOOLS} view -@ 12 -b -F 4 ${SAM} > ${SA}.bam && \
		${SAMTOOLS} calmd -@ 12 -b ${SA}.bam ${MONOMER_DB} > ${SA}.md.bam 2> /dev/null && \
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
    if [[ -e ${SA}.sort.read.sam.gz ]]; then
            rm ${SA}.sam ${SA}.bam
            rm ${SA}.md.bam
            rm ${SA}.sort.read.bam

            # remove for ~600 ...
            rm ${SA}.sort.bam ${SA}.sort.bam.bai
    fi
		echo "Bd. Done."

}; export -f sorted

encode() {

  # This should take the full path of the sam.gz

	SAM=$1
	source ${VENV_ACT}
  RAND=$(cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 12 | head -n 1)

  # TODO: parametrize!
  MONLIST=10mons.lst

  #${PYTHON3} ${ENCODEPY} encode_dp --sam ${SAM} --mons ${MONLIST} --out ${SAM%%.sam.gz}.10mons.$RAND.pickle
  #export MONLIST=${HC}/resource/monomers/s14-for-X.txt
  ${PYTHON3} ${ENCODEPY} encode_dp --sam ${SAM} --mons ${MONLIST} --out ${SAM%%.sam.gz}.10mons.${RAND}.pickle

  echo "Done encoding."

}; export -f encode

# NOTE: HERE IS THE ENTRY POINT OF THE SCRIPT!

CMD=$1
READ_DIR=$2

if [[ $CMD == "all" ]]; then

  echo "$CMD $READ_DIR $TMP_DIR "
  find $READ_DIR | grep .fasta.gz$ | xargs -P 12 -I % bash -c "align %"
  find $READ_DIR | grep .fasta.gz$ | xargs -P 12 -I % bash -c "sorted %"

  exit

  MONREF=Hum14AlpMon
  # TODO do this for 1192 mons later.
  # MONREF=s14-X17
  find $TMP_DIR | grep .sam.gz$ | grep ${MONREF} \
          | xargs -P 12 -I % bash -c "encode %"
  find $TMP_DIR | grep .pickle$ | grep ${MONREF} \
          > $TMP_DIR/${MONREF}.mers.fofn 
  ./hor-encode.sh $TMP_DIR/${MONREF}.mers.fofn > $TMP_DIR/${MONREF}.hers.dat
  ./hor-summary.sh $TMP_DIR/${MONREF}.hers.dat > $TMP_DIR/${MONREF}.hers.summary
fi

if [[ $CMD == "sort" ]]; then
        # TODO
  find $READ_DIR | grep .fasta.gz$ | xargs -P 12 -I % bash -c "sorted %"
fi

if [[ $CMD == "encode" ]]; then
  MONREF=Hum14AlpMon
  # TODO do this for 1192 mons later.
  # MONREF=s14-X17
  find $TMP_DIR | grep .sam.gz$ | grep ${MONREF} \
          | xargs -P 12 -I % bash -c "encode %"
fi

if [[ $CMD == "hor-encode" ]]; then
  MONREF=Hum14AlpMon
  # TODO do this for 1192 mons later.
  # MONREF=s14-X17
  #find $TMP_DIR | grep .pickle$ | grep ${MONREF} \
  #        > $TMP_DIR/${MONREF}.mers.fofn 
  #./hor-encode.sh $TMP_DIR/${MONREF}.mers.fofn > $TMP_DIR/${MONREF}.hers.dat
  ./hor-summary.sh $TMP_DIR/${MONREF}.hers.dat > $TMP_DIR/${MONREF}.hers.summary2
fi
