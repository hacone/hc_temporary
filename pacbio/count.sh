READS_FQ=m160526_170733_42274_c101014002550000001823231410211647_s1_p0.filtered.subreads.fastq

SQUEAKR_DIR=../squeakr/
K=6
S=20

${SQUEAKR_DIR}/squeakr-count -f -k $K -s $S -t 1 -o ./ ../data/monomers/d0.fq \
&& mv d0.fq.ser d0.fq.K$K.S$S.ser

MONS_CQF=./d0.fq.K$K.S$S.ser

if [[ -e CQF_IP.K$K.S$S.dat ]]; then echo "Warning: overwriting CQF_IP.K$K.S$S.dat. Aborting."; fi

cat ${READS_FQ} | head -n 40000 | while read line; do

	# get four lines at a time
	NAME=${line#@}; echo $line > a_read.fq
	read line; READLEN=${#line}; echo $line >> a_read.fq
	for i in {1..2}; do read line; echo $line >> a_read.fq; done

	# calc CQF inner product
	${SQUEAKR_DIR}/squeakr-count -f -k $K -s $S -t 1 -o . a_read.fq  2>&1 > /dev/null
	IP_RAW=$( ${SQUEAKR_DIR}/squeakr-inner-prod -a $MONS_CQF -b a_read.fq.ser 2>&1 | grep "Inner product" )
	IP_RAW=${IP_RAW#Inner product: }
	IP_NORM=$( echo "scale=4; $IP_RAW / $READLEN" | bc -l )

	echo -e "$NAME\t$READLEN\t$IP_RAW\t$IP_NORM" >> CQF_IP.K$K.S$S.dat
done

rm a_read.fq a_read.fq.ser
