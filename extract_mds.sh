#!/bin/bash

### extract mismatch set from SAM (se or pe?) ### TODO: use python

set -o noglob

function extract_mds() {

	sam=$1

	# grep -v "^@" ${sam} > .tmp

	echo -e "QID\tFLAG\tRID\tRPOS\tMAPQ\tCIGAR\tNM\tMD\tXS\tXA"

	while read rec; do
		set $rec

		# TODO maybe comform to sam spec?
		qid=$1; flg=$2;

		if [[ $flg == "4" ]]; then continue; fi

		# assertion
		if [[ "${12:0:5}" != 'NM:i:' ]]; then
			echo "BAD LINE: $rec"
		fi

		rid=$3; rpos=$4; mapq=$5; cigar=$6; nm=${12:5}; md=${13:5};

		xs=${16}; xa=${17} # need annotation ?

		#echo -e "$rid\t$rpos\t$mapq\t$cigar\t$nm\t$md"
		#echo -e "$qid\t$flg\t$rid\t$rpos\t$mapq\t$cigar\t${nm}\t${md}\t${xs}\t${xa}"

	done < .tmp

}

extract_mds ./SRR3189743.un1.single.sam > ./.tmp.mapped2
