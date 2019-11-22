
export HC=/work2/hacone/2018/human_centromeres
export VENV_ACT=${HC}/venv/bin/activate
export PYTHON3=${HC}/venv/bin/python3
export ENCODEPY=${HC}/EncodedRead.py

export UAPY=${HC}/UnitAnalysis.py

function monomer_encode() {

  hor=$1
  sample=$2
  SAM=/data/hacone/filtered_alignments/${sample}.58mons.${hor}.sam.gz
  MONLIST=Monomers/${hor}ers.lst
  FOFN=${SAM%%.sam.gz}.pkl.tmp.fofn
  HER_DAT=${SAM%%.sam.gz}.her.dat
  DEF=${3:-Monomers/HOR-pentamers-chromspec.def}

  # TODO: usually I'd skip this. move the line.
  #echo "m-encoding for $hor : $sample"
  #${PYTHON3} ${ENCODEPY} encode_dp --sam ${SAM} --mons ${MONLIST} --out ${SAM%%.sam.gz}.pickle

  echo ${SAM%%.sam.gz}.pickle > ${FOFN}

  echo "h-encoding for $hor : $sample"
  ./hor-encode.sh ${FOFN} ${DEF} 
  #> ${HER_DAT}
  echo "done hor-encode for $hor : $sample"
  #./hor-summary.sh ${HER_DAT} ${DEF} > ${HER_DAT%%.dat}.summary
  #echo "done hor-summary for $hor : $sample"

  rm ${FOFN}

}; export -f monomer_encode;

function print_var() {

  hor=$1 ; sample=$2
  hortype=$3 ; outsuf=$4
  SAM=/data/hacone/filtered_alignments/${sample}.58mons.${hor}.sam.gz
  HERS=${SAM%%.sam.gz}.hor.pickle

  if [[ -e $HERS ]]; then
          echo "print-var for $hortype among $hor in $sample"
          echo -e "# Variants on:\t${hortype}\t${hor}\t${sample}" > ${HERS%%.hor.pickle}.${outsuf}.vars
          ${PYTHON3} ${UAPY} print-var --hor-reads ${HERS} --hor-type "${hortype}" --all >> ${HERS%%.hor.pickle}.${outsuf}.vars
  else
          echo "$HERS not found."
  fi

}; export -f print_var;


if [[ 1 == 0 ]]; then echo "" ; fi

#find $TMP_DIR | grep .sam.gz$ | grep ${MONREF} \
#        | xargs -P 16 -I % bash -c "encode % ${MONLIST}"

rm all.12m.X-12mW.vars 2> /dev/null

for sample in Ashkenazi CHM13Hifi Dai Esan Finnish Gujarati HG005 Maasai Mende NA12878Hifi Peruvian PuertoRican Toscani CHM13CLR; do
#for sample in Ashkenazi; do
        for hor in 12m; do
        #for hor in 5m 11m 12m 16m; do
                #monomer_encode ${hor} ${sample} "Monomers/HOR-${hor}.def" &
                #print_var ${hor} ${sample} "X-12mW" "X-12mW" &
                grep -v "^$" filtered_alignments/${sample}.58mons.12m.X-12mW.vars \
                | gawk 'BEGIN{OFS="\t"} NR>3{print "'$sample'",$0}'
        done 
        #sleep 10
done | LC_ALL=C sort -k2,2n -k3,3n -k4,4 | column -t > all.12m.X-12mW.vars
