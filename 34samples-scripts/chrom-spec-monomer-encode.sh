
export HC=/work2/hacone/2018/human_centromeres
export VENV_ACT=${HC}/venv/bin/activate
export PYTHON3=${HC}/venv/bin/python3
export ENCODEPY=${HC}/EncodedRead.py

export UAPY=${HC}/UnitAnalysis.py

function monomer_encode() {

  hor=$1 ; sample=$2
  SAM=/data/hacone/filtered_alignments/${sample}.58mons.${hor}.sam.gz
  MONLIST=Monomers/${hor}ers.lst
  FOFN=${SAM%%.sam.gz}.pkl.tmp.fofn

  echo "m-encoding for $hor : $sample"
  ${PYTHON3} ${ENCODEPY} encode_dp --sam ${SAM} --mons ${MONLIST} --out ${SAM%%.sam.gz}.pickle
  echo ${SAM%%.sam.gz}.pickle > ${FOFN}

}; export -f monomer_encode;

function hor_encode() {

  hor=$1 ; sample=$2

  SAM=/data/hacone/filtered_alignments/${sample}.58mons.${hor}.sam.gz
  FOFN=${SAM%%.sam.gz}.pkl.tmp.fofn

  HER_DAT=${SAM%%.sam.gz}.her.dat
  DEF=${3:-Monomers/HOR-pentamers-chromspec.def}

  echo ${SAM%%.sam.gz}.pickle > ${FOFN}

  echo "h-encoding for $hor : $sample"

  if [[ $4 == "print" ]]; then
          ./hor-encode.sh ${FOFN} > ${HER_DAT}
  else
          ./hor-encode.sh ${FOFN} ${DEF} > ${HER_DAT}
  fi
  echo "done h-encoding for $hor : $sample"

  ./hor-summary.sh ${HER_DAT} ${DEF} | column -t > ${HER_DAT%%.dat}.summary
  echo "done hor-summary for $hor : $sample"

}; export -f hor_encode;

function print_var() {

  hor=$1 ; sample=$2
  hortype=$3 ; outsuf=$4
  SAM=/data/hacone/filtered_alignments/${sample}.58mons.${hor}.sam.gz
  HERS=${SAM%%.sam.gz}.hor.pickle

  if [[ -e $HERS ]]; then
          echo "print-var for $hortype among $hor in $sample"
          echo -e "# Variants on:\t${hortype}\t${hor}\t${sample}" > ${HERS%%.hor.pickle}.${outsuf}.vars
          ${PYTHON3} ${UAPY} print-var --hor-reads ${HERS} \
                  --hor-type "${hortype}" --all >> ${HERS%%.hor.pickle}.${outsuf}.vars
  else
          echo "$HERS not found."
  fi

}; export -f print_var;


if [[ 1 == 0 ]]; then
        echo "" ;

#find $TMP_DIR | grep .sam.gz$ | grep ${MONREF} \
#        | xargs -P 16 -I % bash -c "encode % ${MONLIST}"

fi # NOTE: entry

# TODO: next 5m, 11m, 16m for all samples
#for sample in Ashkenazi; do
#for hor in 12m 16m 11m 5m; do
for hor in 12m 16m; do
for sample in Esan CHM13CLR Ashkenazi CHM13Hifi Dai Finnish Gujarati HG005 Maasai Mende NA12878Hifi Peruvian PuertoRican Toscani; do
#        for hor in 5m 11m 12m 16m; do
                hor_encode ${hor} ${sample} "Monomers/HOR-${hor}.def" "print" &
                # "print"
done
sleep 30
done

exit

for sample in Ashkenazi CHM13Hifi Dai Esan Finnish Gujarati HG005 Maasai Mende NA12878Hifi Peruvian PuertoRican Toscani CHM13CLR; do
        print_var 12m ${sample} "X-12mW" "12mW" &
        print_var 5m ${sample} "11-5mW" "5mW" &
        print_var 11m ${sample} "1-11mW" "11mW" &
        print_var 16m ${sample} "17-16mW" "16mW" &
        sleep 30
done 

mkdir -p HOR_SNVs/
for hor in 5m 11m 12m 16m; do
        for sample in Ashkenazi CHM13Hifi Dai Esan Finnish Gujarati HG005 Maasai Mende NA12878Hifi Peruvian PuertoRican Toscani CHM13CLR; do
                grep -v "^$" filtered_alignments/${sample}.58mons.${hor}.${hor}W.vars \
                | gawk 'BEGIN{OFS="\t"} NR>3{print "'$sample'",$0}'
        done | LC_ALL=C sort -k2,2n -k3,3n -k4,4 | column -t > HOR_SNVs/all.${hor}.${hor}W.vars
done 

cd HOR_SNVs/
${PYTHON3} variants.py

for hor in 5m 11m 12m 16m; do

        mv ${hor}.vars.tab ${hor}.vars.tab.tmp

        cat <( head -n 1 ${hor}.vars.tab.tmp ) \
            <( sort -k39,39nr ${hor}.vars.tab.tmp | sed -e "s/\t0\.00/\t-/g" ) \
        | column -t > ${hor}.vars.tabp.p0

        cat <( head -n 1 ${hor}.vars.tab.tmp ) \
            <( sort -k39,39nr ${hor}.vars.tab.tmp ) \
        | column -t > ${hor}.vars.tab.p0

        rm ${hor}.vars.tab.tmp
done

exit # NOTE: exit
