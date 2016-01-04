#!/bin/bash

REG="NA"
CONTROLFILE="NA"

while getopts "f:c:l:o:r:u:v:e:p:" optionName; do
case "$optionName" in

f) CHIPFILE="$OPTARG";;
c) CONTROLFILE="$OPTARG";;
l) STEP="$OPTARG";;
r) RIGHT="$OPTARG";;
o) OUTPUT="$OPTARG";;
u) OUTSTAT="$OPTARG";;
v) GENOME="$OPTARG";;
e) REG="$OPTARG";;
p) PDF="$OPTARG";;

esac
done

LOGTMP=$OUTPUT.log.tmp

#LOGTMP=/data/tmp/log.tmp

#echo $LOGTMP
#echo "$@"

echo "$@" >> $LOGTMP

LOCAL_DIR=`( cd -P $(dirname $0); pwd)`
DOCKER_PATH='/usr/bin/maketssdist_histones'


OUTPUT_DIR=`dirname $OUTPUT`
R_PATH='/usr/bin/R-3.1.0/bin/Rscript --slave '

databasePath=$(find / -type d -name files | grep database)

mkdir -p $databasePath/nebulaAnnotations
mkdir -p $databasePath/nebulaAnnotations/noIdenticalTransc
nebulaAnnotationPath=$databasePath/nebulaAnnotations
[ ! -f $nebulaAnnotationPath/noIdenticalTransc/$GENOME.noIdenticalTransc.txt ] && curl -s -o $nebulaAnnotationPath/noIdenticalTransc/$GENOME.noIdenticalTransc.txt -L https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/noIdenticalTransc/$GENOME.noIdenticalTransc.txt

noIdenticalTranscPath=$nebulaAnnotationPath/noIdenticalTransc

if [ -r $REG ]; then
  echo "1:  perl $DOCKER_PATH/crossBedWithGenes_hist.pl -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -f $CHIPFILE -reg $REG -o $OUTPUT.annotated" >> $LOGTMP
  perl $DOCKER_PATH/crossBedWithGenes_hist.pl -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -f $CHIPFILE -reg $REG -o $OUTPUT.annotated >> $LOGTMP 2>> $LOGTMP
 else
  echo "2:  perl $DOCKER_PATH/crossBedWithGenes_hist.pl -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -f $CHIPFILE -reg $REG -o $OUTPUT.annotated" >> $LOGTMP
  perl $DOCKER_PATH/crossBedWithGenes_hist.pl -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -f $CHIPFILE -o $OUTPUT.annotated >> $LOGTMP 2>> $LOGTMP
fi

if [ -r $CONTROLFILE ]; then
  #create a subset of control peaks (highest peaks, the same number as in the sample) 
   echo "   perl $DOCKER_PATH/createControlPeakSubSet.pl -f $CHIPFILE -c $CONTROLFILE -o ${OUTPUT_DIR}/control.tmp" >>$LOGTMP
   perl $DOCKER_PATH/createControlPeakSubSet.pl -f $CHIPFILE -c $CONTROLFILE -o ${OUTPUT_DIR}/control.tmp >>$LOGTMP 2>>$LOGTMP
  if [ -r $REG ]; then
    perl $DOCKER_PATH/crossBedWithGenes_hist.pl -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -f ${OUTPUT_DIR}/control.tmp -reg $REG -o $OUTPUT.control.annotated >> $LOGTMP 2>> $LOGTMP
  else
    perl $DOCKER_PATH/crossBedWithGenes_hist.pl -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -f ${OUTPUT_DIR}/control.tmp -o $OUTPUT.control.annotated >> $LOGTMP 2>> $LOGTMP
  fi
  echo "3: cat $DOCKER_PATH/makeTSSdist_hist.R | R --slave --args $STEP $RIGHT $OUTPUT.annotated $OUTPUT $OUTSTAT $OUTPUT.control.annotated $PDF" >> $LOGTMP
  cat $DOCKER_PATH/makeTSSdist_hist.R | R --slave --args $STEP $RIGHT $OUTPUT.annotated $OUTPUT $OUTSTAT $OUTPUT.control.annotated $PDF >>$LOGTMP 2>>$LOGTMP
else
  echo "4:  cat $DOCKER_PATH/makeTSSdist_hist.R | R --slave --args $STEP $RIGHT $OUTPUT.annotated $OUTPUT $OUTSTAT $PDF 2>>$LOGTMP" >>$LOGTMP
  cat $DOCKER_PATH/makeTSSdist_hist.R | R --slave --args $STEP $RIGHT $OUTPUT.annotated $OUTPUT $OUTSTAT $PDF 2>>$LOGTMP >>$LOGTMP
fi

if [ -r $LOGTMP ]; then
  rm $LOGTMP
fi

rm $OUTPUT.annotated*

if [ -r $CONTROLFILE ]; then
  rm $OUTPUT.control.annotated*
  rm ${OUTPUT_DIR}/control.tmp
fi

