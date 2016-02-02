#!/bin/bash


while getopts "f:c:m:o:s:u:p:" optionName; do
case "$optionName" in

f) CHIPFILE="$OPTARG";;
c) CONTROLFILE="$OPTARG";;
m) MINVAL="$OPTARG";;
s) MAXVAL="$OPTARG";;
o) OUTPUT="$OPTARG";;
u) OUTSTAT="$OPTARG";;
p) PDF="$OPTARG";;

esac
done

LOGTMP=$OUTPUT.log.tmp

local_path=/usr/bin/ChIPQC
R_PATH='Rscript --slave '

#LOGTMP=/data/tmp/log.tmp

echo "$@" > $LOGTMP
echo "$R_PATH $local_path/ChIPQC.R --args $MINVAL $MAXVAL $CHIPFILE $CONTROLFILE $OUTPUT $OUTSTAT $PDF" >>$LOGTMP
$R_PATH $local_path/ChIPQC.R --args $MINVAL $MAXVAL $CHIPFILE $CONTROLFILE $OUTPUT $OUTSTAT $PDF >>$LOGTMP
if [ -r $LOGTMP ]; then
  rm $LOGTMP
fi
