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

local_path=/usr/bin/ChIPQC
R_PATH='Rscript --slave '

echo "$R_PATH $local_path/ChIPQC.R --args $MINVAL $MAXVAL $CHIPFILE $CONTROLFILE $OUTPUT $OUTSTAT $PDF"
$R_PATH $local_path/ChIPQC.R --args $MINVAL $MAXVAL $CHIPFILE $CONTROLFILE $OUTPUT $OUTSTAT $PDF

