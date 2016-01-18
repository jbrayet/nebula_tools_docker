#!/bin/bash


while getopts "p:w:o:h:" optionName; do
case "$optionName" in

p) PEAKS_FILE="$OPTARG";;
w) WIG_FILE="$OPTARG";;
o) SUB_PEAKS_FILE="$OPTARG";;
h) GALAXY_HOME="$OPTARG";;

esac
done

ROOT_NAME=`basename $PEAKS_FILE .dat`
SUB_PEAKS_FILENAME=${ROOT_NAME}.subpeaks.dat
TMP_DIR=`mktemp -d`

/usr/bin/peaksplitter/PeakSplitter_v1/PeakSplitter -p $PEAKS_FILE -w $WIG_FILE -o $TMP_DIR -f

mv $TMP_DIR/$SUB_PEAKS_FILENAME $SUB_PEAKS_FILE

