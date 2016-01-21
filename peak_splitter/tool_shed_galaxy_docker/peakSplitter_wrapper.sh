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

trackNumber=$(grep -c "track" $WIG_FILE)
if [[ trackNumber > 1 ]]; then
    cp $WIG_FILE $WIG_FILE.tmp
    sed -i "/track/d" $WIG_FILE.tmp
    line=$(head -n 1 $WIG_FILE)
    sed -i "1s/^/$line\n/" $WIG_FILE.tmp
    /usr/bin/peaksplitter/PeakSplitter_v1/PeakSplitter -p $PEAKS_FILE -w $WIG_FILE.tmp -o $TMP_DIR -f
    rm -f $WIG_FILE.tmp
else
    /usr/bin/peaksplitter/PeakSplitter_v1/PeakSplitter -p $PEAKS_FILE -w $WIG_FILE -o $TMP_DIR -f
fi

mv $TMP_DIR/$SUB_PEAKS_FILENAME $SUB_PEAKS_FILE
