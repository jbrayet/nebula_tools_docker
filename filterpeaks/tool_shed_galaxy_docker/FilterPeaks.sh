#!/bin/bash

nameBed=" "

while getopts "f:c:t:v:o:m:w:n:q:l:" optionName; do
case "$optionName" in

f) inputfile="$OPTARG";;
c) controlfile="$OPTARG";;
t) minHeight="$OPTARG";;
v) minRatio="$OPTARG";;
o) output="$OPTARG";;
m) minHeightControl="$OPTARG";;
w) outputControl="$OPTARG";;
n) nameBed="$OPTARG";;
q) nameBedControl="$OPTARG";;
l) lenfile="$OPTARG";;
esac
done

local_path=$(dirname $(readlink -f $0))

#echo "$@" > /data/tmp/tmp.log
echo "java -classpath /bioinfo/http/prod/hosted/nebula.curie.fr/galaxy-dist/tools/curie_tools/peak_calling/ -Xmx6g FilterPeaks -f $inputfile -c $controlfile -t $minHeight -v $minRatio -o $output" >> /dev/null
java -classpath $local_path/ -Xmx6g FilterPeaks -f $inputfile -c $controlfile -t $minHeight -v $minRatio -o $output > /dev/null
java -classpath $local_path/ -Xmx6g FilterPeaks -c $inputfile -f $controlfile -t $minHeightControl -v $minRatio -o $outputControl > /dev/null

if [ ! "$nameBed" == " " ]; then
  perl $local_path/peak2bed.pl -f $output -t 0 -o $output.tmp -g $lenfile -n $nameBed > /dev/null
  perl $local_path/peak2bed.pl -f $outputControl -t 0 -o $outputControl.tmp -g $lenfile -n $nameBedControl > /dev/null
  cp $output.tmp $output
  rm $output.tmp
  cp $outputControl.tmp $outputControl
  rm $outputControl.tmp
fi


