#!/bin/bash

nameBed=" "

while getopts "f:c:t:v:m:r:o:w:b::n::q::l::" optionName; do
case "$optionName" in

f) inputfile="$OPTARG";;
c) controlfile="$OPTARG";;
t) minHeight="$OPTARG";;
v) minRatio="$OPTARG";;
m) minHeightControl="$OPTARG";;
r) ROOT_DIR="$OPTARG";;
o) output="$OPTARG";;
w) outputControl="$OPTARG";;
b) BEDFILE="$OPTARG";;
n) nameBed="$OPTARG";;
q) nameBedControl="$OPTARG";;
l) BUILD="$OPTARG";;
esac
done

local_path=/usr/bin/filterpeaks
DATABASE_PATH=$ROOT_DIR/database/files
mkdir -p $DATABASE_PATH/nebulaAnnotations
mkdir -p $DATABASE_PATH/nebulaAnnotations/$BUILD
nebulaAnnotationPath=$DATABASE_PATH/nebulaAnnotations
nebulaGenomePath=$DATABASE_PATH/nebulaAnnotations/$BUILD

#echo "$@" > /data/tmp/tmp.log
echo "java -classpath local_path/ -Xmx6g FilterPeaks -f $inputfile -c $controlfile -t $minHeight -v $minRatio -o $output"
java -classpath $local_path/ -Xmx6g FilterPeaks -f $inputfile -c $controlfile -t $minHeight -v $minRatio -o $output 
java -classpath $local_path/ -Xmx6g FilterPeaks -c $inputfile -f $controlfile -t $minHeightControl -v $minRatio -o $outputControl

if [ "$BEDFILE" == "yes" ]; then

  ############### Create annotations files ################

  FAIFILE='n'
  LENFILE='n'
  DICTFILE='n'
  CHROFILE='n'
  MAPFILE='n'

  if [ ! -f $nebulaGenomePath/$BUILD.len ]; then
      FAIFILE='y'
      LENFILE='y'
  fi
  
  chrom_info_file=$nebulaGenomePath/$BUILD.len
  FILEPATH=$ROOT_DIR/tool-data
  
  bash /usr/bin/create_annotation_files.sh $FAIFILE $LENFILE $DICTFILE $CHROFILE $FILEPATH $BUILD $MAPFILE $nebulaGenomePath
  
  #################### END ANNOTATION FILES ###########################

  perl $local_path/peak2bed.pl -f $output -t 0 -o $output.tmp -g $chrom_info_file -n $nameBed
  perl $local_path/peak2bed.pl -f $outputControl -t 0 -o $outputControl.tmp -g $chrom_info_file -n $nameBedControl 
  cp $output.tmp $output
  rm $output.tmp
  cp $outputControl.tmp $outputControl
  rm $outputControl.tmp
fi


