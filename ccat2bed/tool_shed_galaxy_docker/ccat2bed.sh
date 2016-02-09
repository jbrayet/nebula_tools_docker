#!/bin/bash

while getopts "f:t:o:g:n:r:" optionName; do
case "$optionName" in

f) inputfile="$OPTARG";;
t) minHeight="$OPTARG";;
o) output="$OPTARG";;
g) BUILD="$OPTARG";;
n) nameBed="$OPTARG";;
r) ROOT_DIR="$OPTARG";;
esac
done

local_path=/usr/bin/ccattobed
DATABASE_PATH=$ROOT_DIR/database/files
mkdir -p $DATABASE_PATH/nebulaAnnotations
mkdir -p $DATABASE_PATH/nebulaAnnotations/$BUILD
nebulaAnnotationPath=$DATABASE_PATH/nebulaAnnotations
nebulaGenomePath=$DATABASE_PATH/nebulaAnnotations/$BUILD

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

perl $local_path/ccat_int2bed.pl -f $inputfile -t $minHeight -o $output -g $chrom_info_file -n $nameBed

