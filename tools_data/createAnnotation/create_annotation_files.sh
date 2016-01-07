#!/usr/bin/env bash

while getopts "s:l:d:c:p:b:m:" optionName; do
case "$optionName" in

s) FAIFILE="$OPTARG";;
l) LENFILE="$OPTARG";;
d) DICTFILE="$OPTARG";
c) CHROFILE="$OPTARG";;
p) FILEPATH="$OPTARG";;
b) BUILD="$OPTARG";;
m) MAPFILE="$OPTARG";;

esac
done

if [[ $FAIFILE == 1 ]]; then

	#Create .fai file
	samtools faidx $FILEPATH/$BUILD/seq/$BUILD.fa

fi

if [[ $LENFILE == 1 ]]; then

	#Create .len file
	while read line
	do      
    	    if [[ $line != "" ]]; then

        	        CHR=`echo $line | awk '{print $1}'`
            	    LENGTH=`echo $line | awk '{print $2}'`
                
                	echo "$CHR	$LENGTH" >> $FILEPATH/$BUILD/seq/$BUILD.len

        	fi

	done < $FILEPATH/$BUILD/seq/$BUILD.fa.fai

fi

if [[ $DICTFILE == 1 ]]; then

	#Create .dict file
	java -jar CreateSequenceDictionary.jar R=$FILEPATH/$BUILD/seq/$BUILD.fa O=$FILEPATH/$BUILD/$BUILD.dict

fi

if [[ $CHROFILE == 1 ]]; then

	#Split by Chromosomes
	perl splitChr.pl $FILEPATH/$BUILD/seq/$BUILD.fa

fi

if [[ $MAPFILE == 1 ]]; then

	#Create genome mappability file
	gem_do_index -i $FILEPATH/$BUILD/seq/$BUILD.fa -o $FILEPATH/$BUILD/gem_index
	gem_mappability -I $FILEPATH/$BUILD/gem_index -l 50 -o $FILEPATH/$BUILD/out50m2_$BUILD.gem.mappability
	
fi







