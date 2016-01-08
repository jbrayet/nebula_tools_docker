#!/usr/bin/env bash

FAIFILE=$1
LENFILE=$2
DICTFILE=$3
CHROFILE=$4
FILEPATH=$5
BUILD=$6
MAPFILE=$7
ANNOPATH=$8

ln -s $FILEPATH/$BUILD/seq/$BUILD.fa $ANNOPATH/$BUILD.fa

if [ "$FAIFILE" == 'y' ]
then

        #Create .fai file
        samtools faidx $ANNOPATH/$BUILD.fa

fi

if [ "$LENFILE" == 'y' ]
then

        #Create .len file
        while read line
        do
            if [[ $line != "" ]]; then

                        CHR=`echo $line | awk '{print $1}'`
                        LENGTH=`echo $line | awk '{print $2}'`
                        echo "$CHR      $LENGTH" >> $ANNOPATH/$BUILD.len

                fi

        done < $ANNOPATH/$BUILD.fa.fai

fi

if [ "$DICTFILE" == 'y' ]
then

        #Create .dict file
        java -jar CreateSequenceDictionary.jar R=$ANNOPATH/$BUILD.fa O=$ANNOPATH/$BUILD.dict

fi

if [ "$CHROFILE" == 'y' ]
then

        #Split by Chromosomes
        perl splitChr.pl $ANNOPATH/$BUILD.fa

fi

if [ "$MAPFILE" == 'y' ]
then

        #Create genome mappability file
        ./gem-indexer -i $ANNOPATH/$BUILD.fa -o $ANNOPATH/gem_index
        ./gem-mappability -I $ANNOPATH/gem_index -l 50 -o $ANNOPATH/out50m2_$BUILD.gem.mappability

fi
