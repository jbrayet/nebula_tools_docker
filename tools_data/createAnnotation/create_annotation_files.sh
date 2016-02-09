#!/usr/bin/env bash

FAIFILE=$1
LENFILE=$2
DICTFILE=$3
CHROFILE=$4
FILEPATH=$5
BUILD=$6
MAPFILE=$7
ANNOPATH=$8

ln -sf $FILEPATH/$BUILD/seq/$BUILD.fa $ANNOPATH/$BUILD.fa

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
                        echo "$CHR\t$LENGTH" >> $ANNOPATH/$BUILD.len

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
        perl /usr/bin/splitChr.pl $BUILD.fa $ANNOPATH $ANNOPATH/chromosomes
        
fi

if [ "$MAPFILE" == 'y' ]
then
	if [ "$BUILD" == 'mm9' ] || [ "$BUILD" == 'mm10' ] || [ "$BUILD" == 'hg18' ] || [ "$BUILD" == 'hg19' ]
	then
		wget http://zerkalo.curie.fr:8080/partage/nebulaAnnotation/mappability/out50m2_$BUILD.gem.mappability -P $ANNOPATH -q
	else
		#Create genome mappability file
        PATH=$PATH:/usr/bin
        wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/createAnnotation/gem-indexer -P /usr/bin
        chmod +x gem-indexer
		wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/createAnnotation/gem-mappability -P /usr/bin
        chmod +x gem-mappability
		wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/createAnnotation/gem-indexer_fasta2meta%2Bcont -P /usr/bin
        chmod +x gem-indexer_fasta2meta+cont
        wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/createAnnotation/gem-indexer_bwt-dna -P /usr/bin
        chmod +x gem-indexer_bwt-dna
		wget https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/createAnnotation/gem-indexer_generate -P /usr/bin
        chmod +x gem-indexer_generate
        ./gem-indexer -i $ANNOPATH/$BUILD.fa -o $ANNOPATH/gem_index.gem
        ./gem-mappability -I $ANNOPATH/gem_index.gem -l 50 -o $ANNOPATH/out50m2_$BUILD.gem
	fi
fi
