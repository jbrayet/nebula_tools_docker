#:t:::::::::::::::::g@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#:t::::::::::::::;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#:::::::::::::z;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#::::::::::::i@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#::::::::::::@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@$@@@@
#:::::::::::3@@@@@@@@@@@@@@@@@@@@@@@@@B@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#::::::::::3@@@@@@@@@@@@@@@@@@@@@BEEESSE5EEEEBBM@@@@@@@@@@@@@@@@@@@@@@@@@@
#::::::::::3@@@@@@@@@@@@@@@@@@@@BEEEEEE35EE55E2355E5SBMB@@@@@@@@@@@@@@@@@$
#::::::::::@@@@@@@@@@@@@@@@@@@EEEE55533t3tttt::::::!!!!7755E755SBBMMM@@@MM
#::::::::::3@@@@@@@@@@@@@@@@@@EEEE2t3ttttt:::::::::::::::::::::::!7?5225EE
#::::::::::3@@@@@@@@@@@@@@@@@@EEEEE31t::::::::::::::::::::::::::::::::3E5@
#::::::::::3@@@@@@@@@@@@@@@@@@EEEEEEtt:::::::::::::::::::::::::::::::::353
#::::::::::3@@@@@@@@@@@@@@@@@@EEEEEE1ttz::::::::::::::::::::::::::::::::35
#:::::::::::@@@@@@@@@@@@@@@@@@EEEEEEEtz1::::::::::::::::::::::::::::::::t:
#:::::::::!3@@@@@@@@@@@@@@@@@@@EEEEEttt::::::::::::::::::::::::::::::::;zz
#::::::::::@@@@@@@@@@@@@@@@@@@@EEEEEttt:::::z;z:::::::::::::::::::::::::13
#::::::::::3B@@@@@@@@@@@@@@@@@@EEEEEEE3tt:czzztti;:::::::::::::::::::::::3
#::::ttt::::3@@@@@@@@@@@@@@@@EEEEE5EE25Ezt1EEEz5Etzzz;;;;:::::::::::::::::
#:::::::::::I9@@@@@@@@@@@@@@@@@@@@@@@@@@EEEEEE@@@@@@@@@@@@@@Ez;:::::::::::
#:::::::::::::E@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@Ez::::::
#::::::::::::::E@@@@@@@@@@@@@@@@@@@@@@@@@@@@@BE5EBB@@@@@@@@@@@@@@@EEE:::::
#:::::::::::::::@@@@@@@@@@@@@@@@@@@@@@@@@@@@E1::35@@@@@@@@@@ME3MMME2::::::
#:::::::::::::::?@@@@@@@@@@@@@@@@@@M@@@@@@@EE:::::3SB@@BBESEEt::::::::::::
#::::::::::::::::J$@@@@@@@B@@@@@@@@@@@@@@@@EE:::::::!35E33t:::::::::::::::
#:::::::::::::::::3@E@@@EE5EESE5EESE@@@@@@@Et::::::::::::tz:::::::::::::::
#:::::::::::::::::J@E$@EEE5133555SE@@@@@@@@Et:::::::::::::::::::::::::::::
#::::::::::::::::::E@E@EEEEtt3523EEE@@@@@@@E::::::::::::::::::::::::::::::
#:t::::::::::::::::JEE3@@@EEEEEEEEEE@@@@@@@E:::::::::t;:::::::::::::::::::
#:t:::::::::::::::::!5ES@EEEEEEEEES@@@@@@@@@E;:::;;;:3Ez::::::::::::::::::
#:t::::::::::::::::::::JE@@EEEEEEE@@@@@@@@@@@@@@@@ME!:::;:::::::::::::::::
#:tz::::::::::::::::::::JE@@@EEEE@@@@@@@@@@@@@@EE!:::::::t::::::::::::::::
#:t::::::::::::::::::::::3@@@@@@@@@@@@@@@@@@ESBE::::::::::::::::::::::::::
#:::::::::::::::::::::::::Q@@@@@@@@@@@@@@@@EE3EE;:::::zzzz::::::::::::::::
#:::::::::::::::::::::::::3@@@@@@@@@@@@@@@@@@@@@@NN@@@@@@Ez:::::::::::::::
#:zt:::::::::::::::::::::::3@@@@EE@@@@@@@@@@EEEEt::;z113E5t:::::::::::::::
#::tt:::::::::::::::::::::::3@@@E@@@@@@@@@@@@@@@@BEt::::::::::::::::t:::::
#:tt:t:::::::::::::::::::::::?S@@@@@@@@@@@BBEEE51!::::::::::::::zzzEt:::::
#::::::::::::::::::::::::::::::3Q@@@@@@@BEEEEEt:::::::::::::;zz@@@EE::::::
#::::::::::::::::::::::::::::::::75B@@@@@EEEtt;:::::::::;zz@@@@BEEEtz:::::
#::::::::::::::::::::::::::::::::::::?9@@@@@@@@@@@E2Ezg@@@@@B@@@EEEE1t::::
#:::::::::::::::::::::::::::::::::::::::3@@@@@@@@@@@@@@@@@@@E@EEEEEEEzzz::
#::::::::::::::::::::::::::::::::::::;@@@@@@@@@@@@@@@@@@@@@@@EEEEEEE5ttttt
#:::::::::::::::::::::::::::::::;g@@@@@@@@@@@@@@@@@@@@@@@@@@EEEEEEEEEEEtzt
#::::::::::::::::::::::::::::;@@@@@@@@@@@@@@@@@@@@@@@@@@E@@EEEEEEEEEEEE@@@
#::::::::::::::::::::::::::g@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@EEEE3EEEE@@@@@@@
#:::::::::::::::::::::;;g@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@EEEt33@@@@@@@@@@
#:::::::::::::::::;g@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@E@@@@@@EEEtg@@@@@@@@@@@@
#::::::::::::::;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@EEEE@@@@@@@@@@@@@@@@@@@@@@@@
#:::::::::::::@@@@@@@@@@@@@@@@@$@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#::::::::::;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
# Copyleft ↄ⃝ 2012 Institut Curie
# Author(s): Valentina Boeva, Alban Lermine (Institut Curie) 2012
# Contact: valentina.boeva@curie.fr, alban.lermine@curie.fr
# This software is distributed under the terms of the GNU General
# Public License, either Version 2, June 1991 or Version 3, June 2007.

#!/bin/bash

REG="NA"
CONTROLFILE="NA"

while getopts "f:c:l:o:r:u:v:e:p:" optionName; do
case "$optionName" in

f) CHIPFILE="$OPTARG";;
c) CONTROLFILE="$OPTARG";;
l) STEP="$OPTARG";;
r) RIGHT="$OPTARG";;
o) OUTPUT="$OPTARG";;
u) OUTSTAT="$OPTARG";;
v) GENOME="$OPTARG";;
e) REG="$OPTARG";;
p) PDF="$OPTARG";;

esac
done

LOCAL_DIR=`( cd -P $(dirname $0); pwd)`
DOCKER_PATH='/usr/bin/maketssdist'


OUTPUT_DIR=`dirname $OUTPUT`
R_PATH='/usr/bin/R-3.1.0/bin/Rscript --slave '

databasePath=$(find / -type d -name files | grep database)

mkdir -p $databasePath/nebulaAnnotations
mkdir -p $databasePath/nebulaAnnotations/noIdenticalTransc
nebulaAnnotationPath=$databasePath/nebulaAnnotations
[ ! -f $nebulaAnnotationPath/noIdenticalTransc/$GENOME.noIdenticalTransc.txt ] && curl -s -o $nebulaAnnotationPath/noIdenticalTransc/$GENOME.noIdenticalTransc.txt -L https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/noIdenticalTransc/$GENOME.noIdenticalTransc.txt

noIdenticalTranscPath=$nebulaAnnotationPath/noIdenticalTransc


if [ -r $REG ]; then
  perl $DOCKER_PATH/crossBedWithGenes.pl -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -f $CHIPFILE -reg $REG -o $OUTPUT.annotated
 else
  perl $DOCKER_PATH/crossBedWithGenes.pl -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -f $CHIPFILE -o $OUTPUT.annotated
fi

if [ -r $CONTROLFILE ]; then
   perl $DOCKER_PATH/createControlPeakSubSet.pl -f $CHIPFILE -c $CONTROLFILE -o ${OUTPUT_DIR}/control.tmp
  if [ -r $REG ]; then
    perl $DOCKER_PATH/crossBedWithGenes.pl -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -f ${OUTPUT_DIR}/control.tmp -reg $REG -o $OUTPUT.control.annotated
  else
    perl $DOCKER_PATH/crossBedWithGenes.pl -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -f ${OUTPUT_DIR}/control.tmp -o $OUTPUT.control.annotated
  fi
    $R_PATH $DOCKER_PATH/makeTSSdist.R --args $STEP $RIGHT $OUTPUT.annotated $OUTPUT $OUTSTAT $OUTPUT.control.annotated $PDF
else
  $R_PATH $DOCKER_PATH/makeTSSdist.R --args $STEP $RIGHT $OUTPUT.annotated $OUTPUT $OUTSTAT $PDF
fi

rm $OUTPUT.annotated*

if [ -r $CONTROLFILE ]; then
  rm $OUTPUT.control.annotated*
  rm ${OUTPUT_DIR}/control.tmp
fi
