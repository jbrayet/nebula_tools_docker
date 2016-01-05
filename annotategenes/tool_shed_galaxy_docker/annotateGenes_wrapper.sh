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
BOOTSTRAP=1

while getopts "f:c:l:o:r:u:v:e:h:d:x:y:p:b:" optionName; do
case "$optionName" in

f) CHIPFILE="$OPTARG";;
c) CONTROLFILE="$OPTARG";;
l) LEFTPROM="$OPTARG";;
r) RIGHTPROM="$OPTARG";;
o) OUTPUT="$OPTARG";;
u) OUTSTAT="$OPTARG";;
v) GENOME="$OPTARG";;
e) REG="$OPTARG";;
h) ENH="$OPTARG";;
d) DOWNGENE="$OPTARG";;
x) CONTROLSTAT="$OPTARG";;
y) LOG="$OPTARG";;
p) PDF="$OPTARG";;
b) BOOTSTRAP="$OPTARG";;

esac
done

LOGTMP=$OUTPUT.log.tmp

echo "$@" >> $LOGTMP

LOCAL_DIR=$(dirname $(readlink -f $0))
DOCKER_PATH='/usr/bin/annotateGenes'

OUTPUT_DIR=`dirname $OUTPUT`
R_PATH='Rscript --slave '

databasePath=$(find / -type d -name files | grep database)

mkdir -p $databasePath/nebulaAnnotations
mkdir -p $databasePath/nebulaAnnotations/noIdenticalTransc
nebulaAnnotationPath=$databasePath/nebulaAnnotations
[ ! -f $nebulaAnnotationPath/noIdenticalTransc/$GENOME.noIdenticalTransc.txt ] && curl -s -o $nebulaAnnotationPath/noIdenticalTransc/$GENOME.noIdenticalTransc.txt -L https://github.com/jbrayet/nebula_tools_docker/raw/master/tools_data/noIdenticalTransc/$GENOME.noIdenticalTransc.txt

noIdenticalTranscPath=$nebulaAnnotationPath/noIdenticalTransc

echo "ChIP:" >$LOG

if [ -r $REG ]; then
  echo "1:  perl $DOCKER_PATH/geneAnnotation.pl -HMmode 0 -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -tf $CHIPFILE -selG $REG -o $OUTSTAT  -lp $LEFTPROM -rightp $RIGHTPROM -enh $ENH -dg $DOWNGENE" >> $LOGTMP
  perl $DOCKER_PATH/geneAnnotation.pl -HMmode 0 -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -tf $CHIPFILE -selG $REG -o $OUTSTAT  -lp $LEFTPROM -rightp $RIGHTPROM -enh $ENH -dg $DOWNGENE  >>$LOG 2>> $LOGTMP
 else
  echo "2: perl $DOCKER_PATH/geneAnnotation.pl -HMmode 0 -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -tf $CHIPFILE -o $OUTSTAT  -lp $LEFTPROM -rightp $RIGHTPROM -enh $ENH -dg $DOWNGENE " >> $LOGTMP
  perl $DOCKER_PATH/geneAnnotation.pl -HMmode 0 -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -tf $CHIPFILE -o $OUTSTAT  -lp $LEFTPROM -rightp $RIGHTPROM -enh $ENH -dg $DOWNGENE  >>$LOG 2>> $LOGTMP
fi



if [ -r $CONTROLFILE ]; then
echo "" >>$LOG
echo "Control:" >>$LOG
  #create a subset of control peaks (highest peaks, the same number as in the sample) 
   echo "   perl $DOCKER_PATH/createControlPeakSubSet.pl -f $CHIPFILE -c $CONTROLFILE -o $CONTROLFILE.tmp -n $BOOTSTRAP " >>$LOGTMP
   perl $DOCKER_PATH/createControlPeakSubSet.pl -f $CHIPFILE -c $CONTROLFILE -o $CONTROLFILE.tmp -n $BOOTSTRAP >>$LOGTMP 2>>$LOGTMP
  if [ -r $REG ]; then
    echo "5: perl $DOCKER_PATH/geneAnnotation.pl -HMmode 0 -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -tf $CONTROLFILE.tmp -selG $REG -o $CONTROLSTAT  -lp $LEFTPROM -rightp $RIGHTPROM -enh $ENH -dg $DOWNGENE ">> $LOGTMP
    perl $DOCKER_PATH/geneAnnotation.pl -HMmode 0 -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -tf $CONTROLFILE.tmp -selG $REG -o $CONTROLSTAT  -lp $LEFTPROM -rightp $RIGHTPROM -enh $ENH -dg $DOWNGENE  >>$LOG 2>> $LOGTMP
	for (( i=2; i<=$BOOTSTRAP; i++ ))
	do
           perl $DOCKER_PATH/geneAnnotation.pl -HMmode 0 -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -tf $CONTROLFILE.tmp$i -selG $REG -o $CONTROLSTAT$i  -lp $LEFTPROM -rightp $RIGHTPROM -enh $ENH -dg $DOWNGENE  >>$LOG 2>> $LOGTMP
	done
  else
    perl $DOCKER_PATH/geneAnnotation.pl -HMmode 0 -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -tf $CONTROLFILE.tmp -o $CONTROLSTAT  -lp $LEFTPROM -rightp $RIGHTPROM -enh $ENH -dg $DOWNGENE  >>$LOG 2>> $LOGTMP
    for (( i=2; i<=$BOOTSTRAP; i++ ))
     do
        perl $DOCKER_PATH/geneAnnotation.pl -HMmode 0 -g $noIdenticalTranscPath/$GENOME.noIdenticalTransc.txt -tf $CONTROLFILE.tmp$i -o $CONTROLSTAT$i  -lp $LEFTPROM -rightp $RIGHTPROM -enh $ENH -dg $DOWNGENE  >>$LOG 2>> $LOGTMP
    done

  fi
  echo "3: $R_PATH $DOCKER_PATH/geneDist.R --args $OUTSTAT $OUTPUT $LOG $CONTROLSTAT $PDF $BOOTSTRAP " >>$LOGTMP
  $R_PATH $DOCKER_PATH/geneDist.R --args $OUTSTAT $OUTPUT $LOG $CONTROLSTAT $PDF $BOOTSTRAP 2>>$LOGTMP >>$LOGTMP
else
  echo "4: $R_PATH $DOCKER_PATH/geneDist.R --args $OUTSTAT $OUTPUT $LOG $PDF" >>$LOGTMP
  $R_PATH $DOCKER_PATH/geneDist.R --args $OUTSTAT $OUTPUT $LOG $PDF 2>>$LOGTMP >>$LOGTMP
fi

if [ -r $LOGTMP ]; then
  rm $LOGTMP
fi

if [ -r $OUTPUT.annotated ]; then
  rm $OUTPUT.annotated*
fi


#if [ -r $CONTROLFILE.tmp ]; then
#  rm $CONTROLFILE.tmp*
#fi

