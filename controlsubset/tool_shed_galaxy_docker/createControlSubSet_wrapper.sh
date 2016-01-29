#!/bin/bash

CONTROLFILE=""

while getopts "f:c:t:o:s:q:p:" optionName; do
case "$optionName" in

f) CHIPFILE="$OPTARG";;
c) CONTROLFILE="$OPTARG";;
t) TYPE="$OPTARG";;
s) CHIPOUTPUT="$OPTARG";;
o) CONTROLOUTPUT="$OPTARG";;
q) IFPROCSAMPLE="$OPTARG";;
p) IFOUTSAM="$OPTARG";;

esac
done

local_path=$(dirname $(readlink -f $0))

LOGTMP=$IFOUTSAM.log.tmp

if [[ $IFPROCSAMPLE == "Yes" ]]
 then
  perl $local_path/createControlSubSet.pl -f $CHIPFILE -c $CONTROLFILE -t $TYPE -s $CHIPOUTPUT -o $CONTROLOUTPUT >$LOGTMP 2>>$LOGTMP
else 
  perl $local_path/createControlSubSet.pl -f $CHIPFILE -c $CONTROLFILE -t $TYPE -o $CONTROLOUTPUT  >>$LOGTMP 2>>$LOGTMP
fi

if [[ $TYPE == "bam" && $IFOUTSAM == "No" ]]
 then

    samtools view -S -b $CONTROLOUTPUT 2>>$LOGTMP | samtools sort -m 4000000000 - $CONTROLOUTPUT 2>>$LOGTMP
    rm $CONTROLOUTPUT
    mv $CONTROLOUTPUT.bam $CONTROLOUTPUT

  if [[ $IFPROCSAMPLE == "Yes" ]]
   then   
    samtools view -S -b $CHIPOUTPUT 2>>$LOGTMP | samtools sort -m 4000000000 - $CHIPOUTPUT 2>>$LOGTMP
    rm $CHIPOUTPUT
    mv $CHIPOUTPUT.bam $CHIPOUTPUT
  fi 
fi

head -2 $LOGTMP
if [ -r $LOGTMP ]; then
  rm $LOGTMP
fi
