#!/bin/bash

: <<'hey'

 check the logs in case of error:
 	- go to : /bioinfo/http/prod/hosted/nebula.curie.fr/galaxy-dist/database/files/XXX
 	- see log1 (.log) : run stdout
	- see log2 (tmp.log): echos, montage log..
	
PWM from the mono version are reversed in summary file compared to what's in the raw log, this is done for compatibility with AhoPro.
diPWM are kept as in the raw log

hey


while getopts "f:n:s:m:v:o:z:i:x:r:t:" optionName; do
case "$optionName" in

f) inputfile="$OPTARG";;
n) NUMBER="$OPTARG";; #nbr of motifs to search
s) VERSION="$OPTARG";;
m) minW="$OPTARG";;
v) maxW="$OPTARG";;
o) LOG="$OPTARG";; #raw chipmunk output file
z) MODE="$OPTARG";; #mask or filter
i) OUTPNG="$OPTARG";;
x) NAME="$OPTARG";;
r) summary="$OPTARG";;
t) seqType="$OPTARG";;
esac
done


# my tmp log for debug needs
TMPLOG1=$OUTPNG.log
TMPLOG2=$OUTPNG.tmp.log
echo "$@" > $TMPLOG2


### Create tmp working dir :

OUTDIR=`mktemp -d`
local_path=/usr/bin/chipmunk/ChIPMunk_6.0a

echo $OUTDIR > $TMPLOG2
#chmod 777 $TMPLOG2

if [ -d $OUTDIR ]; then
	echo "Directory $OUTDIR exists" >> $TMPLOG2
	chmod 777 $OUTDIR
else 
	mkdir $OUTDIR
	
fi 
chmod -R 777 $OUTDIR

#move to tmp dir

cd $OUTDIR


#bluid motif length argument
motifs="$minW:$maxW"
for (( c=1; c<$NUMBER; c++ ))
do
   motifs="$motifs,$minW:$maxW"
done


if [ "$VERSION" == "Mono" ]; then
##run run_chiphordre
# optional option after 's:' <try_limit> <step_limit> <iter_limit> <thread_count>, put only defaut value in order to acess <thread_count> which is set to 4

	
	echo " ruby $local_path/run_chiphorde6.rb $NAME $motifs $MODE yes 1.0 s:$inputfile 100 10 1 4" >>  $TMPLOG2
	ruby $local_path/run_chiphorde6.rb $NAME $motifs $MODE yes 1.0 $seqType:$inputfile 100 10 1 4 1>> $TMPLOG1 2>&1 
	
## 1ST OUTPUT : raw log
	mv $NAME\_chiphorde.log $LOG 
	
			
else 	
##run run_chiphordre


	echo "ruby $local_path/run_dichiphorde6.rb $NAME $motifs $MODE yes 1.0 s:$inputfile 200 20 1 4" >> $TMPLOG2
	ruby $local_path/run_dichiphorde6.rb $NAME $motifs $MODE yes 1.0 $seqType:$inputfile 200 20 1 4 1>> $TMPLOG1 2>&1
	
## 1ST OUTPUT : raw log
	mv $NAME\_dichipmunk.log $LOG #raw log
	
fi 

	## get info line for summary file: INFO|found X motifs, used Y sequences
	INFO=$( awk '/^INFO\|/ {$1=""; print $0}' $LOG )
	
	
	## get nbr of motif X
	N=$( awk '/^INFO\|/ {split($0,tab); print tab[2] }' $LOG )
	
	## get threshold of each motif found : stored in variable Tx / x={1,2,....N}
	cat $LOG |  awk ' /^THRE\|/ { print $0}' | awk -v n=$N 'BEGIN {c=1} { if (NR <=n) print > "threshold_"c++}'

	
## 2ND OUTPUT : build the motif's png : use montage command
	
	#NB: montage -fill_image colxline -mode framing_style input_images.xml.png final_logo.png
	myPng=$NAME\_0.png
	us='_'
	
	for (( i=1; i< $N; i++ ))
	do
   		myPng="$myPng $NAME$us$i.png"
	done
	
	echo " myPng=$myPng" >> $TMPLOG2
	echo "	montage -tile 1x$N -mode Concatenate $myPng $OUTPNG " >> $TMPLOG2
	montage -tile 1x$N -mode Concatenate $myPng $OUTPNG  
	



## 3RD OUTPUT : build summary file (run infos + results summary)

	cat <<-SUM > $summary
	*** Run summary information ***
	
	ChiPMunk version used : $VERSION.ChIPMunk V6 17052014
	Number of different motifs to search : $NUMBER
	Min:Max width used : $minW:$maxW 
	Mode used : $MODE
	
	*** Results summary ***
	
	$VERSION.ChIPMunk found $INFO.
	PWMs of each motif found, with respective threshold :
	
	SUM
	
	#now, write (d)pcm in summary : motif 1, 2, 3..
	
	for (( i=0 ; i< $N ; i++))
		do	
			j=$((i+1))
			echo "#Motif $j:" >> $summary
			## append pwm to summary . In case of mono, first reverse pwm (hori -> verti) for compatibility with AhoPro	
			if [ $VERSION == "Mono" ]; then
				# append  pmw to summary 
				cat $NAME$us$i.pwm >> $summary
			else
				# append dpwm
				awk -v pwm="$NAME$us$i.dpwm" -f $local_path/addCol.awk $local_path/di_nuc.txt >> $summary
				#cat "$NAME$us$i.dpcm" >> $summary
			fi
			# append threshold to summary
			cat threshold_$j >> $summary
			echo -e "\n" >> $summary
		
		done
echo "For detailed results, see ChIPMunk detailed log." >> $summary

rm $NAME*
rm threshold*
#if [ -r $TMPLOG1 ]; then
#  rm $TMPLOG
#fi

#if [ -r $TMPLOG2 ]; then
#  rm $TMPLOG
#fi


