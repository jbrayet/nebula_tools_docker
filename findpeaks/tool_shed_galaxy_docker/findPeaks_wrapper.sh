#!/bin/bash


while getopts "r:p:w:a:m:h:l:i:k:n:s:y:e:" optionName; do
case "$optionName" in

r) OLOG="$OPTARG";;
p) OPEAKS="$OPTARG";;
w) OWIG="$OPTARG";;
a) ALIGNER="$OPTARG";;
m) MED="$OPTARG";;
h) HIGH="$OPTARG";;
l) LOW="$OPTARG";;
i) INPUT="$OPTARG";;
k) MIN="$OPTARG";;
n) NAME="$OPTARG";;
s) SUBPEAKS="$OPTARG";;
y) WIGSTEP="$OPTARG";;
e) PREPEND="$OPTARG";;

esac
done

#echo "$@" > /data/tmp/tmp.log

PATH_FP="/usr/bin/findpeaks/VancouverPackage-4.0.9.2/jars/fp4"

PREPS=" -prepend chr"

if [[ $PREPEND == "No" ]]
 then
   PREPS=""
fi

PATH2READ=`readlink $INPUT`

if [[ $PATH2READ == "" ]]
 then
    PATH2READ=$INPUT
fi


#hardcode outputPath


OUTDIR=`mktemp -d`

if [ -d $OUTDIR ]; then
	echo "Directory $OUTDIR exists" >>/dev/null
	echo "log will be written to $OUTDIR/$NAME.log" >>/dev/null
else 
	mkdir $OUTDIR
fi 


#echo "$@" > /data/tmp/tmp.log


FORMAT=$ALIGNER
# FILE=$INPUT

#sort $INPUT
echo "file $PATH2READ| grep -c gzip >$OUTDIR/ifGZIPed.tmp.tmp.txt" >>/dev/null
file $PATH2READ| grep -c gzip >$OUTDIR/ifGZIPed.tmp.tmp.txt
IFBAM=( $( cat $OUTDIR/ifGZIPed.tmp.tmp.txt))
if [ -r $OUTDIR/ifGZIPed.tmp.tmp.txt ]; then
  rm $OUTDIR/ifGZIPed.tmp.tmp.txt
fi
FILE=$OUTDIR/INPUT.$RANDOM.tmp.sam
TMPSORTED=$OUTDIR/INPUT.$RANDOM.sorted
if [[ $ALIGNER == "sam" && $IFBAM == 0 ]]; then

  #echo "sam to bam" >>/data/tmp/tmp.log
  samtools view -S -h $INPUT -b | samtools sort -m 4000000000 -o - $TMPSORTED 2>> $OUTDIR/$NAME.log | samtools view -h - >$FILE 2>> $OUTDIR/$NAME.log
fi
if [[ $ALIGNER == "sam" && $IFBAM == 1 ]]; then

 #echo "sort bam" >>/data/tmp/tmp.log
 echo "samtools sort -m 4000000000 -o $INPUT $TMPSORTED 2>> $OUTDIR/$NAME.log | samtools view -h - >$FILE 2>> $OUTDIR/$NAME.log" >>/dev/null
 samtools sort -m 4000000000 -o $INPUT $TMPSORTED 2>> $OUTDIR/$NAME.log | samtools view -h - > $FILE 2>> $OUTDIR/$NAME.log
fi
 
echo "java -Xmx20G -jar $PATH_FP/FindPeaks.jar -aligner $FORMAT -duplicatefilter -no_peaks_header$PREPS -dist_type 1 $MED $HIGH $LOW -input $FILE -minimum $MIN -name $NAME -output $OUTDIR -subpeaks $SUBPEAKS -wig_step_size $WIGSTEP 2>> $OUTDIR/$NAME.log >> $OUTDIR/$NAME.log" >>/dev/null
 
java -Xmx20G -jar $PATH_FP/FindPeaks.jar -aligner $FORMAT -duplicatefilter -no_peaks_header$PREPS -dist_type 1 $MED $HIGH $LOW -input $FILE -minimum $MIN -name $NAME -output $OUTDIR -subpeaks $SUBPEAKS -wig_step_size $WIGSTEP 2>> $OUTDIR/$NAME.log >> $OUTDIR/$NAME.log

 
if [ -r $FILE ]; then
  rm $FILE
fi

if [ -r $TMPSORTED ]; then
  rm $TMPSORTED
fi


    #rm $OLOG
mv $OUTDIR/$NAME.log $OLOG
    #rm $OPEAKS
mv $OUTDIR/$NAME\_triangle_subpeaks.peaks $OPEAKS
    #rm $OWIG
mv $OUTDIR/$NAME\_triangle_subpeaks.wig.gz $OWIG

