#!/bin/bash



:<< 'HEYthere'
* check the logs in case of error:
 	- go to : ../galaxy/database/files/XXX
 	- see log.tmp (.tmp): step echos, stderr of all commands... 
 
* working directory for MICSA is : $OUTDIR/wig, contains the ouputs and logs 

* xml command :
micsa_wrapper.sh -i ${chip_file} -c ${control_file} -g $genome['genome_path'] -b $genome.blacklist['blacklist_file'] -r $ratio -n ${project_name} -f ${format_selector} -p $fdr_or_fr['cmd_option'] 

-v $fdr_or_fr['fpp'] ORR -v $fdr_or_fr['fdr'] 

-w ${selected_peaks_wig_gz} -t ${selected_peaks_txt} -l ${motifs_words_out} -m ${micsa_log} -o $outpng


** arguments to add if optional length parameters are added :

in xml : -x ${min_len} -y ${men_len} -z ${med_len}

in wrapper.sh :

x) LOW="$OPTARG";;
y) MED="$OPTARG";;
z) HIGH="$OPTARG";;

add to findpeaks command after -dist_type 1

HEYthere

	
while getopts "i:c:g:b:r:n:f:p:v:w:t:l:m:o:x:y:z:d:a:" optionName; do
case "$optionName" in

i) PATH_CHIP="$OPTARG";;
c) PATH_CONTROL="$OPTARG";;
g) GENOME="$OPTARG";;
b) BLACK_LIST="$OPTARG";;
r) RATIO="$OPTARG";;
n) NAME="$OPTARG";;
f) FORMAT="$OPTARG";;

p) OPTION="$OPTARG";;
v) VALUE="$OPTARG";;

w) selected_peaks_wig_gz="$OPTARG";;
t) selected_peaks_txt="$OPTARG";;
l) motif_words_out="$OPTARG";;
m) LOG="$OPTARG";;
o) OUTPNG="$OPTARG";;

x) LOW="$OPTARG";;
y) MED="$OPTARG";;
z) HIGH="$OPTARG";;

d) ROOT_DIR="$OPTARG";;
a) BUILD="$OPTARG";;
esac
done

PATH_MICSA="/usr/bin/micsa/MICSA"
PATH_FP="/usr/bin/micsa/VancouverPackage-4.0.9.2/jars/fp4/"

DATABASE_PATH=$ROOT_DIR/database/files
mkdir -p $DATABASE_PATH/nebulaAnnotations
mkdir -p $DATABASE_PATH/nebulaAnnotations/$BUILD
nebulaAnnotationPath=$DATABASE_PATH/nebulaAnnotations

############### Create annotations files ################

FAIFILE='n'
LENFILE='n'
DICTFILE='n'
CHROFILE='n'
MAPFILE='n'

if [ ! -d $DATABASE_PATH/nebulaAnnotations/$BUILD/chromosomes ]; then
    mkdir -p $DATABASE_PATH/nebulaAnnotations/$BUILD/chromosomes
    CHROFILE='y'
fi

FILEPATH=$ROOT_DIR/tool-data

bash /usr/bin/create_annotation_files.sh $FAIFILE $LENFILE $DICTFILE $CHROFILE $FILEPATH $BUILD $MAPFILE $nebulaGenomePath

if [ ! -d $DATABASE_PATH/nebulaAnnotations/$BUILD/blacklist ]; then
    mkdir -p $DATABASE_PATH/nebulaAnnotations/$BUILD/blacklist
    if [[ "$BUILD" == "hg18" ]]; then
        wget http://zerkalo.curie.fr:8080/partage/nebulaAnnotation/blacklist/dacDuke_hg19blacklist_lifted_to_hg18.bed -P $DATABASE_PATH/nebulaAnnotations/$BUILD/blacklist
        wget http://zerkalo.curie.fr:8080/partage/nebulaAnnotation/blacklist/dac_hg19blacklist_lifted_to_hg18.bed -P $DATABASE_PATH/nebulaAnnotations/$BUILD/blacklist
        wget http://zerkalo.curie.fr:8080/partage/nebulaAnnotation/blacklist/duke_blacklist_hg18.bed -P $DATABASE_PATH/nebulaAnnotations/$BUILD/blacklist
    fi
    if [[ "$BUILD" == "hg19" ]]; then
        wget http://zerkalo.curie.fr:8080/partage/nebulaAnnotation/blacklist/dacDuke_blacklist_hg19.bed -P $DATABASE_PATH/nebulaAnnotations/$BUILD/blacklist
        wget http://zerkalo.curie.fr:8080/partage/nebulaAnnotation/blacklist/dac_blacklist_hg19.bed -P $DATABASE_PATH/nebulaAnnotations/$BUILD/blacklist
        wget http://zerkalo.curie.fr:8080/partage/nebulaAnnotation/blacklist/duke_blacklist_hg19.bed -P $DATABASE_PATH/nebulaAnnotations/$BUILD/blacklist
    fi
fi

BLACK_LIST_PATH=$DATABASE_PATH/nebulaAnnotations/$BUILD/blacklist

#################### END ANNOTATION FILES ###########################

#### CREATE A SAFE TEMP OUTDIR for  Nebula

OUTDIR=`mktemp -d`

echo $OUTDIR > $LOG.tmp

#hard code the creation of OUTDIR
if [ -d $OUTDIR ]; then
	echo "Directory $OUTDIR exists" >> $LOG.tmp
	#echo "##### before micsa  PATH=$PATH" >> $LOG.tmp
	chmod 777 $OUTDIR
else 
	mkdir $OUTDIR
fi 

mkdir $OUTDIR/wig
mkdir $OUTDIR/sep_chip
mkdir $OUTDIR/sep_control

#new for hg19
mkdir $OUTDIR/chromosomes

chmod -R 777 $OUTDIR


# ***************************************** new : if hg19 was chosen, select only the regular chromosomes

if [[ "$BUILD" == "hg19" ]]; then

	#create symbolic links to 'regular' chromosomes
	ln -s $nebulaAnnotationPath/$BUILD/chromosomes/chr??.fa $nebulaAnnotationPath/$BUILD/chromosomes/chr?.fa $OUTDIR/chromosomes
	GENOME=$OUTDIR/chromosomes

fi
#************************************************************************************************************


############################ DEAL WIH DIFFERENT FORMATS ##### + CALL FINDPEAKS ################################################

TMPBAMSORTED_CHIP=$OUTDIR/tmpbam.chip.$RANDOM.sorted
TMPBAMSORTED_CONTROL=$OUTDIR/tmpbam.control.$RANDOM.sorted

FILE_CHIP=$OUTDIR/input.$RANDOM.sorted.sam
FILE_CONTROL=$OUTDIR/input.$RANDOM.sorted.sam

if [[ "$FORMAT" == "bed" ]]; then
	
	echo "Starting SeparateReads for BED..." >> $LOG.tmp

	java -jar $PATH_FP/SeparateReads.jar $FORMAT $PATH_CHIP $OUTDIR/sep_chip > /dev/null
	java -jar $PATH_FP/SeparateReads.jar $FORMAT $PATH_CONTROL $OUTDIR/sep_control > /dev/null

	echo "Starting Sort reads for BED..." >> $LOG.tmp

	java -jar $PATH_FP/SortFiles.jar $FORMAT $OUTDIR/ $OUTDIR/sep_chip/*.part.bed.gz > /dev/null
	java -jar $PATH_FP/SortFiles.jar $FORMAT $OUTDIR/ $OUTDIR/sep_control/*.part.bed.gz  > /dev/null

	####### PS. sortReads overwrites the separated files


	echo "Starting Findpeaks for chip.bed ..." >> $LOG.tmp
	
	java -Xmx2G -jar $PATH_FP/FindPeaks.jar -aligner $FORMAT -duplicatefilter -input $OUTDIR/sep_chip/*.part.bed.gz -name 'chip' -output $OUTDIR/wig -dist_type 1 -minimum 3 > /dev/null 
	
	echo "Starting Findpeaks for control.bed ..." >> $LOG.tmp
																									
	java -Xmx2G -jar $PATH_FP/FindPeaks.jar -aligner $FORMAT -duplicatefilter -input $OUTDIR/sep_control/*.part.bed.gz -name 'control' -output $OUTDIR/wig -dist_type 1 -minimum 1 > /dev/null

	
fi

if [[ "$FORMAT" != "bed" ]]; then
	
	echo "Entering the SAM/BAM section" >> $LOG.tmp
	
	if [[ "$FORMAT" == "sam" ]]; then
		echo "Starting treatment for .SAM :  SAM -> BAM |sort| BAM_sorted -> SAM" >> $LOG.tmp
		#convert to BAM | sort BAM | convert back to SAM > FILE
		samtools view -S -h $PATH_CHIP -b 2> /dev/null| samtools sort -m 4000000000 -o - $TMPBAMSORTED_CHIP 2> /dev/null | samtools view -h - > $FILE_CHIP 
		samtools view -S -h $PATH_CONTROL -b 2> /dev/null | samtools sort -m 4000000000 -o - $TMPBAMSORTED_CONTROL 2> /dev/null | samtools view -h - > $FILE_CONTROL 
	fi
	
	if [[ "$FORMAT" == "bam" ]]; then
		echo "Starting treatment for .BAM : BAM |sort| BAM_sorted -> SAM" >> $LOG.tmp
		
		#first, get the headers :
		#/bioinfo/local/samtools/samtools view -H $PATH_CHIP > $OUTDIR/headerCHIP.sam
		#/bioinfo/local/samtools/samtools view -H $PATH_CONTROL > $OUTDIR/headerCONTROL.sam 
		
		# sort BAM | reheading sorted BAM | convert back to sam
		#/bioinfo/local/samtools/samtools sort -m 4000000000 -o $PATH_CHIP $TMPBAMSORTED_CHIP | /bioinfo/local/samtools/samtools reheader $OUTDIR/headerCHIP.sam -  | /bioinfo/local/samtools/samtools view -h - > $FILE_CHIP
		#/bioinfo/local/samtools/samtools sort -m 4000000000 -o $PATH_CONTROL $TMPBAMSORTED_CONTROL | /bioinfo/local/samtools/samtools reheader $OUTDIR/headerCONTROL.sam - | /bioinfo/local/samtools/samtools view -h - > $FILE_CONTROL

		
		#sort BAM | convert back to SAM > FILE
		samtools sort -m 4000000000 -o $PATH_CHIP $TMPBAMSORTED_CHIP 2> /dev/null | samtools view -h - > $FILE_CHIP
		samtools sort -m 4000000000 -o $PATH_CONTROL $TMPBAMSORTED_CONTROL 2> /dev/null | samtools view -h - > $FILE_CONTROL
	fi
	
	
	#both cases, call finpeaks -input $FILE_xxxx , format = sam
	echo "runing findpeaks now for SAM/BAM..." >> $LOG.tmp
	java -Xmx2G -jar $PATH_FP/FindPeaks.jar -aligner 'sam' -duplicatefilter -input $FILE_CHIP -name 'chip' -output $OUTDIR/wig -dist_type 1 $MED $HIGH $LOW -minimum 3 >> $LOG.tmp 2>&1
	java -Xmx2G -jar $PATH_FP/FindPeaks.jar -aligner 'sam' -duplicatefilter -input $FILE_CONTROL -name 'control' -output $OUTDIR/wig -dist_type 1 $MED $HIGH $LOW -minimum 3 >> $LOG.tmp 2>&1

echo "FINISHED FIND PEAKS FOR CONTROL/CHIP" >> $LOG.tmp

fi



########## out : chip_triangle_standard.peaks	chip_triangle_standard.wig.gz ############## DELETE REGIONS
echo "Starting Delete Regions..." >> $LOG.tmp

java -Xmx1565m -cp $PATH_MICSA DeleteRegions -f $OUTDIR/wig/chip_triangle_standard.peaks -r $BLACK_LIST_PATH/$BLACK_LIST >> $LOG.tmp #> /dev/null

java -Xmx1565m -cp $PATH_MICSA DeleteRegions -f $OUTDIR/wig/control_triangle_standard.peaks -r $BLACK_LIST_PATH/$BLACK_LIST >> $LOG.tmp # > /dev/null


####### writes in files.peaks ########################################### SUMMARY
echo "Starting Summary..." >> $LOG.tmp

# cd to /wig, because the command likes it ONLY like this

cd $OUTDIR/wig 

java -cp $PATH_MICSA Summary -f chip_triangle_standard.peaks -c control_triangle_standard.peaks -r $RATIO >> $LOG.tmp 2>&1


############################ creates file : $OUTDIR/wig/FindPeaksSummary.txt ########## FILTER peaks (control Vs chip)
echo "Starting filter peaks..." >> $LOG.tmp

java -cp $PATH_MICSA FilterPeaks -f chip_triangle_standard.peaks -c control_triangle_standard.peaks -t 3.5  >> $LOG.tmp 2>&1

########################### modifies only : chip_triangle_standard.peaks ########### MICSA : 
### CWD : $OUTDIR/wig
## NOTE: eventhough -o $OUTDIR is set for micsa, it will output its NAME.log 

echo "## WOKING DIRECTORY FOR MICSA :" >> $LOG.tmp
pwd >> $LOG.tmp

# meme.ben should be added to PATH
export PATH=$PATH:/usr/bin/meme/bin
echo "##### PATH=$PATH" >> $LOG.tmp
echo "Starting Micsa..." >> $LOG.tmp
echo "##### java -Xmx4G -jar $PATH_MICSA/micsa.jar -name $NAME -f $OUTDIR/wig/chip_triangle_standard.peaks $OPTION $VALUE -o $OUTDIR/wig -l $OUTDIR/wig/FindPeaksSummary.txt -g $DATABASE_PATH/nebulaAnnotations -w $OUTDIR/wig/chip_triangle_standard.wig.gz >> $LOG.tmp 2>&1 ">> $LOG.tmp

#double memory 2 to 4
java -Xmx4G -jar $PATH_MICSA/micsa.jar -name $NAME -f $OUTDIR/wig/chip_triangle_standard.peaks $OPTION $VALUE -o $OUTDIR/wig -l $OUTDIR/wig/FindPeaksSummary.txt -g $DATABASE_PATH/nebulaAnnotations -w $OUTDIR/wig/chip_triangle_standard.wig.gz >> $LOG.tmp 2>&1


############# ALL 3 OUTPUTS ARE IN $OUTPUTS/wig
echo "## ls CWD :" >> $LOG.tmp
ls >> $LOG.tmp

################## parse micsa's outputs found in $OUTDIR/wig 

#split motif.txt main output into seperate motifs AND store the number of motifs

motifs=$OUTDIR/wig/motifs.txt
N=`cat $motifs | awk '{   if ( $1 !~ /[*:]/ ) print } ' | awk '/./' | awk '{ if (NR !=1 ) print  }' | awk ' BEGIN{ RS = "MOTIF [0-9]" } { print > "motif_"++c".tmp.tmp.txt" } END {print c}'`

#generates motif_X.tmp.tmp.txt

#remove first line of all motif files, and generate logos AND remove files 
for (( i=1; i<= $N  ; i++))
	
	do

		#remove frist line
		#awk '{  if (NR==1) { {print " " ;} }  else  {print ; }  }' $OUTDIR/wig/motif\_$i.tmp.txt > $OUTDIR/wig/motif\_$i.tmp.tmp.txt
		#rm $OUTDIR/wig/motif\_$i.tmp.txt 1>> $LOG.tmp 2>&1
		
		#generate logs
		ruby /bioinfo/local/build/ChIPMunk_6.0/pmflogo3.rb $OUTDIR/wig/motif\_$i.tmp.tmp.txt $OUTPNG.$i.png 1> /dev/null 2>&1
		
		rm $OUTDIR/wig/motif\_$i.tmp.tmp.txt
		
	done

	
#create FINAL logo from OUTPNG.X.png
us="_"
myPng=$OUTPNG.1.png
for (( c=2; c<= $N; c++ ))
do
   myPng="$myPng $OUTPNG.$c.png"
done

echo "montage -tile 1x$N -mode Concatenate $myPng $OUTPNG.png" >> $LOG.tmp

#montage -fill_image colxline -mode framing_style input_images.xml.png final_logo.png
montage -tile 1x$N -mode Concatenate $myPng $OUTPNG.png >> $LOG.tmp 2>&1

mv $OUTPNG.png $OUTPNG

#remove individual logos
for (( c=1; c<= $N; c++ ))
do
   rm $OUTPNG.$c.png
done

################################## Selected the correct columns from $NAME file (selectedpeaks.txt)
awk 'BEGIN {print "Chromosome Start End CoverageScore PeakLength P-value Strand"} { if (NR !=1) { print $1, $2, $3, $5,  $6, $7, "+"} }' $OUTDIR/wig/$NAME | column -t > $selected_peaks_txt
rm $OUTDIR/wig/$NAME >> $LOG.tmp 2>&1

#### if .wig should be outputed, in xml put format="wig"
#gunzip $OUTDIR/SelectedPeaks.wig.gz
#mv $OUTDIR/SelectedPeaks.wig $selected_peaks_wig

#in xml format="gzip"
mv $OUTDIR/wig/SelectedPeaks.wig.gz $selected_peaks_wig_gz

if [ "$LOG" != "None" ]; then
	mv $OUTDIR/wig/$NAME.log $LOG
else 
	rm $OUTDIR/wig/$NAME.log
fi


if [ "$motif_words_out" != "None" ]; then
	echo "### mv $motifs $motif_words_out" >> $LOG.tmp 2>&1
	mv $motifs $motif_words_out
else 
	rm $motifs			
fi

echo "####removing $OUTDIR" >> $LOG.tmp
#remove created working dir

#OUTDIR is removed by torque ??


