#!usr/bin/bash 


:<<'hey'
Command line looks like this :
ahopro_wrapper.sh ${ahopro_config} ${motif_file} ${function.function_selector} $outfile $function['seq_name'] $function['nbr_motif'] $actualN
# if function = p-value :, commad line will looks like this : 
ahopro_wrapper.sh ${ahopro_config} ${motif_file} ${function.function_selector} $outfile $function['seq_name'] $function['nbr_motif'] $actualN ${letter_freq_file}


hey

#sort out arguments

configfile=$1
motif_file=$2
function=$3
#if function is "p-value", get the freqFile pramater
if [ "$function" == "p-value" ] ; then 
	letterFreqFile=$8
fi
output=$4
SEQ_NAME=$5
NBR_MOTIFS=$6
actual_NBR=$7



PATH_AHOPRO="/usr/bin/ahopro/AhoPro.1.3/"


# create tmp workindg dir and cd

OUTDIR=`mktemp -d`

#hard code the creation of OUTDIR
if [ -d $OUTDIR ]; then
	chmod -R 777 $OUTDIR
else 
	mkdir $OUTDIR
fi 

chmod -R 777 $OUTDIR

cd $OUTDIR



#echo "######configfile out of XML form :" > $output
#cat $configfile >> $output
#echo "motifFile :" >> $output
#cat $motif_file > $output




## ************ First, check if the number of actual motif is met in the entries. If not, warn, and quit ******************************************************

if [[ $NBR_MOTIFS != $actual_NBR ]]; then 
	
cat <<doc >$output
AhoPro error :

The number of motifs set in  "Number of motifs" does not match the the actual number of motif entry in "Motifs" section. They should be the same. 

Number of motifs : $NBR_MOTIFS
Number of actual motif entered : $actual_NBR

Please recheck these parameters.	

doc

exit 0

fi

## ************ END check ******************************************************



### **************** Deal with motif entery ****************************************************************************************************


if [[ $NBR_MOTIFS = 1 ]]; then
	
	#NO PARSING REQUIRED
	
	#remove blank lines, tabs, spaces..
	cat $motif_file | awk '/^[0-9\-ACTGatgc]/ {print}' > motif_tmp_0.txt
	# reverse if the motif is a pmw, reverse it if necessary
	type=`cat motif_tmp_0.txt | awk ' END { if ($1 ~ /[acATGC]/) {print "words"} else {print "pmw"} }'`
	
	if [[ "$type" == "pmw" ]]; then

		NF=`cat motif_tmp_0.txt | awk ' END {print NF}'`
	
		if [[ $NF != 4 ]]; then 	#if it's not vertical, reverse it 
	
			for ((c=1 ; c <= $NF; c++)) 
				do 
					cat motif_tmp_0.txt | awk -v c=$c '{ printf $c" " } END {printf"\n"}' >> motif_0.txt			
				done 
					
			rm motif_tmp_0.txt
			
		
		elif [[ $NF == 0 ]]; then 
			echo "Please check the motif entries.Try removing tabulations, spaces." > $output
			exit 0
			
		else #if vertical 
			mv motif_tmp_0.txt motif_0.txt
		fi
	
	else #if words
			mv motif_tmp_0.txt motif_0.txt

	fi


else # IF there is more than one motif 
	
	#check if it's true + parse motifs + warn
	N=`cat $motif_file | awk '/^[0-9\-ACTGatgc*]/ {print}' | awk ' BEGIN{RS="*"; n=0} {printf $0 > "motif_tmp_"n++} END {print n}'`

	# reverse pmw (or not if it's vertical)

	for ((i=0 ; i< $N ; i++)) 
		do 
			# remove empty lines + check if it's a pmw, motif words.
			cat motif_tmp_$i | awk '/^[0-9\-ACTGatgc]/ {print}' > motif_tmp_tmp_$i
			type=`cat motif_tmp_tmp_$i | awk ' END { if ($1 ~ /[acATGC]/) {print "words"} else {print "pmw"} }'` 
	
			rm motif_tmp_$i
			# if pmw : check if it is horizontal, if so, reverse it (ahopro takes only vertical pmws !)
	
			if [[ "$type" == "pmw" ]]; then

				NF=`cat motif_tmp_tmp_$i | awk ' END {print NF}'`
	
				if [[ $NF != 4 ]]; then 	#if it's not vertical, reverse it 
	
					for ((c=1 ; c <= $NF; c++)) 
						do 
							cat motif_tmp_tmp_$i | awk -v c=$c '{ printf $c" " } END {printf"\n"}' >> motif_$i.txt			
						done 
					rm motif_tmp_tmp_$i
			
		
				elif [[ $NF == 0 ]]; then 
					echo "Please check the motif entries, try removing tabulations, spaces." > $ouput
					exit 0
			
				else #if vertical 
					mv motif_tmp_tmp_$i motif_$i.txt
				fi
	
			else #if words
					mv motif_tmp_tmp_$i motif_$i.txt
	
			fi

		#cat motif_$i.txt 
	done
fi

## ****************************************** END motif entery treatment *********************************************************************************************************************

###motif_0
# next could be improved in bash (sed..)

#Add motif filenames to the configfile byt writing a new_config.txt (in CWD)
##but first, check what kind of correction ? basic : motifnames only ? or LetterFreqFile also ?

if [ "$function" == "p-value" ] ; then 

python $PATH_AHOPRO/replaceLine.py $configfile "name_freq" $letterFreqFile

else

python $PATH_AHOPRO/replaceLine.py $configfile "name" 

fi

## get the new config file
#new_config="$OUTDIR/new_config.txt"
new_config="new_config.txt"

#echo "******** AhoPro Parameters ******** " >> $output
#cat $new_config >> $output
#call ahopro with config file


#echo "******** AhoPro Results ******** " >> $output

$PATH_AHOPRO/ahokocc $new_config >> $output 2>&1

# improve line "Search in sequence .." with user entered sequence name (for more clarity). Do only if seq_name != 'X' 

if [ "$SEQ_NAME" != "X" ] ; then 

	#remove sequence (if it exists)
sed -i '/^Search/d' $output
	# add new line
line="\nSearch in sequence : $SEQ_NAME\n"
sed -i "1s/^/$line/" $output

fi

rm $new_config

