#!/usr/bin/env python

import sys, subprocess, tempfile, os, shutil, glob

# samtools is added to PATH here

def calculate_window_size ( chr_len_file, input_file, format ):
	''' this function calculates the window size used to evaluate copy number variation.
	 
	W = L/(T*CV*CV)
	
	L=genome length
	T=total read count (mapped) present in the input file (bed, bam, sam)
	CV=coefficient of variation (0.05)
	W=window size ''' 
	
	cv=0.05
	#float b/c calculs is done in 'virgule flottante'
	proc1=subprocess.Popen ("awk '{sum=sum+$2} END {print sum}' %s" % chr_len_file, shell=True, stdout=subprocess.PIPE)
	genome_len = float (proc1.communicate()[0])
	if (proc1.returncode !=0):
		raise Exception ("In calculate_window_size(), returncode of proc1 (genome_len) = i%" % proc1.returncode)
	
	if format == "bed" :
	
		#proc_RC_bed=subprocess.Popen ("wc -l < %s > lines;read nbr_reads < lines;echo $nbr_reads" % input_file, shell=True, stdout=subprocess.PIPE)
		proc_RC_bed=subprocess.Popen ("grep -c 'chr' < %s" % input_file, shell=True, stdout=subprocess.PIPE)
		read_count = int (proc_RC_bed.communicate()[0] )
	
		#print ("read_count = %i\n" % read_count)	
		if (proc_RC_bed.returncode !=0):
			raise Exception ("In calculate_window_size(), returncode of proc_RC_bed (read_count) = i%" % proc_RC_bed.returncode)
	
	elif format == "bam":		
		# get index, get read count (can't pipe into index)
		proc_index=subprocess.Popen ("samtools index %s" % input_file, shell=True)
		#print "samtools index "+str(input_file)
		proc_index.wait()
		proc_RC_bam=subprocess.Popen ("samtools idxstats %s | awk '{s=s+$3} END {print s}'" % input_file, shell=True, stdout=subprocess.PIPE )
		read_count=int(proc_RC_bam.communicate()[0])

	else :
		# format = sam : convert to bam, get index, get read count
		proc_SamToBam=subprocess.Popen ("samtools view -bSh %s > %s.bam" % (input_file,input_file), shell=True) 
		proc_SamToBam.wait()
		proc_index=subprocess.Popen ("samtools index %s.bam" % input_file, shell=True)
		proc_index.wait()
		# get read count
		proc_RC_bam=subprocess.Popen ("samtools idxstats %s.bam | awk '{s+=$3} END {print s}'" % input_file, shell=True, stdout=subprocess.PIPE )
		read_count=int(proc_RC_bam.communicate()[0])		
	
	
	# w's type is float, & round returns float.0
	return int(round(genome_len/(read_count *cv*cv )))


def correct_hg19 (hg19Len,corrected_hg19Len ):
	''' this function corrects hg19.len file format : adds tab between the columns instead of space '''	
		
	i=open(hg19Len)
	o=open(corrected_hg19Len, "w")
	for line in i.readlines():
		col=line.split(" ")
		col.insert(1, "    ")
		o.write("".join(col))
	
	i.close()
	o.close()

	return 
	
	
def main():

	input_chip_file = sys.argv[1]
	input_control_file = sys.argv[2]
	hmcan_config_file = sys.argv[3]
	gccount_config_file = sys.argv[4]
	project_name = sys.argv[5]
	output_peaks_file = sys.argv[6]
	output_regions_file = sys.argv[7]
	output_density_file = sys.argv[8]
	output_posterior_proba_file = sys.argv[9]
	hmcan_log_report = sys.argv[10]
	chr_len_file = sys.argv[11]
	format = sys.argv[12]
	genome = sys.argv[13]
	
    hmcanconf=open(hmcan_config_file,'r')
    for line in hmcanconf:
        print(line)
    hmcanconf.close()
    print("sefsedsfs")
    gccount=open(gccount_config_file,'r')
    for line in gccount:
        print line

    #binary files
    GCCOUNT="/usr/bin/HMCan/HMCanV1.20/Utils/GCCount/gccount"
    HMCAN="/usr/bin/HMCan/HMCanV1.20/src/HMCan"

    print(chr_len_file)
    ###### CREATE ANNOTATION FILES #########
    process=subprocess.Popen('find / -type d -name files | grep database', shell=True, stdout=subprocess.PIPE)
    databasePath = (process.communicate()[0]).replace("\n","")
    print(databasePath)

    subprocess.Popen('mkdir -p '+databasePath+'/nebulaAnnotations', shell=True)
    subprocess.Popen('mkdir -p '+databasePath+'/nebulaAnnotations/'+genome, shell=True)

    nebulaGenomePath=databasePath+"/nebulaAnnotations/"+genome
    print(nebulaGenomePath)

    FAIFILE='n'
    LENFILE='n'
    DICTFILE='n'
    CHROFILE='n'
    MAPFILE='n'

    if not os.path.isdir(nebulaGenomePath+"/chromosomes"):
            subprocess.Popen('mkdir -p '+nebulaGenomePath+'/chromosomes', shell=True)
            CHROFILE='y'

    if not os.path.isfile(nebulaGenomePath+"/"+genome+".len"):
            FAIFILE='y'
            LENFILE='y'

    if not os.path.isfile(nebulaGenomePath+"/out50m2_"+genome+".gem.mappability"):
            MAPFILE='y'
    
    print(str(CHROFILE)+" "+str(FAIFILE)+" "+str(LENFILE)+" "+str(MAPFILE))

    FILEPATH=databasePath.replace("database/files","tool-data")

    cmd='bash create_annotation_files.sh '+FAIFILE+" "+LENFILE+" "+DICTFILE+" "+CHROFILE+" "+FILEPATH+" "+genome+" "+MAPFILE+" "+nebulaGenomePath
    print(cmd)
    process=subprocess.Popen(cmd, shell=True)
    process.wait()
    ############### END ##############
	
	# create tmp dir.Returns absolute path to a created tmp_dir  
	tmp_dir=tempfile.mkdtemp()
	#tmp_dir="/data/tmp/HMcanTest"
	
	try:

## GCCOUNT
		
		os.chdir(tmp_dir)
		
		#add 'window size' and 'step' to gccount config file
		window_size=calculate_window_size (chr_len_file, input_control_file, format)

		#window_size=10000 
		#print ("window_size= %i\n" % window_size)


#************** EDIT config files (improved)*****************		
		
		#Frist, correct hg19.len file, in case HMCan is run for hg19 (because hg19.len is not tabular on the server, better change it here)
		#hg19="hg19"
		#if (hg19 in chr_len_file):
		#	correct_hg19 (chr_len_file,"corrected_hg19Len_file.len")	
		#	param_len="corrected_hg19Len_file.len"
		
		#else:
			#to make my life easier..
		param_len= chr_len_file

		#Edit gc count config file
		cmd_step ="sed -i \"/step/c\step = %s\" %s" % ( window_size , gccount_config_file )
		cmd_window ="sed -i \"/window/c\window = %s\" %s" % ( window_size , gccount_config_file )
		cmd_chrLen = "sed \"s~chrLenFile =.*~chrLenFile = "+nebulaGenomePath+"/"+genome+".len~g\" "+gccount_config_file
		cmd_chrFiles = "sed \"s~chrFiles =.*~chrFiles = "+nebulaGenomePath+"/chromosomes~g\" "+gccount_config_file
		cmd_gemMappabilityFile = "sed \"s~gemMappabilityFile =.*~gemMappabilityFile = "+nebulaGenomePath+"/out50m2_"+genome+".gem.mappability~g\" "+gccount_config_file

		
		s=subprocess.Popen(args=cmd_step, shell=True)
		s.wait()
		
		w=subprocess.Popen(args=cmd_window, shell=True)
		w.wait()
			
		cL=subprocess.Popen(args=cmd_chrLen, shell=True)
		cL.wait()

		cF=subprocess.Popen(args=cmd_chrFiles, shell=True)
		cF.wait()
			
		cM=subprocess.Popen(args=cmd_gemMappabilityFile, shell=True)
		cM.wait()

#*********** END edit config files ****************************		

				
		##call gccount  /dev/null
		gc_proc=subprocess.Popen(args= "%s %s >> %s 2>&1" % (GCCOUNT, " ".join(['-conf', gccount_config_file]), hmcan_log_report), shell=True, stderr=subprocess.PIPE)
		gc_stderr_data = gc_proc.communicate()[1] # waits for process to terminate, returns (stdoutdata, stderrdata)
		#if (gc_proc.returncode !=0): 
		#	raise Exception("GCcount returncode = %i\n GCcount stderr : %s\n" % (gc_proc.returncode,gc_stderr_data ))
		
		# Find outputs: .cnp files 
		# glob returns list of absolute paths of files.cnp in tmp_dir
		cnp_list=glob.glob('*.cnp')
		if (len(cnp_list) == 0): 
			raise Exception("Error while running Gccount : program did not generate outputs for HMCan (file.cnp) !\n GCCOUNT stderr : %s\n" % gc_stderr_data) 
		cnp_file=cnp_list[0]
		#print "GC PROFIL file :%s\n" % cnp_file
			
## HMCAN
		# EDIT hmcan config file : add 'GCIndex' and 'largeBinLength' 
		cmd_gc="sed -i \"/GCIndex/c\GCIndex %s\" %s" % ( cnp_file , hmcan_config_file )
		cmd_bin="sed -i \"/largeBinLength/c\largeBinLength %s\" %s" % ( window_size , hmcan_config_file )
		cmd_chrFiles = "sed \"s~chrFiles =.*~chrFiles = "+nebulaGenomePath+"/chromosomes~g\" "+hmcan_config_file

		g=subprocess.Popen(args=cmd_gc, shell=True)
		g.wait()
		
		b=subprocess.Popen(args=cmd_bin, shell=True)
		b.wait()
		
		cF=subprocess.Popen(args=cmd_chrFiles, shell=True)
		cF.wait()
		
		#call HMCan , hmcan_log_report
		hmcan_proc= subprocess.Popen(args= "%s %s >> %s 2>&1" % (HMCAN, " ".join([input_chip_file, input_control_file, hmcan_config_file , project_name]), hmcan_log_report ), shell=True, stderr=subprocess.PIPE)
		hmcan_stderr_data=hmcan_proc.communicate()[1]
		if (hmcan_proc.returncode !=0):
			raise Exception	("Error occured, HMCan returncode = %i\n HMCan stderr : %s\n" % ( hmcan_proc.returncode, hmcan_stderr_data ) )
	
		#deal with outputs : copy outputs from tmp_dir to dataset_XXX :
		
		#shutil.move( os.path.join (tmp_dir,"%s_peaks.bed" % project_name ), output_peaks_file) 
		shutil.move( os.path.join (tmp_dir,"%s_peaks.narrowPeak" % project_name ), output_peaks_file)
		shutil.move( os.path.join( tmp_dir, "%s_regions.bed" % project_name ) ,output_regions_file )
		if (output_posterior_proba_file !='None') :
			shutil.move( os.path.join( tmp_dir, "%s_posterior.wig" % project_name ) ,output_posterior_proba_file ) 
		if (output_density_file != 'None'):
			shutil.move( os.path.join( tmp_dir, "%s.wig" % project_name ) ,output_density_file )
	
	#Handle exceptions
	except OSError:
		sys.stderr.write ( "Error occured while runing HMCan : trying to execute a non existing file.\n" )
		#remove tmp_dir created:
		if os.path.exists( tmp_dir ):
			shutil.rmtree(tmp_dir)
		raise
		
	except ValueError:
		sys.stderr.write ( "Error occured while runing HMCan :  Popen is called with invalid arguments.\n" )
		if os.path.exists( tmp_dir ):
			shutil.rmtree(tmp_dir)
		raise
	
	#Handle other	
	except Exception:
		sys.stderr.write ( "Error occured while runing HMCan\n" )		
		if os.path.exists( tmp_dir ):
			shutil.rmtree(tmp_dir)
		raise
		
	#clean up my mess 
	if os.path.exists( tmp_dir ):
			shutil.rmtree(tmp_dir)


#main

if __name__=='__main__': 
	main()
 
