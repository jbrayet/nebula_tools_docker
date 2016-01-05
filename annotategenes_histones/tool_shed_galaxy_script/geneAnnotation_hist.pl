#!/usr/bin/perl -w
#outputs statistics for all genes in the list

#different boundaries
#no motif p-value for binding sites
#read directly No-resp/Up/down category

#all isoforms from the file with genes

#RNApolII sites on junctions

#for HISTONES!!!!

use strict;
use POSIX;
use warnings;

my $usage = qq{
    $0

    -----------------------------
    mandatory parameters:
    
    -g   filename 	file with all genes
    -tf filename 	file with sites of TF1
    
    
    -----------------------------
    optional parameters:
    -k36 filename 	file with sites of K36me3 
    -rp  filename 	file with sites of RNApolII   
    -i   filename       file with a table where to add colomnes
    -add values         which colomns to add

    -o 		filename 	output filename (defaut "genes.output.txt") 
    -v 				verbose mode
    -mir 	filename	file with positions of miRNA
    -k9		filename 	ile with sites of K9me3     

    -c_rp	value		cutoff for -rp
    -c_k9	value		cutoff for -k9
    -c_k36	value		cutoff for -k36

    -selG       filename        selected genes (up-down-regulated)
    -fluo	filename	file with fluorescence
    -gc         filename        file with gc-islands
    
    -long                       for each gene take the longest isoform


};

if(scalar(@ARGV) == 0){
    print $usage;
    exit(0);
}

## mandatory arguments

my $RNApolFilename = "";
my $H3K36Me3polFilename = "";
my $H3K9Me3polFilename = "";
my $TF1Filename = "";
my $TF2Filename = "";
my $GenesFilename = "";
my $MirFilename = "";
my $TF1FilenameALL = "";
my $TF2FilenameALL  = "";
my $SelectedGenesFilename = "";
my $fluoFile = "";
my $initialTable = "";
my $colomnesToAdd = "";

my $enhLeft = -30000;
my $enhRight = -1500;
my $immediateDownstream = 2000;
my $K9dist = 5000;
my $kb5 = 5000;
my $INFINITY = 10000000000;
my $jonctionSize = 50;
## optional arguments
my $outname = "genes.output.txt";
my $verbose = 0;
my $GCislands = "";


my $longest = 0;

#my $cutoff_tf1 = 0;
#my $cutoff_tf2 = 0;
my $cutoff_tf1All = 0;
my $cutoff_tf2All = 0;
my $cutoff_rp = 0;
my $cutoff_k9 = 0;
my $cutoff_k36 = 0;
my $ifTFcoord = 0;

my $ifHMmode = 0;

## parse command line arguments

while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;
    if ( $this_arg eq '-h') {print "$usage\n"; exit; }

    elsif ( $this_arg eq '-selG') {$SelectedGenesFilename = shift @ARGV;}  

    elsif ( $this_arg eq '-g') {$GenesFilename = shift @ARGV;}
    elsif ( $this_arg eq '-rp') {$RNApolFilename = shift @ARGV;}
    elsif ( $this_arg eq '-k36') {$H3K36Me3polFilename = shift @ARGV;}
    elsif ( $this_arg eq '-k9') {$H3K9Me3polFilename = shift @ARGV;}

    elsif ( $this_arg eq '-tf') {$TF1Filename = shift @ARGV;}

    elsif ( $this_arg eq '-v') {$verbose = 1;}
    
    elsif ( $this_arg eq '-long') {$longest = 1;}

    
    elsif ( $this_arg eq '-o') {$outname = shift @ARGV;}
    elsif ( $this_arg eq '-mir') {$MirFilename = shift @ARGV;}
  
    elsif ( $this_arg eq '-c_rp') {$cutoff_rp = shift @ARGV;}
    elsif ( $this_arg eq '-c_k9') {$cutoff_k9 = shift @ARGV;}
    elsif ( $this_arg eq '-c_k36') {$cutoff_k36 = shift @ARGV;}
    
    elsif ( $this_arg eq '-fluo') {$fluoFile = shift @ARGV;}
    elsif ( $this_arg eq '-gc') {$GCislands = shift @ARGV;}  
 
    elsif ( $this_arg eq '-i') {$initialTable = shift @ARGV;}
    elsif ( $this_arg eq '-add') {$colomnesToAdd = shift @ARGV;}  

    elsif ( $this_arg eq '-lp') {$enhRight = shift @ARGV;}
    elsif ( $this_arg eq '-rightp') {$immediateDownstream = shift @ARGV;}  
    elsif ( $this_arg eq '-enh') {$enhLeft = shift @ARGV;} 
    elsif ( $this_arg eq '-dg') {$kb5 = shift @ARGV;}      

    

    elsif ( $this_arg =~ m/^-/ ) { print "unknown flag: $this_arg\n";}
}

my $maxDistFromGene = max(abs($enhLeft),abs($kb5));

if ( $GenesFilename eq "" ){
    die "you should specify a file with genes \n";
}
if(( $RNApolFilename eq "")&&($H3K36Me3polFilename eq "")&&($TF1Filename eq "")&&($H3K9Me3polFilename eq "")){  
    die "you should specify at least one file with peaks\n";
}


#-----------read selected genes----------------
my %selectedGenes;
my %selectedGenesFoldChange;
if ( $SelectedGenesFilename ne "") {
    open (FILE, "<$SelectedGenesFilename ") or die "Cannot open file $SelectedGenesFilename !!!!: $!";
    while (<FILE>) {
        s/\R//g;
        my @a = split/\t/;
        $selectedGenes{$a[1]} = $a[3];
	$selectedGenesFoldChange{$a[1]} = $a[2];
	#print "gene:$a[1],reg:$selectedGenes{$a[1]},FC:$selectedGenesFoldChange{$a[1]}\n"; 
    }
    
    close FILE;
    print "\t\t$SelectedGenesFilename is read!\n" if ($verbose);
}

#-----------read genes with fluorescence---------
my %fluoGenes;
if ( $fluoFile ne "") {
    open (FILE, "<$fluoFile ") or die "Cannot open file $fluoFile !!!!: $!";
    my $gene = "";
    my $med = 0;
    my %h;
    while (<FILE>) {
        s/\R//g;
        my @a = split/\t/;     

        next if (scalar(@a)<5);
        next if ($a[0] eq "probesets");
        next unless ($a[0] =~m/\S/);
        next unless ($a[4] =~m/\S/);
        if ($gene ne "" && $gene ne $a[2]) {
        	#calcMed        	
        	$med = med(keys %h);        	        	
        	$fluoGenes{$gene} = $med;
        	$med=0;
        	delete @h{keys %h};
        } else {
        	#$h{$a[4]} = 1;
        }
        $gene = $a[2];
        $h{$a[4]} = 1;        
       	#print "keys ", scalar(keys %h),"\t",keys %h,"\n";
    }
    if ($gene ne "") {
        $med = med(keys %h);
        $fluoGenes{$gene} = $med;
    }
    close FILE;
    print "\t\t$fluoFile is read!\n" if ($verbose);;
}
#-----------read GC-islands----------------
my %GCislands;
if ($GCislands ne "") {
    open (FILE, "<$GCislands ") or die "Cannot open file $GCislands !!!!: $!";
    
     while (<FILE>) {
        s/\R//g;
        my @a = split/\t/;
        #bin	chrom	chromStart	chromEnd	name	length	cpgNum	gcNum	perCpg	perGc	obsExp
        #107	chr1	36568608	36569851	CpG: 128	1243	128	766	20.6	61.6	1.09
        my $chr = $a[1];
        my $start = $a[2];
        my $end = $a[3];
        $GCislands{$chr}->{$start}=$end;
     }
    close FILE;
    if ($verbose) {
        print "$GCislands is read\n";
    }
} elsif ($verbose) {
    print "you did not specify a file with GC-islands\n";    
}

#-----------read genes----------------

my %genes;

my $count = 0;

open (GENES, "<$GenesFilename") or die "Cannot open file $GenesFilename!!!!: $!";
<GENES>;
while (<GENES>) {
    s/\R//g;	
    if (/(chr.*)\s([+-])\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\S+)\s(\S+)\s\S+\s(\S+)/){
        my $name = $10;
        my $chr = $1;
        my $strand = $2;
        if ($strand eq '+') {
                $strand = 1;
        }
        else {
                $strand = -1;	
        }
        my $leftPos = $3;
        my $rightPos = $4;
        my $cdsStart= $5;
        my $cdsEnd= $6;
        my $exonCount= $7;
        my $exonStarts= $8;
        my $exonEnds= $9;
        my $ID = "$name\t$chr:$leftPos-$rightPos\t$count";
        my $foldChange = 1;
        my $reg = "NA";
        my $fluo = "NA";
        if ( $SelectedGenesFilename ne "") {
                #print "$name\t";
                if (exists($selectedGenes{$name})) {
                        $reg = $selectedGenes{$name};                        
                        $foldChange = $selectedGenesFoldChange{$name};
                }
        }
        if ( $fluoFile ne "") {                
                if (exists($fluoGenes{$name})) {
                         $fluo = $fluoGenes{$name};
                }
        }        
        unless (exists($genes{$chr})) {
                my %h;
                $genes{$chr} = \%h;
        }
        
        my $RNAlength = 0;
        my $skip = 0;
        
        #print "$ID\n";
        if($longest) {
            $RNAlength = getRNAlength($exonStarts,$exonEnds);
            for my $IDgene (keys %{$genes{$chr}}) {
                my $nameGene= (split('\t', $IDgene))[0];               
                if ($nameGene eq $name && $RNAlength > $genes{$chr}->{$IDgene}{'RNAlength'}) {
                    #print  "found longer isofome: $ID longer than $IDgene\n";
                    #                    print "$RNAlength > ".$genes{$chr}->{$IDgene}{'RNAlength'}."\n";

                    $ID=$IDgene;
                } elsif ($nameGene eq $name && $RNAlength <= $genes{$chr}->{$IDgene}{'RNAlength'}) {
                    #print  "found shorter isofome: $ID shorted than $IDgene\nwill skip it\n";
                    #print "$RNAlength <= ".$genes{$chr}->{$IDgene}{'RNAlength'}."\n";
                    $skip = 1;
                }          
            }            
        }
        
        
        unless ($skip) {    
        
            unless (exists($genes{$chr}->{$ID})) {
                    my %h1;
                    $genes{$chr}->{$ID} = \%h1;	
                    $count++;
            }
            
            $genes{$chr}->{$ID}{'name'} = $name ;
            $genes{$chr}->{$ID}{'left'} = $leftPos ;
            $genes{$chr}->{$ID}{'right'} = $rightPos ;
            $genes{$chr}->{$ID}{'cdsStart'} = $cdsStart;
            $genes{$chr}->{$ID}{'cdsEnd'} = $cdsEnd;			
            $genes{$chr}->{$ID}{'strand'} = $strand; 
            $genes{$chr}->{$ID}{'exonCount'} = $exonCount;			
            $genes{$chr}->{$ID}{'exonStarts'} = $exonStarts; 
            $genes{$chr}->{$ID}{'exonEnds'} = $exonEnds; 			
            $genes{$chr}->{$ID}{'TSS'} = ($strand == 1) ? $leftPos :$rightPos ;
            $genes{$chr}->{$ID}{'TE'} = ($strand == -1) ? $leftPos :$rightPos ;		
            $genes{$chr}->{$ID}{'reg'} = $reg;
            $genes{$chr}->{$ID}{'foldChange'} = $foldChange;
            $genes{$chr}->{$ID}{'length'} = abs ($leftPos-$rightPos); 
            $genes{$chr}->{$ID}{'RNAlength'} = $RNAlength ;

     	    my $enhStart = min($genes{$chr}->{$ID}{'TSS'}-$enhLeft*$genes{$chr}->{$ID}{'strand'},$genes{$chr}->{$ID}{'TSS'}-$enhRight*$genes{$chr}->{$ID}{'strand'});
     	    my $enhEnd = max($genes{$chr}->{$ID}{'TSS'}-$enhLeft*$genes{$chr}->{$ID}{'strand'},$genes{$chr}->{$ID}{'TSS'}-$enhRight*$genes{$chr}->{$ID}{'strand'});
            $genes{$chr}->{$ID}{'enhStart'} =$enhStart;   
            $genes{$chr}->{$ID}{'enhEnd'} =$enhEnd;   

	    my $IntraStart = min($genes{$chr}->{$ID}{'TSS'}+$immediateDownstream*$genes{$chr}->{$ID}{'strand'},$genes{$chr}->{$ID}{'TE'});
     	    my $IntraEnd = max($genes{$chr}->{$ID}{'TSS'}+$immediateDownstream*$genes{$chr}->{$ID}{'strand'},$genes{$chr}->{$ID}{'TE'});
            $genes{$chr}->{$ID}{'IntraStart'} =$IntraStart;   
            $genes{$chr}->{$ID}{'IntraEnd'} =$IntraEnd;   

            $genes{$chr}->{$ID}{'fluo'} = $fluo;
          
            $genes{$chr}->{$ID}{'RNApolScore'} = 0;
            $genes{$chr}->{$ID}{'RNApolDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'RNApol_junctionScore'} = 0;
            $genes{$chr}->{$ID}{'RNApol_junctionDist'} = $INFINITY;            
            
            $genes{$chr}->{$ID}{'normalizedGBodyScore'} = 0;
            $genes{$chr}->{$ID}{'K9promScore'} = 0;
            $genes{$chr}->{$ID}{'K9promDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'K9largeScore'} = 0;
            $genes{$chr}->{$ID}{'K9largeDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'largePromScore'} = 0;
            $genes{$chr}->{$ID}{'promDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'EnhScore'} = 0;
            $genes{$chr}->{$ID}{'EnhDistTSS'} = $INFINITY;
            $genes{$chr}->{$ID}{'intraScore'} = 0;
            $genes{$chr}->{$ID}{'intraDistTSS'} = $INFINITY;
            $genes{$chr}->{$ID}{'allScore'} = 0;
            $genes{$chr}->{$ID}{'allDistTSS'} = $INFINITY;
            
            
            $genes{$chr}->{$ID}{'TFFirstIntronAndIntraScore'} = 0;
            $genes{$chr}->{$ID}{'TFFirstIntronAndIntraDist'} = $INFINITY;        
            $genes{$chr}->{$ID}{'firstIntronScore'} = 0;
            $genes{$chr}->{$ID}{'firstIntronDistTSS'} = $INFINITY;
            $genes{$chr}->{$ID}{'TFintraMinusFirstIntronScore'} = 0;
            $genes{$chr}->{$ID}{'TFintraMinusFirstIntronDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'immDownScore'} = 0;
            $genes{$chr}->{$ID}{'immDownDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'promSimpleScore'} = 0;
            $genes{$chr}->{$ID}{'promSimpleDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'TFenh60kbScore'} = 0;
            $genes{$chr}->{$ID}{'TFenh60kbDist'} = $INFINITY;
            
            $genes{$chr}->{$ID}{'TF_FirstExonScore'} = 0;
            $genes{$chr}->{$ID}{'TF_FirstExonDist'} = $INFINITY;            
            $genes{$chr}->{$ID}{'TF_junctionScore'} = 0;
            $genes{$chr}->{$ID}{'TF_junctionDist'} = $INFINITY;
            
            $genes{$chr}->{$ID}{'TF_junctionAndIntraScore'} = 0;
            $genes{$chr}->{$ID}{'TF_junctionAndIntraDist'} = $INFINITY;            
            $genes{$chr}->{$ID}{'TF_FirstExonAndIntraScore'} = 0;
            $genes{$chr}->{$ID}{'TF_FirstExonAndIntraDist'} = $INFINITY;            
            
            $genes{$chr}->{$ID}{'TF_OtherExonsScore'} = 0;
            $genes{$chr}->{$ID}{'TF_OtherExonsDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'TF_OtherExonsAndIntraScore'} = 0;
            $genes{$chr}->{$ID}{'TF_OtherExonsAndIntraDist'} = $INFINITY;

            
            $genes{$chr}->{$ID}{'TF_OtherIntronsScore'} = 0;
            $genes{$chr}->{$ID}{'TF_OtherIntronsDist'} = $INFINITY;   
            $genes{$chr}->{$ID}{'TF_OtherIntronsAndIntraScore'} = 0;
            $genes{$chr}->{$ID}{'TF_OtherIntronsAndIntraDist'} = $INFINITY;            
            
            $genes{$chr}->{$ID}{'geneDownstreamScore'} = 0;
            $genes{$chr}->{$ID}{'geneDownstreamDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'K9enhScore'} = 0;
            $genes{$chr}->{$ID}{'K9EnhDistTSS'} = $INFINITY;   
            

            $genes{$chr}->{$ID}{'score_GeneBody'} = 0;
            $genes{$chr}->{$ID}{'GeneBodyPeak_DistToTSS'} = $INFINITY;
            
            $genes{$chr}->{$ID}{'GCisland'} = 0;
            if ($GCislands ne "") {
                $genes{$chr}->{$ID}{'GCisland'} = checkIfGC ($genes{$chr}->{$ID}{'TSS'},$strand,2000,$GCislands{$chr});
            }
        }
    }
}

print "Total genes (including isoforms) : $count\n" ;
close GENES;
print "\t\t$GenesFilename is read!\n" if ($verbose);;
#for my $gName (sort keys %{$genes{'chr18'}}) {

 #    print "$gName\t$genes{'chr18'}->{$gName}{'TSS'}\n";
#}


#-----------read file with sites of TF1-----
my $numberOfAllSites = 0;

if ($TF1Filename eq "") {
    print "No file with peaks of TF1!\n" if ($verbose) ;
} else {
    open (FILE, "<$TF1Filename ") or die "Cannot open file $TF1Filename !!!!: $!";
    $_ = <FILE>;
    my $correction = 0;
    my @a = split /\t/;
    if (scalar(@a)>1 && $a[1] =~ m/chr/ ) {
            $correction = 1;		
    }

    while (<FILE>) {
        s/\R//g;		
            
        my @a = split /\t/;	
            
        my $chr = $a[0+$correction];
        my $firstPos = $a[1+$correction];
        my $LastPos = $a[2+$correction];

	my $maxPos = int(($firstPos+$LastPos)/2);
	if (scalar(@a)>(3+$correction)) {
		$maxPos = $a[3+$correction];
	}	
	if ($maxPos =~ m/\D/) {
		$maxPos = int(($firstPos+$LastPos)/2);
	}
	my $score = 0;
	if (scalar(@a)>(4+$correction)) {
		$score = $a[4+$correction];
	}
        for my $ID (keys %{$genes{$chr}}) {
	    my $distTSS;
	    if ($firstPos<=$genes{$chr}->{$ID}{'TSS'} && $LastPos>=$genes{$chr}->{$ID}{'TSS'}) {
		 $distTSS = 0 ;
	    } else {
		$distTSS = ($firstPos - $genes{$chr}->{$ID}{'TSS'})*$genes{$chr}->{$ID}{'strand'};
	    	$distTSS = ($LastPos - $genes{$chr}->{$ID}{'TSS'})*$genes{$chr}->{$ID}{'strand'} if (abs($LastPos - $genes{$chr}->{$ID}{'TSS'})<abs($distTSS));
	    }

            my $distTE;
	    if ($firstPos<=$genes{$chr}->{$ID}{'TE'} && $LastPos>=$genes{$chr}->{$ID}{'TE'}) {
 		$distTE = 0;
	    } else {
               $distTE = ($firstPos - $genes{$chr}->{$ID}{'TE'})*$genes{$chr}->{$ID}{'strand'};
	       $distTE = ($LastPos - $genes{$chr}->{$ID}{'TE'})*$genes{$chr}->{$ID}{'strand'} if (abs($LastPos - $genes{$chr}->{$ID}{'TE'})<abs($distTE));
	    }

            if (($distTSS>= $enhLeft)&&($distTE<=$kb5)) {
                if ($genes{$chr}->{$ID}{'allScore'}<$score) {
                    $genes{$chr}->{$ID}{'allScore'}=$score;
                    $genes{$chr}->{$ID}{'allDistTSS'} = $distTSS;
                }   
            } else {next;}	 

	    my $hmCorDist = -1;
            if ($LastPos>=$genes{$chr}->{$ID}{'left'} && $firstPos <= $genes{$chr}->{$ID}{'right'} ) {
                  $hmCorDist = 0;            		            	
            }                         

            if ($hmCorDist == 0) {
			if ($genes{$chr}->{$ID}{'score_GeneBody'}<$score) {
		           $genes{$chr}->{$ID}{'score_GeneBody'}=$score;
		           $genes{$chr}->{$ID}{'GeneBodyPeak_DistToTSS'} = $distTSS;
                	}		
				
            
			my $scoreToadd = 0;

			if (($firstPos>=$genes{$chr}->{$ID}{'left'})&&($LastPos<=$genes{$chr}->{$ID}{'right'})) {
				$scoreToadd = $score/2.*($LastPos-$firstPos+1)/($genes{$chr}->{$ID}{'right'}-$genes{$chr}->{$ID}{'left'}+1);
			} elsif (($firstPos>=$genes{$chr}->{$ID}{'left'})&&($firstPos<$genes{$chr}->{$ID}{'right'})&&($LastPos>$genes{$chr}->{$ID}{'right'})) {
				    $scoreToadd = $score/2.*($genes{$chr}->{$ID}{'right'}-$firstPos+1)/($genes{$chr}->{$ID}{'right'}-$genes{$chr}->{$ID}{'left'}+1);
			}
			elsif (($firstPos<$genes{$chr}->{$ID}{'left'})&&($LastPos>$genes{$chr}->{$ID}{'left'})&&($LastPos<=$genes{$chr}->{$ID}{'right'})) {
				    $scoreToadd = $score/2.*($LastPos-$genes{$chr}->{$ID}{'left'}+1)/($genes{$chr}->{$ID}{'right'}-$genes{$chr}->{$ID}{'left'}+1);
			}
			elsif (($firstPos<$genes{$chr}->{$ID}{'left'})&&($LastPos>$genes{$chr}->{$ID}{'right'})) {
				    $scoreToadd = $score/2.;
			}			    		
			$genes{$chr}->{$ID}{'normalizedGBodyScore'} += $scoreToadd;    
	    } 	    

            if ( $LastPos>=$genes{$chr}->{$ID}{'enhStart'} && $firstPos <=$genes{$chr}->{$ID}{'enhEnd'} ) {  
                if ($genes{$chr}->{$ID}{'EnhScore'}<$score) {
                    $genes{$chr}->{$ID}{'EnhScore'}=$score;
                    $genes{$chr}->{$ID}{'EnhDistTSS'} = $distTSS;
                }
            } 
	    if (($distTSS>= $enhRight)&&($distTSS<=$immediateDownstream)) {
                if ($genes{$chr}->{$ID}{'largePromScore'}<$score) {
                    $genes{$chr}->{$ID}{'largePromScore'}=$score;
                    $genes{$chr}->{$ID}{'promDist'} = $distTSS;
                }
            } 
            if ( $LastPos>=$genes{$chr}->{$ID}{'IntraStart'} && $firstPos <= $genes{$chr}->{$ID}{'IntraEnd'} ) { 
                if ($genes{$chr}->{$ID}{'intraScore'}<$score) {
                    $genes{$chr}->{$ID}{'intraScore'}=$score;
                    $genes{$chr}->{$ID}{'intraDistTSS'} = $distTSS;
                }
            }

            
            if (($distTSS>= 0)&&($distTSS<=$immediateDownstream)) {
                if ($genes{$chr}->{$ID}{'immDownScore'}<$score) {
                    $genes{$chr}->{$ID}{'immDownScore'}=$score;
                    $genes{$chr}->{$ID}{'immDownDist'} = $distTSS;
                }
            }
            if (($distTSS<= 0)&&($distTSS>=$enhRight)) {
                if ($genes{$chr}->{$ID}{'promSimpleScore'}<$score) {
                    $genes{$chr}->{$ID}{'promSimpleScore'}=$score;
                    $genes{$chr}->{$ID}{'promSimpleDist'} = $distTSS;
                }
            }             

            if (($distTE>=0)&&($distTE<=$kb5)) {
                if ($genes{$chr}->{$ID}{'geneDownstreamScore'}<$score) {
                    $genes{$chr}->{$ID}{'geneDownstreamScore'}=$score;
                    $genes{$chr}->{$ID}{'geneDownstreamDist'} = $distTSS;
                }
            }  
            
        }
        $numberOfAllSites++;
    }

    close FILE;
    print "\t$TF1Filename is read!\n" if ($verbose) ;
    print "$numberOfAllSites sites\n" ;
}




open (OUT , ">$outname") or die "Cannot open file $outname!!!!: $!";

print OUT "name\tchr\tstart\tend\tstrand\tReg\tfoldChange\t";

if ($GCislands ne "") {
   print OUT "GC-island\t";
}

if ( $fluoFile ne "") {
   print OUT "fluorescence\t";
}

print OUT "GeneBodyNormalizedScore\t";

if ($TF1Filename ne "") {
	print OUT "score_Gene\tdistTSS_Gene\t";
	print OUT "score_Promoter\tdistTSS_Promoter\t";
	print OUT "score_ImmDown\tdistTSS_ImmDown\t";
	print OUT "score_PromoterORImmDown\tdistTSS_PromoterORImmDown\t";
	print OUT "score_Enhancer\tdistTSS_Enhancer\t";
	print OUT "score_Intragenic\tdistTSS_Intragenic\t";
	print OUT "score_GeneDownstream\tdistTSS_GeneDownstream\t";


	print OUT "score_GeneBody\tGeneBodyPeak_DistToTSS"
}

print OUT "\n";

for my $chr (keys %genes) {
    for my $ID (keys %{$genes{$chr}}) {
        print OUT "$genes{$chr}->{$ID}{'name'}\t$chr\t$genes{$chr}->{$ID}{'left'}\t$genes{$chr}->{$ID}{'right'}\t$genes{$chr}->{$ID}{'strand'}\t$genes{$chr}->{$ID}{'reg'}\t$genes{$chr}->{$ID}{'foldChange'}\t";

	if ($GCislands ne "") {
	   print OUT "$genes{$chr}->{$ID}{'GCisland'}\t";
	}
	if ( $fluoFile ne "") {
           print OUT "$genes{$chr}->{$ID}{'fluo'}\t"; 
	}

	print OUT "$genes{$chr}->{$ID}{'normalizedGBodyScore'}\t";	
	
	if ($TF1Filename ne "") {  
		print OUT "$genes{$chr}->{$ID}{'allScore'}\t$genes{$chr}->{$ID}{'allDistTSS'}\t";
		print OUT "$genes{$chr}->{$ID}{'promSimpleScore'}\t$genes{$chr}->{$ID}{'promSimpleDist'}\t";
		print OUT "$genes{$chr}->{$ID}{'immDownScore'}\t$genes{$chr}->{$ID}{'immDownDist'}\t";
		print OUT "$genes{$chr}->{$ID}{'largePromScore'}\t$genes{$chr}->{$ID}{'promDist'}\t";
		print OUT "$genes{$chr}->{$ID}{'EnhScore'}\t$genes{$chr}->{$ID}{'EnhDistTSS'}\t";
		print OUT "$genes{$chr}->{$ID}{'intraScore'}\t$genes{$chr}->{$ID}{'intraDistTSS'}\t";
		print OUT "$genes{$chr}->{$ID}{'geneDownstreamScore'}\t$genes{$chr}->{$ID}{'geneDownstreamDist'}\t";

		print OUT "$genes{$chr}->{$ID}{'score_GeneBody'}\t$genes{$chr}->{$ID}{'GeneBodyPeak_DistToTSS'}"

        }
        print OUT "\n";
    }
}

close OUT;

###################################
sub med {
	my @arr = @_;
	my $med = 0;
	@arr = sort {$a <=> $b} @arr; 
	if ((scalar(@arr)/2) =~ m/[\.\,]5/) {
		return $arr[floor(scalar(@arr)/2)];
	} else {
		return ($arr[scalar(@arr)/2]+$arr[scalar(@arr)/2-1])/2;
	}
	$med;	
}

sub checkIfGC {
    my ($TSS,$strand,$dist,$GCislandsChr)=@_; 
    my $ifGC = 0;
    my $leftProm=$TSS-$dist;
    my $rightProm = $TSS;    
    if ($strand== -1) {
	my $leftProm=$TSS;
	my $rightProm = $TSS+$dist;	
    } #print "$leftProm\t"; print "$rightProm\n";
    for my $leftGC (keys %{$GCislandsChr}) {
	my $rightGC = $GCislandsChr->{$leftGC}; 
	if ($leftGC>=$leftProm&&$leftGC<=$rightProm || $rightGC>=$leftProm&&$rightGC<=$rightProm) {
	    return "GC-island";
	}
    }       
    return $ifGC ;
}

sub getFirstIntron {
	my ($exonCount,$exonStarts,$exonEnds,$strand) = @_;
	my ($left,$right);
	if ($exonCount == 1) {
		return (0,0);
	}
	if ($strand == 1) {
		$left = (split ",", $exonEnds)[0];
		$right = (split (",", $exonStarts))[1];
	} else {
		$left = (split (",", $exonEnds))[$exonCount-2];
		$right = (split (",", $exonStarts))[$exonCount-1];	
	}
	($left,$right);
}

sub getFirstExon {
	my ($exonCount,$exonStarts,$exonEnds,$strand) = @_;
	my ($left,$right);
	if ($exonCount == 1) {
		return (0,0);
	}
	if ($strand == 1) {
		$left = (split ",", $exonStarts)[0];
		$right = (split (",", $exonEnds))[0]-$jonctionSize;
	} else {
		$left = (split (",", $exonStarts))[$exonCount-1]+$jonctionSize;
		$right = (split (",", $exonEnds))[$exonCount-1];	
	}
	($left,$right);
}


sub getIntronExon {
	my ($pos,$exonCount,$exonStarts,$exonEnds,$strand) = @_;
	my (@left,@right);
	@left = (split ",", $exonStarts);
	@right = (split (",", $exonEnds));

	for (my $i = 0; $i<$exonCount;$i++) {
		#print "$left[$i] <= $pos ? $pos <= $right[$i]\n";
		if (($left[$i]+$jonctionSize < $pos) && ($pos < $right[$i]-$jonctionSize)) {
			#print "URA!\n";
			return "exon";
		} elsif (($i+1<$exonCount)&&($right[$i]+$jonctionSize < $pos) && ($pos < $left[$i+1]-$jonctionSize)) {
			return "intron";
		}
	}
	return "jonction";	
}


sub getTypeIntra {
    
    my ($geneEntry, $pos) = @_;
    my $type;
    
    if (($pos >= $geneEntry->{'firstIntronStart'})&&($pos <=$geneEntry->{'firstIntronEnd'})) {	   
	    return "f_intron";
    }
    if (($pos >= $geneEntry->{'firstExonStart'})&&($pos <=$geneEntry->{'firstExonEnd'})) {	   
	    return "f_exon";
    }
    $type = getIntronExon ($pos, $geneEntry->{'exonCount'},$geneEntry->{'exonStarts'},$geneEntry->{'exonEnds'},$geneEntry->{'strand'});
    return $type;
}

sub getRNAlength {
    my ($exonStarts,$exonEnds) = @_;
    my (@left,@right);
    @left = (split ",", $exonStarts);
    @right = (split (",", $exonEnds));
    my $length = 0;
    for (my $i = 0; $i<scalar(@right);$i++) {
	$length+=$right[$i]-$left[$i];
    }
    #print STDERR "length = $length\n";
    return $length;
}

sub min {
	my @a = sort {$a <=> $b} @_; 
	$a[0];	
}

sub max {
	my @a = sort {$b <=> $a} @_; 
	$a[0];	
}
