#!/usr/bin/perl -w
#outputs statistics for all genes in the list

#different boundaries
#no motif p-value for binding sites
#read directly No-resp/Up/down category

#all isoforms from the file with genes

#RNApolII sites on junctions

#has Histone Modification Mode: covering TSS and overlaping genes

use strict;
use POSIX;

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
    -HMmode	1/0		for histone modification type H3K27me3: look for overlap with TSS or gene body (default 0)

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
my $longEnhLeft = -60000;
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
    elsif ( $this_arg eq '-HMmode') {$ifHMmode = shift @ARGV; if ($ifHMmode eq "no" || $ifHMmode eq "0") {$ifHMmode = 0;} else {$ifHMmode=1;}}      
    
    

    elsif ( $this_arg =~ m/^-/ ) { print "unknown flag: $this_arg\n";}
}


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
    
            $genes{$chr}->{$ID}{'fluo'} = $fluo;
          
            $genes{$chr}->{$ID}{'RNApolScore'} = 0;
            $genes{$chr}->{$ID}{'RNApolDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'RNApol_junctionScore'} = 0;
            $genes{$chr}->{$ID}{'RNApol_junctionDist'} = $INFINITY;            
            
            $genes{$chr}->{$ID}{'K36score'} = 0;
            $genes{$chr}->{$ID}{'K9promScore'} = 0;
            $genes{$chr}->{$ID}{'K9promDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'K9largeScore'} = 0;
            $genes{$chr}->{$ID}{'K9largeDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'TFpromScore'} = 0;
            $genes{$chr}->{$ID}{'TFpromDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'TFenhScore'} = 0;
            $genes{$chr}->{$ID}{'TFenhDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'TFintraScore'} = 0;
            $genes{$chr}->{$ID}{'TFintraDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'TFallScore'} = 0;
            $genes{$chr}->{$ID}{'TFallDist'} = $INFINITY;
            
            
            $genes{$chr}->{$ID}{'TFFirstIntronAndIntraScore'} = 0;
            $genes{$chr}->{$ID}{'TFFirstIntronAndIntraDist'} = $INFINITY;        
            $genes{$chr}->{$ID}{'TFFirstIntronScore'} = 0;
            $genes{$chr}->{$ID}{'TFFirstIntronDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'TFintraMinusFirstIntronScore'} = 0;
            $genes{$chr}->{$ID}{'TFintraMinusFirstIntronDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'TFImmDownScore'} = 0;
            $genes{$chr}->{$ID}{'TFImmDownDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'TFpromSimpleScore'} = 0;
            $genes{$chr}->{$ID}{'TFpromSimpleDist'} = $INFINITY;
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
            
            $genes{$chr}->{$ID}{'TF5kbDownScore'} = 0;
            $genes{$chr}->{$ID}{'TF5kbDownDist'} = $INFINITY;
            $genes{$chr}->{$ID}{'K9enhScore'} = 0;
            $genes{$chr}->{$ID}{'K9enhDist'} = $INFINITY;   
            
            $genes{$chr}->{$ID}{'K27TSS_Score'} = 0;
            $genes{$chr}->{$ID}{'K27TSS_Dist'} = $INFINITY;
            $genes{$chr}->{$ID}{'K27GeneB_Score'} = 0;
            $genes{$chr}->{$ID}{'K27GeneB_Dist'} = $INFINITY;
            
            ($genes{$chr}->{$ID}{'firstIntronStart'},$genes{$chr}->{$ID}{'firstIntronEnd'}) = getFirstIntron ($exonCount,$exonStarts,$exonEnds,$strand);
            ($genes{$chr}->{$ID}{'firstExonStart'},$genes{$chr}->{$ID}{'firstExonEnd'}) = getFirstExon ($exonCount,$exonStarts,$exonEnds,$strand);
            
            $genes{$chr}->{$ID}{'exonCount'} = $exonCount;			
            $genes{$chr}->{$ID}{'exonStarts'} = $exonStarts; 
            $genes{$chr}->{$ID}{'exonEnds'} = $exonEnds; 
            
            
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

#-----------read file with sites miRNA, store as genes-----

if ( $MirFilename eq ""){
    print "you did not specify file with miRNA\n" if ($verbose);;
}
else {
    $count = 0;
    open (MIR, "<$MirFilename ") or die "Cannot open file $MirFilename !!!!: $!";
    #chr1	20669090	20669163	mmu-mir-206	960	+
    while (<MIR>) {
        s/\R//g;	
        my ($name, $chr, $leftPos, $rightPos, $strand );
#1	.	miRNA	20669091	20669163	.	+	.	ACC="MI0000249"; ID="mmu-mir-206";
        if (/([0-9XYM]+)\s.\smiRNA\s(\d+)\s(\d+)\s.\s([+-])\s.\sACC=.*ID=\"(.*)\"/) {
                $name = $5;
                $chr = $1;
                $leftPos = $2;
                $rightPos = $3;
                $strand = $4;
        }
        elsif (/(.*)\s(\d+)\s(\d+)\s(.*)\s(.*)\s(.*)/){
                $name = $4;
                $chr = $1;
                $leftPos = $2;
                $rightPos = $3;
                $strand = $6;
        } else {
                next;
        }

        unless ($chr =~ m/chr/) {
                $chr = "chr".$chr;
        }
        my $ID = "$name\t$chr:$leftPos-$rightPos\t$count";

        if ($strand eq '+') {
                $strand = 1;
        }
        else {
                $strand = -1;	
        }

        unless (exists($genes{$chr})) {
                my %h;
                $genes{$chr} = \%h;
        }
        unless (exists($genes{$chr}->{$ID})) {
                my %h1;
                $genes{$chr}->{$ID} = \%h1;	
                $count++;
        }
        $genes{$chr}->{$ID}{'name'} = $name ;
        $genes{$chr}->{$ID}{'left'} = $leftPos ;
        $genes{$chr}->{$ID}{'right'} = $rightPos ;
        $genes{$chr}->{$ID}{'cdsStart'} = $leftPos ;
        $genes{$chr}->{$ID}{'cdsEnd'} = $rightPos ;			
        $genes{$chr}->{$ID}{'strand'} = $strand;
        $genes{$chr}->{$ID}{'length'} = abs ($leftPos-$rightPos); 
        $genes{$chr}->{$ID}{'exonCount'} = 1;			
        $genes{$chr}->{$ID}{'exonStarts'} = $leftPos ; 
        $genes{$chr}->{$ID}{'exonEnds'} = $rightPos ; 			
        $genes{$chr}->{$ID}{'TSS'} = ($strand == 1) ? $leftPos :$rightPos ;
        $genes{$chr}->{$ID}{'TE'} = ($strand == -1) ? $leftPos :$rightPos ;		
        $genes{$chr}->{$ID}{'reg'} = "miRNA";
        $genes{$chr}->{$ID}{'foldChange'} = 1;	

        $genes{$chr}->{$ID}{'fluo'} = "N/A";


        $genes{$chr}->{$ID}{'RNApolScore'} = 0;
        $genes{$chr}->{$ID}{'RNApolDist'} = $INFINITY;
        $genes{$chr}->{$ID}{'RNApol_junctionScore'} = 0;
        $genes{$chr}->{$ID}{'RNApol_junctionDist'} = $INFINITY;           
        $genes{$chr}->{$ID}{'K36score'} = 0;
        $genes{$chr}->{$ID}{'K9promScore'} = 0;
        $genes{$chr}->{$ID}{'K9promDist'} = $INFINITY;
        $genes{$chr}->{$ID}{'K9largeScore'} = 0;
        $genes{$chr}->{$ID}{'K9largeDist'} = $INFINITY;
        $genes{$chr}->{$ID}{'TFpromScore'} = 0;
        $genes{$chr}->{$ID}{'TFpromDist'} = $INFINITY;
        $genes{$chr}->{$ID}{'TFenhScore'} = 0;
        $genes{$chr}->{$ID}{'TFenhDist'} = $INFINITY;
        $genes{$chr}->{$ID}{'TFintraScore'} = 0;
        $genes{$chr}->{$ID}{'TFintraDist'} = $INFINITY;
        $genes{$chr}->{$ID}{'TFallScore'} = 0;
        $genes{$chr}->{$ID}{'TFallDist'} = $INFINITY;
        
        $genes{$chr}->{$ID}{'TFFirstIntronAndIntraScore'} = 0;
        $genes{$chr}->{$ID}{'TFFirstIntronAndIntraDist'} = $INFINITY;        
        $genes{$chr}->{$ID}{'TFFirstIntronScore'} = 0;
        $genes{$chr}->{$ID}{'TFFirstIntronDist'} = $INFINITY;
        $genes{$chr}->{$ID}{'TFintraMinusFirstIntronScore'} = 0;
        $genes{$chr}->{$ID}{'TFintraMinusFirstIntronDist'} = $INFINITY;
        $genes{$chr}->{$ID}{'TFImmDownScore'} = 0;
        $genes{$chr}->{$ID}{'TFImmDownDist'} = $INFINITY;
        $genes{$chr}->{$ID}{'TFpromSimpleScore'} = 0;
        $genes{$chr}->{$ID}{'TFpromSimpleDist'} = $INFINITY;
        
        $genes{$chr}->{$ID}{'TF_FirstExonScore'} = 0;
        $genes{$chr}->{$ID}{'TF_FirstExonDist'} = $INFINITY;
        $genes{$chr}->{$ID}{'TF_OtherExonsScore'} = 0;
        $genes{$chr}->{$ID}{'TF_OtherExonsDist'} = $INFINITY;
        $genes{$chr}->{$ID}{'TF_OtherExonsAndIntraScore'} = 0;
        $genes{$chr}->{$ID}{'TF_OtherExonsAndIntraDist'} = $INFINITY;
        $genes{$chr}->{$ID}{'TF_junctionScore'} = 0;
        $genes{$chr}->{$ID}{'TF_junctionDist'} = $INFINITY;
        $genes{$chr}->{$ID}{'TF_OtherIntronsScore'} = 0;
        $genes{$chr}->{$ID}{'TF_OtherIntronsDist'} = $INFINITY;    
        $genes{$chr}->{$ID}{'TF_OtherIntronsAndIntraScore'} = 0;
        $genes{$chr}->{$ID}{'TF_OtherIntronsAndIntraDist'} = $INFINITY;
        
        $genes{$chr}->{$ID}{'TF_junctionAndIntraScore'} = 0;
        $genes{$chr}->{$ID}{'TF_junctionAndIntraDist'} = $INFINITY;            
        $genes{$chr}->{$ID}{'TF_FirstExonAndIntraScore'} = 0;
        $genes{$chr}->{$ID}{'TF_FirstExonAndIntraDist'} = $INFINITY;        
        
        $genes{$chr}->{$ID}{'TFenh60kbScore'} = 0;
        $genes{$chr}->{$ID}{'TFenh60kbDist'} = $INFINITY;
        
        $genes{$chr}->{$ID}{'TF5kbDownScore'} = 0;
        $genes{$chr}->{$ID}{'TF5kbDownDist'} = $INFINITY;
        $genes{$chr}->{$ID}{'K9enhScore'} = 0;
        $genes{$chr}->{$ID}{'K9enhDist'} = $INFINITY;  
        
        $genes{$chr}->{$ID}{'K27TSS_Score'} = 0;
	$genes{$chr}->{$ID}{'K27TSS_Dist'} = $INFINITY;
	$genes{$chr}->{$ID}{'K27GeneB_Score'} = 0;
        $genes{$chr}->{$ID}{'K27GeneB_Dist'} = $INFINITY;
        
        ($genes{$chr}->{$ID}{'firstIntronStart'},$genes{$chr}->{$ID}{'firstIntronEnd'}) = (0,0);
        ($genes{$chr}->{$ID}{'firstExonStart'},$genes{$chr}->{$ID}{'firstExonEnd'}) = (0,0);

        $genes{$chr}->{$ID}{'GCisland'} = 0;
        
        $genes{$chr}->{$ID}{'exonCount'} = 1;			
	$genes{$chr}->{$ID}{'exonStarts'} = $leftPos ; 
	$genes{$chr}->{$ID}{'exonEnds'} = $rightPos ; 	

            
    }


    close MIR;
    print "\t\t$MirFilename is read!\n" if ($verbose) ;
    print "$count miRNA\n" if ($verbose);;
}

#-----------read file with sites of TF1-----
my $numberOfAllSites = 0;

if ($TF1Filename eq "") {
    print "No file with peaks of TF1!\n" if ($verbose) ;
} else {
    open (FILE, "<$TF1Filename ") or die "Cannot open file $TF1Filename !!!!: $!";
    $_ = <FILE>;
    my $correction = 0;
    my @a = split /\t/;
    if ( $a[1] =~ m/chr/ ) {
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
	if ($maxPos=~/\D/) {
		$maxPos = int(($firstPos+$LastPos)/2);
	}
	my $score = 0;
	if (scalar(@a)>(4+$correction)) {
		$score = $a[4+$correction];
	}
        for my $ID (keys %{$genes{$chr}}) {

            my $distTSS = ($maxPos - $genes{$chr}->{$ID}{'TSS'})*$genes{$chr}->{$ID}{'strand'};
            my $distTE = ($maxPos - $genes{$chr}->{$ID}{'TE'})*$genes{$chr}->{$ID}{'strand'};
            
            my $hmDist;
            my $hmCorDist = -1;
            if ($ifHMmode) {
            	if ($genes{$chr}->{$ID}{'strand'}==1) {
            		if ($LastPos <= $genes{$chr}->{$ID}{'TSS'}) {
            			$hmDist = $LastPos-$genes{$chr}->{$ID}{'TSS'};
            		} elsif ($firstPos >= $genes{$chr}->{$ID}{'TSS'}) {
            			$hmDist = $firstPos - $genes{$chr}->{$ID}{'TSS'};
            		} else {$hmDist = 0;}
            		
            		if ($hmDist >= 0 && $genes{$chr}->{$ID}{'TE'} >= $firstPos) {
            			$hmCorDist = 0;
            		}           		

            	} else {            	
            	
            		if ($LastPos <= $genes{$chr}->{$ID}{'TSS'}) {
		            	$hmDist = $genes{$chr}->{$ID}{'TSS'}- $LastPos;
		        } elsif ($firstPos >= $genes{$chr}->{$ID}{'TSS'}) {
		            	$hmDist = $genes{$chr}->{$ID}{'TSS'}-$firstPos;
		        } else {$hmDist = 0;}
		            		
		        if ($hmDist >= 0 && $genes{$chr}->{$ID}{'TE'} <= $LastPos) {
		            	$hmCorDist = 0;
            		}            	
            	} 
            	
            	
            	if (($hmDist>= $enhRight)&&($hmDist<=$immediateDownstream)) {
		        if ($genes{$chr}->{$ID}{'K27TSS_Score'}<$score) {
		           $genes{$chr}->{$ID}{'K27TSS_Score'}=$score;
		           $genes{$chr}->{$ID}{'K27TSS_Dist'} = $hmDist;
                	}
                }
                
                if ($hmDist>$immediateDownstream && $hmCorDist == 0) {
			if ($genes{$chr}->{$ID}{'K27GeneB_Score'}<$score) {
		           $genes{$chr}->{$ID}{'K27GeneB_Score'}=$score;
		           $genes{$chr}->{$ID}{'K27GeneB_Dist'} = $hmDist;
                	}                
                }
            	
            }
            
            
            
            if (($distTSS>= $enhLeft)&&($distTSS<$enhRight)) {
                if ($genes{$chr}->{$ID}{'TFenhScore'}<$score) {
                    $genes{$chr}->{$ID}{'TFenhScore'}=$score;
                    $genes{$chr}->{$ID}{'TFenhDist'} = $distTSS;
                }
            } elsif (($distTSS>= $enhRight)&&($distTSS<=$immediateDownstream)) {
                if ($genes{$chr}->{$ID}{'TFpromScore'}<$score) {
                    $genes{$chr}->{$ID}{'TFpromScore'}=$score;
                    $genes{$chr}->{$ID}{'TFpromDist'} = $distTSS;
                }
            } elsif (($distTSS >= $immediateDownstream)&&($distTE<=0)) {
                if ($genes{$chr}->{$ID}{'TFintraScore'}<$score) {
                    $genes{$chr}->{$ID}{'TFintraScore'}=$score;
                    $genes{$chr}->{$ID}{'TFintraDist'} = $distTSS;
                }
            }
            if (($distTSS>= $enhLeft)&&($distTE<=$kb5)) {
                if ($genes{$chr}->{$ID}{'TFallScore'}<$score) {
                    $genes{$chr}->{$ID}{'TFallScore'}=$score;
                    $genes{$chr}->{$ID}{'TFallDist'} = $distTSS;
                }   
            }
            
            if (($distTSS>= 0)&&($distTSS<=$immediateDownstream)) {
                if ($genes{$chr}->{$ID}{'TFImmDownScore'}<$score) {
                    $genes{$chr}->{$ID}{'TFImmDownScore'}=$score;
                    $genes{$chr}->{$ID}{'TFImmDownDist'} = $distTSS;
                }
            }
            if (($distTSS<= 0)&&($distTSS>=$enhRight)) {
                if ($genes{$chr}->{$ID}{'TFpromSimpleScore'}<$score) {
                    $genes{$chr}->{$ID}{'TFpromSimpleScore'}=$score;
                    $genes{$chr}->{$ID}{'TFpromSimpleDist'} = $distTSS;
                }
            }     
            
            my ($firstIntronStart,$firstIntronEnd)=($genes{$chr}->{$ID}{'firstIntronStart'},$genes{$chr}->{$ID}{'firstIntronEnd'});
            ($firstIntronStart,$firstIntronEnd)= ($firstIntronEnd,$firstIntronStart) if ($firstIntronStart>$firstIntronEnd) ;
            
            my ($firstExonStart,$firstExonEnd) = ($genes{$chr}->{$ID}{'firstExonStart'},$genes{$chr}->{$ID}{'firstExonEnd'}) ;
            ($firstExonStart,$firstExonEnd)= ($firstExonEnd,$firstExonStart) if ($firstExonStart>$firstExonEnd) ;
            
            if ($maxPos>=$firstIntronStart && $maxPos <= $firstIntronEnd) {
                if ($genes{$chr}->{$ID}{'TFFirstIntronScore'}<$score) {
                    $genes{$chr}->{$ID}{'TFFirstIntronScore'}=$score;
                    $genes{$chr}->{$ID}{'TFFirstIntronDist'} = $distTSS;
                }
                
                if (($distTSS >= $immediateDownstream)&&($distTE<=0)) {
                    if ($genes{$chr}->{$ID}{'TFFirstIntronAndIntraScore'}<$score) {
                        $genes{$chr}->{$ID}{'TFFirstIntronAndIntraScore'}=$score;
                        $genes{$chr}->{$ID}{'TFFirstIntronAndIntraDist'} = $distTSS;
                    }
                }
            } elsif (($distTSS >= $immediateDownstream)&&($distTE<=0)) {
                if ($genes{$chr}->{$ID}{'TFintraMinusFirstIntronScore'}<$score) {
                    $genes{$chr}->{$ID}{'TFintraMinusFirstIntronScore'}=$score;
                    $genes{$chr}->{$ID}{'TFintraMinusFirstIntronDist'} = $distTSS;
                }
            }     
            if (($distTSS>= $longEnhLeft)&&($distTSS<$enhRight)) {
                if ($genes{$chr}->{$ID}{'TFenh60kbScore'}<$score) {
                    $genes{$chr}->{$ID}{'TFenh60kbScore'}=$score;
                    $genes{$chr}->{$ID}{'TFenh60kbDist'} = $distTSS;
                }
            }
            if (($distTE>=0)&&($distTE<=$kb5)) {
                if ($genes{$chr}->{$ID}{'TF5kbDownScore'}<$score) {
                    $genes{$chr}->{$ID}{'TF5kbDownScore'}=$score;
                    $genes{$chr}->{$ID}{'TF5kbDownDist'} = $distTSS;
                }
            }
            if ($distTSS>=0 && $distTE<=0) {
                my $typeIntra = &getTypeIntra($genes{$chr}->{$ID},$maxPos);
                if ($typeIntra eq "f_exon") {
                    if ($genes{$chr}->{$ID}{'TF_FirstExonScore'}<$score) {
                        $genes{$chr}->{$ID}{'TF_FirstExonScore'}=$score;
                        $genes{$chr}->{$ID}{'TF_FirstExonDist'} = $distTSS;
                    }
                    
                    if ($genes{$chr}->{$ID}{'TF_FirstExonAndIntraScore'}<$score && ($distTSS >= $immediateDownstream)&&($distTE<=0)) {
                        $genes{$chr}->{$ID}{'TF_FirstExonAndIntraScore'}=$score;
                        $genes{$chr}->{$ID}{'TF_FirstExonAndIntraDist'} = $distTSS;
                    }   
                    
                } else {
                    if ($typeIntra eq "exon") {
                        if ($genes{$chr}->{$ID}{'TF_OtherExonsScore'}<$score) {
                            $genes{$chr}->{$ID}{'TF_OtherExonsScore'}=$score;
                            $genes{$chr}->{$ID}{'TF_OtherExonsDist'} = $distTSS;
                        }
                        
                        if (($distTSS >= $immediateDownstream)&&($distTE<=0)) {
                            if ($genes{$chr}->{$ID}{'TF_OtherExonsAndIntraScore'}<$score) {
                                $genes{$chr}->{$ID}{'TF_OtherExonsAndIntraScore'}=$score;
                                $genes{$chr}->{$ID}{'TF_OtherExonsAndIntraDist'} = $distTSS;
                            }
                        }
                    
                    } elsif ($typeIntra eq "intron") {
                        if ($genes{$chr}->{$ID}{'TF_OtherIntronsScore'}<$score) {
                            $genes{$chr}->{$ID}{'TF_OtherIntronsScore'}=$score;
                            $genes{$chr}->{$ID}{'TF_OtherIntronsDist'} = $distTSS;
                        }
                        
                        if (($distTSS >= $immediateDownstream)&&($distTE<=0)) {
                            if ($genes{$chr}->{$ID}{'TF_OtherIntronsAndIntraScore'}<$score) {
                                $genes{$chr}->{$ID}{'TF_OtherIntronsAndIntraScore'}=$score;
                                $genes{$chr}->{$ID}{'TF_OtherIntronsAndIntraDist'} = $distTSS;
                            }
                        }                    
                    } elsif ($typeIntra eq "jonction") {
                        if ($genes{$chr}->{$ID}{'TF_junctionScore'}<$score) {
                            $genes{$chr}->{$ID}{'TF_junctionScore'}=$score;
                            $genes{$chr}->{$ID}{'TF_junctionDist'} = $distTSS;
                        }
                        
                        if ($genes{$chr}->{$ID}{'TF_junctionAndIntraScore'}<$score && ($distTSS >= $immediateDownstream)&&($distTE<=0)) {
                            $genes{$chr}->{$ID}{'TF_junctionAndIntraScore'}=$score;
                            $genes{$chr}->{$ID}{'TF_junctionAndIntraDist'} = $distTSS;
                        }  
                        
                    }
                    
                }
            }

            
        }
        $numberOfAllSites++;
    }

    close FILE;
    print "\t$TF1Filename is read!\n" if ($verbose) ;
    print "$numberOfAllSites sites\n" ;
}

#-----------read file with sites RNApolII-----
$numberOfAllSites = 0;

if ($RNApolFilename eq "") {
    print "No file with peaks of RNA pol II!\n" if ($verbose) ;
} else {
    open (FILE, "<$RNApolFilename ") or die "Cannot open file $RNApolFilename !!!!: $!";
    $_ = <FILE>;
    my @a = split /\t/;	
    my $correction = 0;
    
    if ($a[0]=~m/chr/) {
            $correction = -1;
    }
    #seek (FILE, 0, 0);
    while (<FILE>) {
        s/\R//g;			
        
        my @a = split /\t/;	
        
        my $chr = $a[1+$correction];
        my $firstPos = $a[2+$correction];
        my $LastPos = $a[3+$correction];
        my $maxPos =  $a[4+$correction];
        my $score = $a[5+$correction];
        #print "$numberOfAllSites: $score\n";
        next if ($score < $cutoff_rp); 
        
        for my $ID (keys %{$genes{$chr}}) {

            my $distTSS = ($maxPos - $genes{$chr}->{$ID}{'TSS'})*$genes{$chr}->{$ID}{'strand'};
            my $distTE = ($maxPos - $genes{$chr}->{$ID}{'TE'})*$genes{$chr}->{$ID}{'strand'};
            
            if (($distTSS>= $enhRight)&&($distTSS<=$immediateDownstream)) {
                if ($genes{$chr}->{$ID}{'RNApolScore'}<$score) {
                    $genes{$chr}->{$ID}{'RNApolScore'}=$score;
                    $genes{$chr}->{$ID}{'RNApolDist'} = $distTSS;
                }
            }
            if ($distTSS>=0 && $distTE<=0) {
                my $typeIntra = &getTypeIntra($genes{$chr}->{$ID},$maxPos);
                if ($typeIntra eq "jonction") {
                    if ($genes{$chr}->{$ID}{'RNApol_junctionScore'}<$score) {
                        $genes{$chr}->{$ID}{'RNApol_junctionScore'}=$score;
                        $genes{$chr}->{$ID}{'RNApol_junctionDist'} = $distTSS;                         
                    }                    
                }
            }   
        }

        $numberOfAllSites++;
    }
    close FILE;
    print "\t$RNApolFilename is read!\n$numberOfAllSites sites\n" if ($verbose) ;

}

#-----------read file with sites K36me3-----
$numberOfAllSites = 0;
my @K36Score;

if ($H3K36Me3polFilename eq "") {
	print "No file with peaks of H3K36me3!\n" if ($verbose) ;
} else {
    open (FILE, "<$H3K36Me3polFilename ") or die "Cannot open file $H3K36Me3polFilename !!!!: $!";
    
    $_ = <FILE>;
    my @a = split /\t/;	
    my $correction = 0;
    
    if ($a[0]=~m/chr/) {
            $correction = -1;
    }
    
    while (<FILE>) {
        s/\R//g;			
        
        my @a = split /\t/;	
        
        my $chr = $a[1+$correction];
        my $firstPos = $a[2+$correction];
        my $lastPos = $a[3+$correction];
	my $maxPos = int(($firstPos+$lastPos)/2);
	if (scalar(@a)>(4+$correction)) {
		$maxPos = $a[4+$correction];
	}	
	if ($maxPos =~ m/\D/) {
		$maxPos = int(($firstPos+$lastPos)/2);
	}
	my $score = 0;
	if (scalar(@a)>(5+$correction)) {
		$score = $a[5+$correction];
	}
        for my $ID (keys %{$genes{$chr}}) {

            my $distTSS = ($maxPos - $genes{$chr}->{$ID}{'TSS'})*$genes{$chr}->{$ID}{'strand'};
            my $distTE = ($maxPos - $genes{$chr}->{$ID}{'TE'})*$genes{$chr}->{$ID}{'strand'};

            my $scoreToadd = 0;

            if (($firstPos>=$genes{$chr}->{$ID}{'left'})&&($lastPos<=$genes{$chr}->{$ID}{'right'})) {
                $scoreToadd = $score/2.*($lastPos-$firstPos+1)/($genes{$chr}->{$ID}{'right'}-$genes{$chr}->{$ID}{'left'}+1);
            }
            if (($firstPos>=$genes{$chr}->{$ID}{'left'})&&($firstPos<$genes{$chr}->{$ID}{'right'})&&($lastPos>$genes{$chr}->{$ID}{'right'})) {
                    $scoreToadd = $score/2.*($genes{$chr}->{$ID}{'right'}-$firstPos+1)/($genes{$chr}->{$ID}{'right'}-$genes{$chr}->{$ID}{'left'}+1);
            }
            if (($firstPos<$genes{$chr}->{$ID}{'left'})&&($lastPos>$genes{$chr}->{$ID}{'left'})&&($lastPos<=$genes{$chr}->{$ID}{'right'})) {
                    $scoreToadd = $score/2.*($lastPos-$genes{$chr}->{$ID}{'left'}+1)/($genes{$chr}->{$ID}{'right'}-$genes{$chr}->{$ID}{'left'}+1);
            }
            if (($firstPos<$genes{$chr}->{$ID}{'left'})&&($lastPos>$genes{$chr}->{$ID}{'right'})) {
                    my $scoreToadd = $score/2.*($lastPos-$firstPos+1);
            }
            $genes{$chr}->{$ID}{'K36score'} += $scoreToadd;          
        }
        $numberOfAllSites++;
    }
    close FILE;
    print "\t$H3K36Me3polFilename is read!\n$numberOfAllSites sites\n" if ($verbose) ;
}

#-----------read file with sites H3K9me3-----
$numberOfAllSites = 0;


if ($H3K9Me3polFilename eq "") {
	print "No file with peaks of H3K9me3!\n" if ($verbose) ;
} else {
    open (FILE, "<$H3K9Me3polFilename ") or die "Cannot open file $H3K9Me3polFilename !!!!: $!";
    
    $_ = <FILE>;
    my @a = split /\t/;	
    my $correction = 0;
    
    if ($a[0]=~m/chr/) {
            $correction = -1;
    }
    
    while (<FILE>) {
        s/\R//g;	
        my @a = split /\t/;			
        my $chr = $a[1+$correction];
        my $firstPos = $a[2+$correction];
        my $LastPos = $a[3+$correction];
        my $maxPos =  $a[4+$correction];
        my $score = $a[5+$correction];
        
        next if ($score < $cutoff_k9); 

        for my $ID (keys %{$genes{$chr}}) {

            my $distTSS = ($maxPos - $genes{$chr}->{$ID}{'TSS'})*$genes{$chr}->{$ID}{'strand'};
            my $distTE = ($maxPos - $genes{$chr}->{$ID}{'TE'})*$genes{$chr}->{$ID}{'strand'};
            
            if (($distTSS>= $enhRight)&&($distTSS<=$immediateDownstream)) {
                if ($genes{$chr}->{$ID}{'K9promScore'}<$score) {
                    $genes{$chr}->{$ID}{'K9promScore'}=$score;
                    $genes{$chr}->{$ID}{'K9promDist'} = $distTSS;
                }
            } elsif (($distTSS >= -$K9dist)&&($distTSS<=$K9dist)) {
                if ($genes{$chr}->{$ID}{'K9largeScore'}<$score) {
                    $genes{$chr}->{$ID}{'K9largeScore'}=$score;
                    $genes{$chr}->{$ID}{'K9largeDist'} = $distTSS;
                }
            }
            if (($distTSS>= $enhLeft)&&($distTSS<$enhRight)) {
                if ($genes{$chr}->{$ID}{'K9enhScore'}<$score) {
                    $genes{$chr}->{$ID}{'K9enhScore'}=$score;
                    $genes{$chr}->{$ID}{'K9enhDist'} = $distTSS;
                }
            }
        }
        $numberOfAllSites++;
    }
    close FILE;
    print "\t$H3K9Me3polFilename is read!\n$numberOfAllSites sites\n" if ($verbose) ;
}

#-----------output all-----
#unless($initialTable eq "") {}


open (OUT , ">$outname") or die "Cannot open file $outname!!!!: $!";

print OUT "name\tchr\tstart\tend\tstrand\tReg\tfoldChange\t";

if ($GCislands ne "") {
   print OUT "GC-island\t";
}

if ( $fluoFile ne "") {
   print OUT "fluorescence\t";
}

if ($RNApolFilename ne "") {
   print OUT "RNApolII_score\tRNApolII_distTSS\tRNApol_junctionScore\tRNApol_junctionDist\t";
}
if ($H3K36Me3polFilename ne "") {
   print OUT "H3K36me3_score\t";
}
if ($H3K9Me3polFilename ne "") {
   print OUT "H3K9me3_score_prom\tH3K9me3_distTSS_prom\tH3K9me3_score_large\tH3K9me3_distTSS_large\tH3K9me3_score_enh\tH3K9me3_distTSS_enh\t";
}
if ($TF1Filename ne "") {
	print OUT "TF_score_Gene\tTF_distTSS_Gene\t";
	print OUT "TF_score_Promoter\tTF_distTSS_Promoter\t";
	print OUT "TF_score_ImmDown\tTF_distTSS_ImmDown\t";
	print OUT "TF_score_PromoterORImmDown\tTF_distTSS_PromoterORImmDown\t";
	print OUT "TF_score_Enhancer\tTF_distTSS_Enhancer\t";
	print OUT "TF_score_Intragenic\tTF_distTSS_Intragenic\t";
	print OUT "TF_score_GeneDownstream\tTF_distTSS_GeneDownstream\t";
	print OUT "TF_score_FirstExon\tTF_distTSS_FirstExon\t";
	print OUT "TF_score_FisrtIntron\tTF_distTSS_FisrtIntron\t";
	print OUT "TF_score_FirstExonAND>$immediateDownstream\tTF_distTSS_FirstExonAND>$immediateDownstream\t";
	print OUT "TF_score_FisrtIntronAND>$immediateDownstream\tTF_distTSS_FisrtIntronAND>$immediateDownstream\t";
	#print OUT "TF_score_IntraMinusFisrtIntron\tTF_distTSS_IntraMinusFisrtIntron\t";

	#print OUT "TF_score_enh60kb\tTF_distTSS_enh60kb\t";

	print OUT "TF_score_Exons2,3,4,etc\tTF_distTSS_Exons2,3,4,etc\t";
	print OUT "TF_score_Exons2,3,4,etcAND>$immediateDownstream\tTF_distTSS_Exons2,3,4,etcAND>$immediateDownstream\t";

	print OUT "TF_score_Introns2,3,4,etc\tTF_distTSS_Introns2,3,4,etc\t";
	print OUT "TF_score_Introns2,3,4,etcAND>$immediateDownstream\tTF_distTSS_Introns2,3,4,etcAND>$immediateDownstream\t";

	print OUT "TF_score_EIjunction\tTF_distTSS_EIjunction\t";
	print OUT "TF_score_EIjunctionAND>$immediateDownstream\tTF_distTSS_EIjunctionAND>$immediateDownstream";
	
	if ($ifHMmode) {print OUT "\tK27TSS_Score\tK27TSS_Dist\tK27GeneB_Score\tK27GeneB_Dist"}
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

	if ($RNApolFilename ne "") {
	     print OUT "$genes{$chr}->{$ID}{'RNApolScore'}\t$genes{$chr}->{$ID}{'RNApolDist'}\t$genes{$chr}->{$ID}{'RNApol_junctionScore'}\t$genes{$chr}->{$ID}{'RNApol_junctionDist'}\t";
	}
	if ($H3K36Me3polFilename ne "") {
	   print OUT "$genes{$chr}->{$ID}{'K36score'}\t";
	}
	if ($H3K9Me3polFilename ne "") {
           print OUT "$genes{$chr}->{$ID}{'K9promScore'}\t$genes{$chr}->{$ID}{'K9promDist'}\t$genes{$chr}->{$ID}{'K9largeScore'}\t$genes{$chr}->{$ID}{'K9largeDist'}\t$genes{$chr}->{$ID}{'K9enhScore'}\t$genes{$chr}->{$ID}{'K9enhDist'}\t";
	}
	if ($TF1Filename ne "") {  
		print OUT "$genes{$chr}->{$ID}{'TFallScore'}\t$genes{$chr}->{$ID}{'TFallDist'}\t";
		print OUT "$genes{$chr}->{$ID}{'TFpromSimpleScore'}\t$genes{$chr}->{$ID}{'TFpromSimpleDist'}\t";
		print OUT "$genes{$chr}->{$ID}{'TFImmDownScore'}\t$genes{$chr}->{$ID}{'TFImmDownDist'}\t";
		print OUT "$genes{$chr}->{$ID}{'TFpromScore'}\t$genes{$chr}->{$ID}{'TFpromDist'}\t";
		print OUT "$genes{$chr}->{$ID}{'TFenhScore'}\t$genes{$chr}->{$ID}{'TFenhDist'}\t";
		print OUT "$genes{$chr}->{$ID}{'TFintraScore'}\t$genes{$chr}->{$ID}{'TFintraDist'}\t";
		print OUT "$genes{$chr}->{$ID}{'TF5kbDownScore'}\t$genes{$chr}->{$ID}{'TF5kbDownDist'}\t";
		print OUT "$genes{$chr}->{$ID}{'TF_FirstExonScore'}\t$genes{$chr}->{$ID}{'TF_FirstExonDist'}\t";
		print OUT "$genes{$chr}->{$ID}{'TF_FirstExonAndIntraScore'}\t$genes{$chr}->{$ID}{'TF_FirstExonAndIntraDist'}\t";
		print OUT "$genes{$chr}->{$ID}{'TFFirstIntronScore'}\t$genes{$chr}->{$ID}{'TFFirstIntronDist'}\t";
		print OUT "$genes{$chr}->{$ID}{'TFFirstIntronAndIntraScore'}\t$genes{$chr}->{$ID}{'TFFirstIntronAndIntraDist'}\t";
		#print OUT "$genes{$chr}->{$ID}{'TFintraMinusFirstIntronScore'}\t$genes{$chr}->{$ID}{'TFintraMinusFirstIntronDist'}\t";

		#print OUT "$genes{$chr}->{$ID}{'TFenh60kbScore'}\t$genes{$chr}->{$ID}{'TFenh60kbDist'}\t";

		print OUT "$genes{$chr}->{$ID}{'TF_OtherExonsScore'}\t$genes{$chr}->{$ID}{'TF_OtherExonsDist'}\t";
		print OUT "$genes{$chr}->{$ID}{'TF_OtherExonsAndIntraScore'}\t$genes{$chr}->{$ID}{'TF_OtherExonsAndIntraDist'}\t";

		print OUT "$genes{$chr}->{$ID}{'TF_OtherIntronsScore'}\t$genes{$chr}->{$ID}{'TF_OtherIntronsDist'}\t";
		print OUT "$genes{$chr}->{$ID}{'TF_OtherIntronsAndIntraScore'}\t$genes{$chr}->{$ID}{'TF_OtherIntronsAndIntraDist'}\t";

		print OUT "$genes{$chr}->{$ID}{'TF_junctionScore'}\t$genes{$chr}->{$ID}{'TF_junctionDist'}\t";
		print OUT "$genes{$chr}->{$ID}{'TF_junctionAndIntraScore'}\t$genes{$chr}->{$ID}{'TF_junctionAndIntraDist'}";
		
		if ($ifHMmode) {print OUT "\t$genes{$chr}->{$ID}{'K27TSS_Score'}\t$genes{$chr}->{$ID}{'K27TSS_Dist'}\t$genes{$chr}->{$ID}{'K27GeneB_Score'}\t$genes{$chr}->{$ID}{'K27GeneB_Dist'}"}

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
