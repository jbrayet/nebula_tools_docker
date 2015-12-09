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

#!/usr/bin/perl -w

#for a complete list of genes and a list of sites, sorts genes close to sites
#printDistance

#read only gene files with refSeqGenes
#counts only once overlapping transcripts or genes

#calculates some stats about locations of peaks

#uses hierarchie : promoter, imUpstream, intragenic,enh, 5kbdownstream, intergenic + for genes: fExon,exon,fIntron,intron,junction+-50kb ONLY FOR "findClosestGene"
#outputs fold change of expression is available

#outputname - corrected
#all isoforms are considered, even those that start and end at the same coordinates
#uses fluorescence values

#only one location per gene in mode "all" (but promoter+intragenic if both)

#read directly Noresp/down/up genes

#prints distance to TE

use strict;
use POSIX;

my $SitesFilename ;
my $GenesFilename ;
my $SelectedGenesFilename = "";
my $MirFilename = "";
my $minScore = 0;
my $BIGNUMBER = 10000000000;	
my $enhLeft = -30000;
my $enhRight = -1500;
my $promLeft = $enhRight ;
my $immediateDownstream = 2000;
my $maxScore = $BIGNUMBER ;
my $verbose = 0;
my $maxDistance = $BIGNUMBER ;
my $downStreamDist = 5000;
my $jonctionSize = 50;
my $all = 0;
my $fluoFile = "";
my $regType = "";
my $GCislands = "";
my $ResFilename = "";

my $usage = qq{
    $0

    -----------------------------
    mandatory parameters:
    
    -g   filename 	file with all genes   
    -tf filename 	file with sites of TF 
    
    -----------------------------
    optional parameters:
    -v				for verbose
    -mir 	filename	file with positions of miRNA
    -maxDist	value		maximal Distance to genes (def. Inf)
    -minScore	value		minimal Score (def. 0)
    -maxScore	value		maximal Score (def. Inf)  
    -selG	filename	file with selected genes and their expression levels
    -all 			to output all genes intersecting with the peak
    -fluo	filename	file with fluorescence
    -regType	valeu		down-regulated, up-regulated, no-response
    -gc         filename        file with gc-islands

    -lp         valeu        	upstream of TSS region to define promoter
    -rp         valeu        	downstream of TSS region to define immediate downstream
    -enh         valeu        	upstream of TSS region to define enhancer
    -dg         valeu        	downstream of transcription end region to define gene downstream

    -o		filename	output filename

};

if(scalar(@ARGV) == 0){
    print $usage;
    exit(0);
}


while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;
    if ( $this_arg eq '-h') {print "$usage\n"; exit; }

    elsif ( $this_arg eq '-g') {$GenesFilename = shift @ARGV;}
    elsif ( $this_arg eq '-tf') {$SitesFilename = shift @ARGV;}
    elsif ( $this_arg eq '-mir') {$MirFilename = shift @ARGV;}   
    elsif ( $this_arg eq '-selG') {$SelectedGenesFilename = shift @ARGV;}  
    elsif ( $this_arg eq '-maxDist') {$maxDistance = 1;}
    elsif ( $this_arg eq '-maxScore') {$maxScore = shift @ARGV;}   
    elsif ( $this_arg eq '-minScore') {$minScore = shift @ARGV;}  
    elsif ( $this_arg eq '-v') {$verbose = 1;}  
    elsif ( $this_arg eq '-all') {$all = 1;}
    elsif ( $this_arg eq '-regType') {$regType = shift @ARGV;}

    elsif ( $this_arg eq '-fluo') {$fluoFile = shift @ARGV;}
    elsif ( $this_arg eq '-gc') {$GCislands = shift @ARGV;}  
    elsif ( $this_arg eq '-o') {$ResFilename = shift @ARGV;}     

    elsif ( $this_arg eq '-lp') {$enhRight = shift @ARGV;$promLeft = $enhRight;}
    elsif ( $this_arg eq '-rp') {$immediateDownstream = shift @ARGV;}  
    elsif ( $this_arg eq '-enh') {$enhLeft = shift @ARGV;}    
    elsif ( $this_arg eq '-dg') {$downStreamDist = shift @ARGV;}      

    elsif ( $this_arg =~ m/^-/ ) { print "unknown flag: $this_arg\n";}
}

if ($maxDistance<0) {
	die "\tThe maximal distance must be positive!\n";
}
my $tmpn = $SitesFilename;
$tmpn =~ s/.*[\/\\]//;
if ($ResFilename eq "") {
	$ResFilename = $SelectedGenesFilename.$tmpn."_".$minScore."_dists.txt";

	if ($maxScore != $BIGNUMBER) {	
		$ResFilename = $SelectedGenesFilename.$SitesFilename."_".$minScore."_$maxScore"."_dists.txt";
	}

	if ($all == 1) {
	    $ResFilename =~ s/_dists/_all_dists/;
	}
}


#-----------read selected genes----------------
my %selectedGenes;
my %selectedGenesFoldChange;
if ( $SelectedGenesFilename ne "") {
    open (FILE, "<$SelectedGenesFilename ") or die "Cannot open file $SelectedGenesFilename !!!!: $!";
    while (<FILE>) {
        chomp;
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
        chomp;
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
    print "\t\t$fluoFile is read!\n" if ($verbose);
}
#-----------read GC-islands----------------
my %GCislands;
if ($GCislands ne "") {
    open (FILE, "<$GCislands ") or die "Cannot open file $GCislands !!!!: $!";
    
     while (<FILE>) {
        chomp;
        my @a = split/\t/;
        #bin	chrom	chromStart	chromEnd	name	length	cpgNum	gcNum	perCpg	perGc	obsExp
        #107	chr1	36568608	36569851	CpG: 128	1243	128	766	20.6	61.6	1.09
        my $chr = $a[1];
        my $start = $a[2];
        my $end = $a[3];
        $GCislands{$chr}->{$start}=$end;
     }
    close FILE;
} elsif ($verbose) {
    print "you did not specify a file with GC-islands\n";    
}

#-----------read genes----------------

my %genes;

my $count = 0;

open (GENES, "<$GenesFilename") or die "Cannot open file $GenesFilename!!!!: $!";
<GENES>;
while (<GENES>) {
	chomp;	
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
		my $ID = "$name\t$chr:$leftPos-$rightPos;$exonStarts-$exonEnds";
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
		unless (exists($genes{$chr})) {
			my %h;
			$genes{$chr} = \%h;
		}
		unless (exists($genes{$chr}->{$ID})) {
			my %h1;
			$genes{$chr}->{$ID} = \%h1;	
			$count++;
		}
		if ( $fluoFile ne "") {                
		    if (exists($fluoGenes{$name})) {
			$fluo = $fluoGenes{$name};
		    }
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
		($genes{$chr}->{$ID}{'firstIntronStart'},$genes{$chr}->{$ID}{'firstIntronEnd'}) = getFirstIntron ($exonCount,$exonStarts,$exonEnds,$strand);
		($genes{$chr}->{$ID}{'firstExonStart'},$genes{$chr}->{$ID}{'firstExonEnd'}) = getFirstExon ($exonCount,$exonStarts,$exonEnds,$strand);
	
		$genes{$chr}->{$ID}{'fluo'} = $fluo;
		
		$genes{$chr}->{$ID}{'GCisland'} = 0;
		if ($GCislands ne "") {
		    $genes{$chr}->{$ID}{'GCisland'} = checkIfGC ($genes{$chr}->{$ID}{'TSS'},$strand,2000,$GCislands{$chr});
		}
			
	}
}
print "Total genes : $count\n";
close GENES;
print "\t\t$GenesFilename is read!\n" if ($verbose);
#for my $gName (sort keys %{$genes{'chr18'}}) {

 #    print "$gName\t$genes{'chr18'}->{$gName}{'TSS'}\n";
#}

#-----------read file with sites miRNA, store as genes-----

if ( $MirFilename eq ""){
    print "you did not specify file with miRNA\n" if ($verbose);;
}
else {


	open (MIR, "<$MirFilename ") or die "Cannot open file $MirFilename !!!!: $!";
	#chr1	20669090	20669163	mmu-mir-206	960	+
	while (<MIR>) {
		chomp;	
		my ($name, $chr, $leftPos, $rightPos, $strand );
		my $fluo = "NA";

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
		my $ID = "$name\t$chr:$leftPos-$rightPos";

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
		($genes{$chr}->{$ID}{'firstIntronStart'},$genes{$chr}->{$ID}{'firstIntronEnd'}) = (0,0);
		($genes{$chr}->{$ID}{'firstExonStart'},$genes{$chr}->{$ID}{'firstExonEnd'}) = (0,0);

		$genes{$chr}->{$ID}{'fluo'} = $fluo;
		$genes{$chr}->{$ID}{'GCisland'} = 0;


		
	}


	close MIR;
	print "\t\t$MirFilename is read!\n" if ($verbose) ;
}
#-----------read file with sites and find N closest genes-----
my $numberOfAllSites = 0;

open (FILE, "<$SitesFilename") or die "Cannot open file $SitesFilename!!!!: $!";
open (OUT , ">$ResFilename") or die "Cannot open file $ResFilename!!!!: $!";
print OUT "Chromosome\tStart\tEnd\tMax\tScore\tDistTSS\tType\tTypeIntra\tReg\tFoldChange\tDistTE\tGeneName\tGeneCoordinates\n" ;
my $header = <FILE>; #read header
my @a = split (/\t/,$header );
my $correction = 0;
if ($a[1] =~m/chr/) {
	$correction=1;
}

unless ($header=~m/Chromosome/ || $header=~m/track/|| $header=~m/Start/i) {
	close FILE;
	open (FILE, "<$SitesFilename") or die "Cannot open file $SitesFilename!!!!: $!";
}

$count = 0;

my %typehash;
my %typehashIntra;

$typehashIntra{"f_exon"} = 0;
$typehashIntra{"f_intron"} = 0;
$typehashIntra{"jonction"} = 0;
$typehashIntra{"exon"} = 0;
$typehashIntra{"intron"} = 0;


while (<FILE>) {
	chomp;		
	$numberOfAllSites++;	
	my @a = split /\t/;	
	my $chro = $a[0+$correction]; 
	my $fpos = $a[1+$correction];
	my $lpos = $a[2+$correction];
	my $maxpos = $a[3+$correction];
	if ($maxpos =~ m/\D/) {
		$maxpos = int(($lpos+$fpos)/2);
	}
	my $score = $a[4+$correction];
	next if ($score < $minScore);
	next if ($score >= $maxScore);
	$score = $score +1 -1;
    	#print  "$score" ;	

	#my $site = "$chro:$fpos-$lpos\t$maxpos\t$score";
	#print $site ;
	#my $motifNumber = $a[3];
	
	#my $geneSet = &search_genes($N,$chro,($fpos+$lpos)/2,\%genes);
	my @b;
	unless ($all) {
		@b = &findClosestGene($fpos,$lpos,$score,$chro,$maxpos,\%genes); #findAllClosestGene for all genes overlaping a peak
	} else {
		@b = &findAllClosestGeneOneLocation ($fpos,$lpos,$score,$chro,$maxpos,\%genes); #findAllClosestGene
	}
		#print "$distTSS\t$distTE\n" unless ($distTSS == $BIGNUMBER);	
	for my $type (@b) {
		$typehash{$type} = 0 unless (exists($typehash{$type}));
		$typehash{$type}++;
	}
	$count ++;	

}
close FILE;
close OUT ;

print "Total Sites: $count\n";
my $nEntries = 0;
for my $type (sort keys %typehash) {
	$nEntries += $typehash{$type};
}
print "Type\tCount\tFrequency\n";
for my $type (sort keys %typehash) {
	print $type,"\t",$typehash{$type},"\t",$typehash{$type}/$nEntries,"\n";
}
for my $type (sort keys %typehashIntra) {
	print $type,"\t",$typehashIntra{$type},"\t",$typehashIntra{$type}/$nEntries,"\n";
}
 
 
#my @arr = @{$genes{'chr'}};

sub findClosestGene {	
	my ($fpos,$lpos,$score,$chro,$pos,$genes) = @_;
	my ($distTSS,$distTE,$type,$reg,$geneName,$foldChange) = ($BIGNUMBER,$BIGNUMBER,"intergenic",0,0,0);
	my $locDistTSS = 0;
	my $locDistTE  = 0;
	my %hashForOverlapingGenes;
	my @array2return;
	my $typeIntra = "NA";
	for my $gName (keys %{$genes->{$chro}}) {
		my $TSS = $genes->{$chro}{$gName}{'TSS'};	
		my $strand = $genes->{$chro}{$gName}{'strand'};
	    	my $TE = $genes->{$chro}{$gName}{'TE'};	
	    	
	    	$locDistTE = ($pos-$TE)*$strand;
		$locDistTSS = ($pos-$TSS)*$strand;
		#print "$pos-$genes->{$chro}{$gName}{'TSS'})*$genes->{$chro}{$gName}{'strand'} = $locDistTSS\n";
		#my $locDistTE = ($pos-$genes->{$chro}{$gName}{'TE'})*$genes->{$chro}{$gName}{'strand'};
		
		#print "$gName $maxDistance $pos-$genes->{$chro}{$gName}{'TSS'} $genes->{$chro}{$gName}{'strand'}\n";
		#print "$locDistTSS $enhLeft $enhRight\n";
				
		if (($locDistTSS>= $enhLeft)&&($locDistTSS<$enhRight)) {
		    $type = "enhancer";			
		} elsif (($locDistTSS>= $enhRight)&&($locDistTSS<=0)) {
		    $type = "promoter";
		} elsif (($locDistTSS> 0)&&($locDistTSS<=$immediateDownstream)&&($locDistTE<=0)) {
		    $type = "immediateDownstream";			
		} elsif (($locDistTSS >= $immediateDownstream)&&($locDistTE<=0)) {
		    $type = "intragenic"; 
		} elsif (($locDistTE >= 0)&&($locDistTE<=$downStreamDist)) {
			$type = "5kbDownstream";		
		} elsif (abs($locDistTSS)<abs($distTSS)) {
			$distTSS = $locDistTSS;
			$distTE = $locDistTE;
			($reg,$geneName,$foldChange) = ($genes->{$chro}{$gName}{'reg'},$gName,$genes->{$chro}{$gName}{'foldChange'});
		}
		
		if (($type ne "intergenic")&&($type ne "NA")) {
			unless (exists ($hashForOverlapingGenes{$type})) {
				my %h;
				$hashForOverlapingGenes{$type} = \%h;
			}
			$hashForOverlapingGenes{$type}->{$locDistTSS}=$gName;

			#print OUT "$chro\t$fpos\t$lpos\t$pos\t$score\t$locDistTSS\t$type\t$genes->{$chro}{$gName}{'reg'}\t$gName\n" ;
			$type = "NA";
			#unless ($gName =~ m/$TSS/) {
			 #   print "$gName\t$TSS\n";
			#}
		}
		
	}
	if ($type eq "intergenic") {
		print OUT "$chro\t$fpos\t$lpos\t$pos\t$score\t$distTSS\t$type\t$typeIntra\t$reg\t$foldChange\t$distTE\t$geneName\n" ;
		#print  "$chro\t$fpos\t$lpos\t$pos\t$score\t$distTSS\t$type\t$reg\t$geneName\n" ;
		push(@array2return,$type);

	} else {
		my $selectedType;
		if (exists ($hashForOverlapingGenes{"promoter"})) {
			$selectedType = "promoter";
		} elsif (exists ($hashForOverlapingGenes{"immediateDownstream"})) {
			$selectedType = "immediateDownstream";
		}  elsif (exists ($hashForOverlapingGenes{"intragenic"})) {
			$selectedType = "intragenic";	
		}   elsif (exists ($hashForOverlapingGenes{"enhancer"})) {
			$selectedType = "enhancer";			
		} elsif (exists ($hashForOverlapingGenes{"5kbDownstream"})) {
			$selectedType = "5kbDownstream";		
		}
		my $bestGene = printBest($hashForOverlapingGenes{$selectedType},$chro,$fpos,$lpos,$pos,$score,$genes,$selectedType); #print "$type\n";		
		push(@array2return,$selectedType);		
		my %otherGenes;
		for my $type ("promoter","immediateDownstream","intragenic","enhancer","5kbDownstream") {  #"firstIntron","exon","intron",
			for my $locDistTSS (sort {$a<=>$b} keys %{$hashForOverlapingGenes{$type}}) {
				my $otherGene = $hashForOverlapingGenes{$type}->{$locDistTSS};
				next if ($genes->{$chro}{$bestGene}{'name'} eq $genes->{$chro}{$otherGene}{'name'});
				next if (exists ($otherGenes{$genes->{$chro}{$otherGene}{'name'}}));
				if (overlapGenes($otherGene,$bestGene,$chro,$genes)) {
				    
					if (($type eq "intragenic")||($type eq "immediateDownstream")) {
					    $typeIntra = &getTypeIntra($genes->{$chro}{$otherGene},$pos);
					    $typehashIntra{$typeIntra}++;
					}else {
					    $typeIntra = "NA";					}				    
					print OUT "$chro\t$fpos\t$lpos\t$pos\t$score\t$locDistTSS\t$type\t$typeIntra\t$genes->{$chro}{$otherGene}{'reg'}\t$genes->{$chro}{$otherGene}{'foldChange'}\t",($pos-$genes->{$chro}{$otherGene}{'TE'})*$genes->{$chro}{$otherGene}{'strand'},"\t$otherGene\n" ;
					push(@array2return,$type);
					
					#print "$chro\t$fpos\t$lpos\t$pos\t$score\t$locDistTSS\t$type\t$genes->{$chro}{$otherGene}{'reg'}\t$otherGene\n" ;
					#print $type,"\n";
					$otherGenes{$genes->{$chro}{$otherGene}{'name'}} = 1;
				}
			}	
		}
		
	}
	return @array2return;
}

sub findAllClosestGeneOneLocation {	
	my ($fpos,$lpos,$score,$chro,$pos,$genes) = @_;
	my ($distTSS,$type,$reg,$geneName,$foldChange) = ($BIGNUMBER,"intergenic",0,0,0);
	my $locDistTSS = 0;
	my $locDistTE  = 0;
	my %hashForOverlapingGenes;
	my @array2return;
	my $typeIntra = "NA";
	for my $gName (keys %{$genes->{$chro}}) {
		my $TSS = $genes->{$chro}{$gName}{'TSS'};	
		my $strand = $genes->{$chro}{$gName}{'strand'};
	    	my $TE = $genes->{$chro}{$gName}{'TE'};	
	    	
	    	$locDistTE = ($pos-$TE)*$strand;
		$locDistTSS = ($pos-$TSS)*$strand;
		#print "$pos-$genes->{$chro}{$gName}{'TSS'})*$genes->{$chro}{$gName}{'strand'} = $locDistTSS\n";
		#my $locDistTE = ($pos-$genes->{$chro}{$gName}{'TE'})*$genes->{$chro}{$gName}{'strand'};
		
		#print "$gName $maxDistance $pos-$genes->{$chro}{$gName}{'TSS'} $genes->{$chro}{$gName}{'strand'}\n";
		#print "$locDistTSS $enhLeft $enhRight\n";
				
		if (($locDistTSS>= $enhLeft)&&($locDistTSS<$enhRight)) {
		    $type = "enhancer";			
		} elsif (($locDistTSS>= $enhRight)&&($locDistTSS<=0)) {
		    $type = "promoter";
		} elsif (($locDistTSS> 0)&&($locDistTSS<=$immediateDownstream)&&($locDistTE<=0)) {
		    $type = "immediateDownstream";			
		} elsif (($locDistTSS >= $immediateDownstream)&&($locDistTE<=0)) {
		    $type = "intragenic"; 
		} elsif (($locDistTE >= 0)&&($locDistTE<=$downStreamDist)) {
			$type = "5kbDownstream";		
		} elsif (abs($locDistTSS)<abs($distTSS)) {
			$distTSS = $locDistTSS;
			($reg,$geneName,$foldChange) = ($genes->{$chro}{$gName}{'reg'},$gName,$genes->{$chro}{$gName}{'foldChange'});
		}
		
		if (($type ne "intergenic")&&($type ne "NA")) {
			unless (exists ($hashForOverlapingGenes{$type})) {
			    my %h;
			    $hashForOverlapingGenes{$type} = \%h;			
			    $hashForOverlapingGenes{$type}->{$gName}=$locDistTSS;
			    
			    if (($type eq "intragenic")||($type eq "immediateDownstream")) {
				$typeIntra = &getTypeIntra($genes->{$chro}{$gName},$pos);
				$typehashIntra{$typeIntra}++;
			    }	
			    print OUT "$chro\t$fpos\t$lpos\t$pos\t$score\t$locDistTSS\t$type\t$typeIntra\t$genes->{$chro}{$gName}{'reg'}\t$genes->{$chro}{$gName}{'fluo'}\t$genes->{$chro}{$gName}{'foldChange'}\t$genes{$chro}->{$gName}{'GCisland'}\t",($pos-$genes->{$chro}{$gName}{'TE'})*$genes->{$chro}{$gName}{'strand'},"\t$gName\n" ; 
			    push(@array2return,$type);	
			}
			$typeIntra = "NA";
			$type = "NA";
			#unless ($gName =~ m/$TSS/) {
			 #   print "$gName\t$TSS\n";
			#}
		}
		
	}
	if ($type eq "intergenic") {
		print OUT "$chro\t$fpos\t$lpos\t$pos\t$score\t$locDistTSS\t$type\t$typeIntra\t$genes->{$chro}{$geneName}{'reg'}\t$genes->{$chro}{$geneName}{'fluo'}\t$genes->{$chro}{$geneName}{'foldChange'}\t$genes{$chro}->{$geneName}{'GCisland'}\t",($pos-$genes->{$chro}{$geneName}{'TE'})*$genes->{$chro}{$geneName}{'strand'},"\t$geneName\n" ; 

		#print  "$chro\t$fpos\t$lpos\t$pos\t$score\t$distTSS\t$type\t$reg\t$geneName\n" ;
		push(@array2return,$type);

	} 					
	return @array2return;
}

sub findAllClosestGene {	
	my ($fpos,$lpos,$score,$chro,$pos,$genes) = @_;
	my ($distTSS,$type,$reg,$geneName,$foldChange) = ($BIGNUMBER,"intergenic",0,0,0);
	my $locDistTSS = 0;
	my $locDistTE  = 0;
	my %hashForOverlapingGenes;
	my @array2return;
	my $typeIntra = "NA";
	for my $gName (keys %{$genes->{$chro}}) {
		my $TSS = $genes->{$chro}{$gName}{'TSS'};	
		my $strand = $genes->{$chro}{$gName}{'strand'};
	    	my $TE = $genes->{$chro}{$gName}{'TE'};	
	    	
	    	$locDistTE = ($pos-$TE)*$strand;
		$locDistTSS = ($pos-$TSS)*$strand;
		#print "$pos-$genes->{$chro}{$gName}{'TSS'})*$genes->{$chro}{$gName}{'strand'} = $locDistTSS\n";
		#my $locDistTE = ($pos-$genes->{$chro}{$gName}{'TE'})*$genes->{$chro}{$gName}{'strand'};
		
		#print "$gName $maxDistance $pos-$genes->{$chro}{$gName}{'TSS'} $genes->{$chro}{$gName}{'strand'}\n";
		#print "$locDistTSS $enhLeft $enhRight\n";
				
		if (($locDistTSS>= $enhLeft)&&($locDistTSS<$enhRight)) {
		    $type = "enhancer";			
		} elsif (($locDistTSS>= $enhRight)&&($locDistTSS<=0)) {
		    $type = "promoter";
		} elsif (($locDistTSS> 0)&&($locDistTSS<=$immediateDownstream)&&($locDistTE<=0)) {
		    $type = "immediateDownstream";			
		} elsif (($locDistTSS >= $immediateDownstream)&&($locDistTE<=0)) {
		    $type = "intragenic"; 
		} elsif (($locDistTE >= 0)&&($locDistTE<=$downStreamDist)) {
			$type = "5kbDownstream";		
		} elsif (abs($locDistTSS)<abs($distTSS)) {
			$distTSS = $locDistTSS;
			($reg,$geneName,$foldChange) = ($genes->{$chro}{$gName}{'reg'},$gName,$genes->{$chro}{$gName}{'foldChange'});
		}
		
		if (($type ne "intergenic")&&($type ne "NA")) {
			unless (exists ($hashForOverlapingGenes{$type})) {
				my %h;
				$hashForOverlapingGenes{$type} = \%h;
			}
			#$hashForOverlapingGenes{$type}->{$gName}=$locDistTSS;
			
			$typeIntra = "NA";		
			if (($type eq "intragenic")||($type eq "immediateDownstream")) {
				$typeIntra = &getTypeIntra($genes->{$chro}{$gName},$pos);
				$typehashIntra{$typeIntra}++;
			}	
			print OUT "$chro\t$fpos\t$lpos\t$pos\t$score\t$locDistTSS\t$type\t$typeIntra\t$genes->{$chro}{$gName}{'reg'}\t$genes->{$chro}{$gName}{'fluo'}\t$genes->{$chro}{$gName}{'foldChange'}\t$genes{$chro}->{$gName}{'GCisland'}\t",($pos-$genes->{$chro}{$gName}{'TE'})*$genes->{$chro}{$gName}{'strand'},"\t$gName\n" ; 
			push(@array2return,$type);
			
			$type = "NA";
			#unless ($gName =~ m/$TSS/) {
			 #   print "$gName\t$TSS\n";
			#}
		}
		
	}
	if ($type eq "intergenic") {
		print OUT "$chro\t$fpos\t$lpos\t$pos\t$score\t$locDistTSS\t$type\t$typeIntra\t$genes->{$chro}{$geneName}{'reg'}\t$genes->{$chro}{$geneName}{'fluo'}\t$genes->{$chro}{$geneName}{'foldChange'}\t$genes{$chro}->{$geneName}{'GCisland'}\t",($pos-$genes->{$chro}{$geneName}{'TE'})*$genes->{$chro}{$geneName}{'strand'},"\t$geneName\n" ; 

		#print  "$chro\t$fpos\t$lpos\t$pos\t$score\t$distTSS\t$type\t$reg\t$geneName\n" ;
		push(@array2return,$type);

	} 					
	return @array2return;
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

sub overlapGenes {
	my ($otherGene,$bestGene,$chro,$genes)= @_;
	my $a1 = $genes->{$chro}{$otherGene}{'left'};
	my $a2 = $genes->{$chro}{$bestGene}{'left'};
	my $e1 = $genes->{$chro}{$otherGene}{'right'};
	my $e2 = $genes->{$chro}{$bestGene}{'right'};
	if (($a1 >= $a2)&&($a1 <= $e2)) {
		return 1;
	}
	if (($a2 >= $a1)&&($a2 <= $e1)) {
		return 1;
	}
	return 0;
}

sub printBest {
	my ($hash,$chro,$fpos,$lpos,$pos,$score,$genes,$type) = @_;
	for my $locDistTSS (sort {$a<=>$b} keys %{$hash}) {
		my $gName = $hash->{$locDistTSS};
		my $typeIntra = "NA";		
		if (($type eq "intragenic")||($type eq "immediateDownstream")) {
			$typeIntra = &getTypeIntra($genes->{$chro}{$gName},$pos);
			$typehashIntra{$typeIntra}++;
		}	
		print OUT "$chro\t$fpos\t$lpos\t$pos\t$score\t$locDistTSS\t$type\t$typeIntra\t$genes->{$chro}{$gName}{'reg'}\t$genes->{$chro}{$gName}{'foldChange'}\t",($pos-$genes->{$chro}{$gName}{'TE'})*$genes->{$chro}{$gName}{'strand'},"\t$gName\n" ;		
		return $gName;
	}
}

#sub findAllClosestGene {	
#	my ($fpos,$lpos,$score,$chro,$pos,$genes) = @_;
#	my ($distTSS,$type,$reg,$geneName,$foldChange) = ($BIGNUMBER,"intergenic",0,0,0);
#	my $locDistTSS = 0;
#	my $locDistTE = 0;
#	my @array2return;
#	for my $gName (keys %{$genes->{$chro}}) {
#		my $TSS = $genes->{$chro}{$gName}{'TSS'};	
#		my $strand = $genes->{$chro}{$gName}{'strand'};
#	    	my $TE = $genes->{$chro}{$gName}{'TE'};	
#		$locDistTSS = ($pos-$TSS)*$strand;
#		$locDistTE = ($pos-$TE)*$strand;
#		#print "$pos-$genes->{$chro}{$gName}{'TSS'})*$genes->{$chro}{$gName}{'strand'} = $locDistTSS\n";
#		#my $locDistTE = ($pos-$genes->{$chro}{$gName}{'TE'})*$genes->{$chro}{$gName}{'strand'};
#		
#		#print "$gName $maxDistance $pos-$genes->{$chro}{$gName}{'TSS'} $genes->{$chro}{$gName}{'strand'}\n";
#		#print "$locDistTSS $enhLeft $enhRight\n";
#				
#		if (($locDistTSS>= $enhLeft)&&($locDistTSS<$enhRight)) {
#			$type = "enhancer";
#		} elsif (($locDistTSS>= $enhRight)&&($locDistTSS<=0)) {
#			$type = "promoter";
#		} elsif (($locDistTSS> 0)&&($locDistTSS<=$immediateDownstream)) {
#			$type = "immediateDownstream";			
#		} elsif (($pos >= $genes{$chro}->{$gName}{'firstIntronStart'})&&($pos <= $genes{$chro}->{$gName}{'firstIntronEnd'})) {
#			$type = "firstIntron";						
#		} elsif (($locDistTSS> $immediateDownstream)&&($locDistTSS<=$genes{$chro}->{$gName}{'length'})) {
#			#$type = "intragenic";
#			$type = &getIntronExon ($pos, $genes{$chro}->{$gName}{'exonCount'},$genes{$chro}->{$gName}{'exonStarts'},$genes{$chro}->{$gName}{'exonEnds'},$genes{$chro}->{$gName}{'strand'});
#		} elsif (($locDistTE >= 0)&&($locDistTE<=$downStreamDist)) {
#			$type = "5kbDownstream";		
#		} elsif (abs($locDistTSS)<abs($distTSS)) {
#			$distTSS = $locDistTSS;
#			($reg,$geneName,$foldChange) = ($genes->{$chro}{$gName}{'reg'},$gName,$genes->{$chro}{$gName}{'foldChange'});
#		}
#		if (($type ne "intergenic")&&($type ne "NA")) {
#			print OUT "$chro\t$fpos\t$lpos\t$pos\t$score\t$locDistTSS\t$type\t$genes->{$chro}{$gName}{'reg'}\t$genes->{$chro}{$gName}{'foldChange'}\t$gName\n" ;
#			push(@array2return,$type);
#			$type = "NA";
#			#unless ($gName =~ m/$TSS/) {
#			 #   print "$gName\t$TSS\n";
#			#}
#		}
#		
#	}
#	if ($type eq "intergenic") {
#		print OUT "$chro\t$fpos\t$lpos\t$pos\t$score\t$distTSS\t$type\t$reg\t$foldChange\t$geneName\n" ;
#		push(@array2return,$type);
#		#print  "$chro\t$fpos\t$lpos\t$pos\t$score\t$distTSS\t$type\t$reg\t$geneName\n" ;
#
#	}
#	return @array2return;
#}

sub getFirstIntron {
	my ($exonCount,$exonStarts,$exonEnds,$strand) = @_;
	my ($left,$right);
	if ($exonCount == 1) {
		return (0,0);
	}
	if ($strand == 1) {
		$left = (split ",", $exonEnds)[0]+$jonctionSize;
		$right = (split (",", $exonStarts))[1]-$jonctionSize;
	} else {
		$left = (split (",", $exonEnds))[$exonCount-2]+$jonctionSize;
		$right = (split (",", $exonStarts))[$exonCount-1]-$jonctionSize;	
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

sub min {
	return $_[0] if ($_[0]<$_[1]);
	$_[1];
}
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
    }
    for my $leftGC (keys %{$GCislandsChr}) {
	my $rightGC = $GCislandsChr->{$leftGC};
	if ($leftGC>=$leftProm&&$leftGC<=$rightProm || $rightGC>=$leftProm&&$rightGC<=$rightProm) {
	    return "GC-island";
	}
    }       
    return $ifGC ;
}
