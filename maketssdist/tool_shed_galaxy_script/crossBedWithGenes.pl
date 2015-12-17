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

use strict;

my $SitesFilename ;
my $SitesFilename2 ="";
my $GenesFilename ;
my $MirFilename = "";
my $BIGNUMBER = 10000000;	
my $verbose = 0;
my $header = 0;
my $noXY = 0;
my $minscore = 0;

my $usage = qq{
    $0

    -----------------------------
    mandatory parameters:
    
    -g   filename 	file with all genes   
    -f filename 	file with sites in BED format
    
    -----------------------------
    optional parameters:
    -v				for verbose
    -minScore			minimal score
    -mir 	filename	file with positions of miRNA
    -head			if there is a header
    -reg	filename	file with FoldChanges or expression values and annotation (should in colomn 4)
    -noXY
    -o		filename	output file
};

if(scalar(@ARGV) == 0){
    print $usage;
    exit(0);
}

my $ResFilename = "";

my $regFilename = "";

while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;
    if ( $this_arg eq '-h') {print "$usage\n"; exit; }

    elsif ( $this_arg eq '-g') {$GenesFilename = shift @ARGV;}
    elsif ( $this_arg eq '-f') {$SitesFilename = shift @ARGV;}
    elsif ( $this_arg eq '-mir') {$MirFilename = shift @ARGV;}   
    elsif ( $this_arg eq '-v') {$verbose = 1;}  
    elsif ( $this_arg eq '-head') {$header = 1;}
    elsif ( $this_arg eq '-noXY') {$noXY = 1;}
    elsif ( $this_arg eq '-o') {$ResFilename = shift @ARGV;}
    elsif ( $this_arg eq '-reg') {$regFilename = shift @ARGV;}   
    elsif ( $this_arg eq '-minScore') {$minscore = shift @ARGV;}   
    elsif ( $this_arg =~ m/^-/ ) { print "unknown flag: $this_arg\n";}
}

if ($ResFilename eq "") {
    $ResFilename = $SitesFilename.".withGenes.txt";
}



#------------read expression/foldChange annotation of genes---------

my %geneReg;
if ( $regFilename eq ""){
    print "you did not specify file with expression of fold change values\n";
}
else {
	open (REG, "<$regFilename ") or die "Cannot open file $regFilename !!!!: $!";
	#Gene	median	mad	flag
	#RPL41	13.60755187	0.074233841	EXPRESSED
        #A	1110013L07Rik	0.26	up-regulated				

	while (<REG>) {
		chomp;	
		my ($name, $value) = (split)[1,3];
		$geneReg{$name} = $value;
	}
	close REG;
}
#-----------read genes----------------

my %genes;

my $count = 0;

open (GENES, "<$GenesFilename") or die "Cannot open file $GenesFilename!!!!: $!";
<GENES>;

my ($name,$chr,$strand,$leftPos,$rightPos);

while (<GENES>) {
	chomp;
	if (/(chr\S*)\s(\d+)\s(\d+)\s([-]?1)\s\S+\s\S+\s(\S+)/){
		$name = $5;
		$chr = $1;	
		
		 $strand = $4;
		if ($strand eq '1') {
			$strand = 1;
		}
		else {
			$strand = -1;	
		}
		 $leftPos = $2;
		 $rightPos = $3;	
		next if (($chr =~ m/[XY]/)&&($noXY));
		next if ($name =~ m/MIR/);		
	} elsif (/(chr\S*)\s([+-])\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\d+)\s(\S+)\s(\S+)\s\S+\s(\S+)/) {
		$name = $10;
		$chr = $1;
		$strand = $2;
		$leftPos = $3;
		$rightPos = $4;
		next if (($chr =~ m/[XY]/)&&($noXY));
		next if ($name =~ m/MIR/);		
		
		if ($strand eq '+') {
			$strand = 1;
		}
		else {
			$strand = -1;	
		}	
	} else {
	    print "Wrong type: ",$_,"\n";
	}
	my $ID = "$name\t$chr:$leftPos-$rightPos";					
	my $reg = "NA";	
	if (exists($geneReg{$name})) {
		$reg = $geneReg{$name};
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
				
	$genes{$chr}->{$ID}{'strand'} = $strand; 
				
	$genes{$chr}->{$ID}{'TSS'} = ($strand == 1) ? $leftPos :$rightPos ; 
	$genes{$chr}->{$ID}{'TE'} = ($strand == -1) ? $leftPos :$rightPos ;		
	$genes{$chr}->{$ID}{'reg'} = $reg;
	$genes{$chr}->{$ID}{'length'} = abs ($leftPos-$rightPos); 
	$genes{$chr}->{$ID}{'closestPicDist'} = $BIGNUMBER;
	$genes{$chr}->{$ID}{'closestPositivePicDist'} = $BIGNUMBER;
	$genes{$chr}->{$ID}{'closestNegativePicDist'} = -$BIGNUMBER;

	$genes{$chr}->{$ID}{'closestPicDistTE'} = $BIGNUMBER;	
	
}
print "Total genes : $count\n";



close GENES;
print "\t\t$GenesFilename is read!\n";
#for my $gName (sort keys %{$genes{'chr18'}}) {

 #    print "$gName\t$genes{'chr18'}->{$gName}{'TSS'}\n";
#}

#-----------read file with sites miRNA, store as genes-----

if ( $MirFilename eq ""){
    print "you did not specify file with miRNA\n";
}
else {


	open (MIR, "<$MirFilename ") or die "Cannot open file $MirFilename !!!!: $!";
	#chr1	20669090	20669163	mmu-mir-206	960	+
	while (<MIR>) {
		chomp;	
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
		($genes{$chr}->{$ID}{'firstIntronStart'},$genes{$chr}->{$ID}{'firstIntronEnd'}) = (0,0);	
		$genes{$chr}->{$ID}{'closestPicDist'} = $BIGNUMBER;
		$genes{$chr}->{$ID}{'closestPicDistTE'} = $BIGNUMBER;
		$genes{$chr}->{$ID}{'closestPositivePicDist'} = $BIGNUMBER;
		$genes{$chr}->{$ID}{'closestNegativePicDist'} = -$BIGNUMBER;	

		
	}


	close MIR;
	print "\t\t$MirFilename is read!\n" if ($verbose) ;
}



#-----------read file with sites and find overlapping genes-----
my $numberOfAllSites = 0;

my $lengthLosses = 0;
my $lengthGains = 0;

open (FILE, "<$SitesFilename") or die "Cannot open file $SitesFilename!!!!: $!";
open (OUT , ">$ResFilename") or die "Cannot open file $ResFilename!!!!: $!";
print OUT "chrom	start	end	max_coord	score	Dist	DistTE	geneName	geneCoord	Reg\n";
open (GENESOUT , ">$ResFilename.genes") or die "Cannot open file $ResFilename.genes!!!!: $!";

print GENESOUT "Name\tCoord\tDist\tDistTE\tReg\n";

my $correction = 0;

if ($header) {
	   my $headerString = <FILE>; #read header
	   my @a = split (/\t/,$header );
	    
	    if ($a[1] =~m/chr/) {
		    $correction=1;
	    }
}

$count = 0;
print $SitesFilename, "\n";
while (<FILE>) {
	chomp; 
	next if (/track/);
	my $string = $_; 
	chomp $string; 
	$numberOfAllSites++;	
	my @a = split /\t/;
	if (($a[1] =~m/chr/)&&($correction==0)) {
	    $correction=1;
	}
	my $chro = $a[0+$correction];
	next if (($chro =~ m/[XY]/)&&($noXY));
		unless ($chro =~ m/chr/) {
			$chro = "chr".$chro;
		}

	my $fpos = $a[1+$correction];
	my $lpos = $a[2+$correction];
	my $score = $a[4+$correction];
	if ($minscore>0) {
		next if ($score<$minscore);
	}
	next if ($chro=~m/M/);
	
	print OUT $string, "\t"; 
	
	findClosestTSS ($fpos,$lpos,$chro,\%genes);
	
	$count ++;	
}
close FILE;
close OUT ;
close GENESOUT;
print "Total Sites: $count\n";

open (GENESOUT , ">$ResFilename.genes.ClosestPeakDist") or die "Cannot open file $ResFilename.genes.ClosestPeakDist!!!!: $!";
print GENESOUT "Name	Coord	Dist	DistTE	Reg	posDist	negDist\n";
for my $chro (keys %genes) {
	for my $gName (keys %{$genes{$chro}}) {
		print GENESOUT "$gName\t",$genes{$chro}->{$gName}{'closestPicDist'},"\t",$genes{$chro}->{$gName}{'closestPicDistTE'},"\t",$genes{$chro}->{$gName}{'reg'},"\t",$genes{$chro}->{$gName}{'closestPositivePicDist'},"\t",$genes{$chro}->{$gName}{'closestNegativePicDist'},"\n";
	}
}
close GENESOUT; 

#my @arr = @{$genes{'chr'}};

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

sub overlap {
	my ($a1,$a2,$e1,$e2)= @_;	#e for genes
	
	if ($a1 > $a2) {
	    ($a1,$a2) = ($a2,$a1);
	}

	if ($e1 > $e2) {
	    ($e1,$e2) = ($e2,$e1);
	}

	
	if (($a1 >= $e1)&&($a1 <= $e2)) {
		return 2;
	}
	if (($a2 >= $e1)&&($a2 <= $e2)) {
		return 2;
	}
	if (($e1 >= $a1)&&($e2 <= $a2)) {
		return 1;
	}

	return 0;
}


sub findClosestTSS {
    my ($fpos,$lpos,$chro,$genes) = @_;
    my $pos = ($fpos+$lpos)/2;
    my $gene =  "";
    my $tssDist=$BIGNUMBER; 
    my $teDist=$BIGNUMBER; 

    for my $gName (keys %{$genes->{$chro}}) {
	    my $TSS = $genes->{$chro}{$gName}{'TSS'};	
	    #print  $gName,":",$TSS," " if ($chro eq "chr15"); 
	    if (abs($pos-$TSS)<abs($tssDist)) {
		$gene=$gName;
		$tssDist = ($pos-$TSS)*$genes->{$chro}{$gName}{'strand'};
		$teDist = ($pos-$genes->{$chro}{$gName}{'TE'})*$genes->{$chro}{$gName}{'strand'};
	    }	
	    if (abs($pos-$TSS)<$BIGNUMBER) {
		print GENESOUT "$gName\t",($pos-$TSS)*$genes->{$chro}{$gName}{'strand'},"\t",($pos-$genes->{$chro}{$gName}{'TE'})*$genes->{$chro}{$gName}{'strand'},"\t",$genes->{$chro}{$gName}{'reg'},"\n";

		if (abs($pos-$TSS)<abs($genes->{$chro}{$gName}{'closestPicDist'})) {
			$genes->{$chro}{$gName}{'closestPicDist'} = ($pos-$TSS)*$genes->{$chro}{$gName}{'strand'};
			$genes->{$chro}{$gName}{'closestPicDistTE'} = ($pos-$genes->{$chro}{$gName}{'TE'})*$genes->{$chro}{$gName}{'strand'};
		}
		if (($pos-$TSS)*$genes->{$chro}{$gName}{'strand'}>=0 && abs($pos-$TSS)<$genes{$chro}->{$gName}{'closestPositivePicDist'} ) {
			$genes->{$chro}{$gName}{'closestPositivePicDist'} = ($pos-$TSS)*$genes->{$chro}{$gName}{'strand'};
		}
		if (($pos-$TSS)*$genes->{$chro}{$gName}{'strand'}<=0 && -1*abs($pos-$TSS)>$genes{$chro}->{$gName}{'closestNegativePicDist'} ) {
			$genes->{$chro}{$gName}{'closestNegativePicDist'} = ($pos-$TSS)*$genes->{$chro}{$gName}{'strand'};
		}
		
	    }	
    }
    print OUT "$tssDist\t$teDist\t$gene\t$genes{$chro}->{$gene}{'reg'}\n" if ($gene ne "");  
    if ($gene eq "") {
	 print OUT "$tssDist	$teDist	NA	NA	NA\n";
    } 

}

sub printOvelapingGenes {	
	my ($fpos,$lpos,$chro,$genes,$globalH,$type) = @_;	
	#print $fpos,$lpos,$chro,$genes,$globalH,$type,"\n";
	my %hash1;
	my %hash2;	
	my @array2return;
	for my $gName (keys %{$genes->{$chro}}) {
		my $TSS = $genes->{$chro}{$gName}{'TSS'};	
	    	my $TE = $genes->{$chro}{$gName}{'TE'};	
		
		if (overlap($fpos,$lpos,$TSS,$TE)==1) {
		    $hash1{$gName} = overlap($fpos,$lpos,$TSS,$TE);
		}
		if (overlap($fpos,$lpos,$TSS,$TE)==2) {
		    $hash2{$gName} = overlap($fpos,$lpos,$TSS,$TE);
		} 
	}
	print OUT join(",", sort keys %hash1),"\t",join(",", sort keys %hash2),"\n";
	if ($type eq "loss") {
	    for my $key (sort keys %hash1) {
		$globalH->{$key} = 1;
	    }
	    for my $key (sort keys %hash2) {
		$globalH->{$key} = 1;
	    }	    
	}else {
	    for my $key (sort keys %hash1) {
		$globalH->{$key} = 1;
	    }	    
	}
	
}

sub min {
	return $_[0] if ($_[0]<$_[1]);
	$_[1];
}
