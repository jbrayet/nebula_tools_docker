#!/usr/bin/perl

#filter out dulpicates from SAMPLE (optional) and create a control dataset w/o duplicates with the same number of reads as in the SAMPLE

use strict;
use warnings;
use diagnostics;

my $usage = qq{
    $0

    -----------------------------
    mandatory parameters:
    
    -f peaks              	
    -t min peak height
    -o output file

    -----------------------------
    optional parameters:
    -n name   
    none			
};

if(scalar(@ARGV) == 0){
    print $usage;
    exit(0);
}

## mandatory arguments

my $filename = "";
my $output_fname = "";

my $minPeakH = 0;

my $chromLengthsFile="";
my $expName = "User Track";

## optional arguments

## parse command line arguments

while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;
    if ( $this_arg eq '-h') {print "$usage\n"; exit; }

    elsif ( $this_arg eq '-f') {$filename = shift @ARGV;}
    elsif ( $this_arg eq '-t') {$minPeakH = shift @ARGV;}
    elsif ( $this_arg eq '-o') {$output_fname = shift @ARGV;}
    elsif ( $this_arg eq '-g') {$chromLengthsFile = shift @ARGV;}
    elsif ( $this_arg eq '-n') {$expName = shift @ARGV;}
    elsif ( $this_arg =~ m/^-/ ) { print "unknown flag: $this_arg\n";}
}

if ( $filename eq ""){
    die "you should specify chip file\n";
}
if( $output_fname eq ""){
    die "you should specify output filename\n";
}

$minPeakH-=0.5 unless ($minPeakH=~m/0\.5/);

#read chromosome lengths if provided:
my %max;
if ($chromLengthsFile ne "") {
	open FILE, "< $chromLengthsFile " || die "$chromLengthsFile : $!\n";
	while(<FILE>){
    		chomp;		
		if (/(chr\S+)\s(\d+)/) {
			$max{$1}=$2;
		}
	}
	close FILE;
}

######

print "\n-----------------\n\n";

my %hash;
my $chipCount = 0;
my @header;


open FILE, "< $filename " || die "$filename : $!\n";
open OUT, "> $output_fname" || die "$output_fname: $!\n";
print OUT "track name=\'$expName\' description=\'$expName\'\n";
my $count = 0;
my $scount = 0;

while(<FILE>){
    chomp;
    next if (/max/);
    next if (/track/);
    next if (/^\#/);
    my @fields = split(/\t/,$_);
    my $entry = $fields[0]."\t".$fields[2]."\t".$fields[3];
    $count++; 
    if ($fields[4]>=$minPeakH) {		    
	$scount ++;
	$fields[0]= "chr".$fields[0] unless ($fields[0]=~m/chr/);

	if ($chromLengthsFile ne "") {
		my $maxV = $max{$fields[0]}; 
		$fields[2]= min($fields[2],$maxV);
		$fields[3]=min($fields[3],$maxV);
		$fields[1]=min($fields[1],$maxV); 
	}
	print OUT join("\t",$fields[0],$fields[2],$fields[3],$fields[1],$fields[4],"+",$fields[2],$fields[3],"255,120,11","1",$fields[3]-$fields[2],0,"\n");
    }
}
 
close FILE;
close OUT;
print "read: $count peaks; selected: $scount\n";

sub min {
	my ($a,$b) = @_;
	if($a<$b) {
		return $a;
	}
	$b;
}

