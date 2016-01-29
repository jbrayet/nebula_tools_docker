#!/usr/bin/perl

#filter out dulpicates from SAMPLE (optional) and create a control dataset w/o duplicates with the same number of reads as in the SAMPLE

use strict;
use warnings;
use diagnostics;

my $usage = qq{
    $0

    -----------------------------
    mandatory parameters:
    
    -f CHiP_file               	
    -c control_file
    -t type [bam, sam, eland]
    -o output file
    -----------------------------
    optional parameters:
    
    none			
};

if(scalar(@ARGV) == 0){
    print $usage;
    exit(0);
}

## mandatory arguments

my $filename = "";
my $output_fname = "";

my $controlFilename = "";
my $type = "";
my $sampleOutput = "";


## optional arguments

## parse command line arguments

while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;
    if ( $this_arg eq '-h') {print "$usage\n"; exit; }

    elsif ( $this_arg eq '-f') {$filename = shift @ARGV;}
    elsif ( $this_arg eq '-c') {$controlFilename = shift @ARGV;}
    elsif ( $this_arg eq '-t') {$type = shift @ARGV;}
    elsif ( $this_arg eq '-o') {$output_fname = shift @ARGV;}
    elsif ( $this_arg eq '-s') {$sampleOutput = shift @ARGV;}


    elsif ( $this_arg =~ m/^-/ ) { print "unknown flag: $this_arg\n";}
}

if ( $filename eq ""){
    die "you should specify chip file\n";
}
if( $controlFilename eq ""){
    die "you should specify control file\n";
}
if( $type eq ""){
    die "you should specify file type (bam, sam or eland)\n";
}
if( $output_fname eq ""){
    die "you should specify output filename\n";
}


print "\n-----------------\n\n";

my %hash;
my $chipCount = 0;
my @header;


if ($type eq "eland") {
    open FILE, "< $filename " || die "$filename : $!\n";    
    while(<FILE>){
	    my @fields = split(/\t/,$_);
	    my $entry = $fields[6].":".$fields[7]."-".$fields[8];
	    unless (exists($hash{$entry})) {
		    $hash{$entry} = $_;
		    $chipCount++;
	    }
    }
} elsif ($type eq "sam") {
    open FILE, "< $filename " || die "$filename : $!\n";    
    while(<FILE>){
	    if (m/^@/) {		
		push(@header,$_);
		next;
	    }
	    my @fields = split(/\t/,$_);
	    next if (scalar(@fields)<10);
	    my $entry = $fields[2].":".$fields[3]."-".$fields[1];
	    unless (exists($hash{$entry})) {
		    $hash{$entry} = $_;
		    $chipCount++;
	    }
    }
} elsif ($type eq "bam") {
    open(FILE, "samtools view -h $filename |") or die "$0: can't open ".$filename.":$!\n";
    while(<FILE>){
	    if (m/^@/) {		
		push(@header,$_);
		next;
	    }
	    my @fields = split(/\t/,$_);
	    next if (scalar(@fields)<10);
	    my $entry = $fields[2].":".$fields[3]."-".$fields[1];
	    unless (exists($hash{$entry})) {
		    $hash{$entry} = $_;
		    $chipCount++;
	    }
    }
}
close FILE;
print "ChIP: $chipCount\n";

if ($sampleOutput ne "") {
    
    open OUT, "> $sampleOutput" || die "$sampleOutput: $!\n";

    if ($type eq "bam" || $type eq "sam") { #print header
	for my $headerLine (@header) {
	    print OUT $headerLine;
	}
    }
    for my $line (values %hash) {
	print OUT $line;	  
    }    
    close OUT;
}

delete @hash{keys %hash};
@header = ();

my $controlCount = 0;
if ($type eq "eland") {
    open FILE, "< $controlFilename " || die "$controlFilename : $!\n";    
    while(<FILE>){
	    my @fields = split(/\t/,$_);
	    my $entry = $fields[6].":".$fields[7]."-".$fields[8];
	    unless (exists($hash{$entry})) {
		    $hash{$entry} = $_;
		    $controlCount++;
	    }
    }
} elsif ($type eq "sam") {
    open FILE, "< $controlFilename " || die "$controlFilename : $!\n";    
    while(<FILE>){
	    if (m/^@/) {		
		push(@header,$_);
		next;
	    }
	    my @fields = split(/\t/,$_);
	    my $entry = $fields[2].":".$fields[3]."-".$fields[1];
	    unless (exists($hash{$entry})) {
		    $hash{$entry} = $_;
		    $controlCount++;
	    }
    }
} elsif ($type eq "bam") {
    open(FILE, "samtools view -h $controlFilename |") or die "$0: can't open ".$controlFilename.":$!\n";
    while(<FILE>){
	    if (m/^@/) {		
		push(@header,$_);
		next;
	    }
	    my @fields = split(/\t/,$_);
	    my $entry = $fields[2].":".$fields[3]."-".$fields[1];
	    unless (exists($hash{$entry})) {
		    $hash{$entry} = $_;
		    $controlCount++;
	    }
    }
}
close FILE;
print "Control: $controlCount\n";
my $prob = $chipCount/$controlCount;

open OUT, "> $output_fname" || die "$output_fname: $!\n";

if ($type eq "bam" || $type eq "sam") { #print header
    for my $headerLine (@header) {
	print OUT $headerLine;
    }
} 
my $count = 0;

for my $line (values %hash) {
	my $rand = rand();
	
	if ($rand < $prob) {
		print OUT $line;
		$count ++;	
	}
	last if ($count == $chipCount);	   
}


if ($count < $chipCount) {

	$prob = ($chipCount-$count)/$controlCount*1.1;

	for my $line (values %hash) {
		my $rand = rand();
	
		if ($rand < $prob) {
			print OUT $line;
			$count ++;	
		}
		last if ($count == $chipCount);	  
	} 	
}

print "count = $count\n";
close OUT;

