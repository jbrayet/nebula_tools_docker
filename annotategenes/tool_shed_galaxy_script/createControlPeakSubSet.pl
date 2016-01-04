#!/usr/bin/perl

#create a control dataset with the same number of reads as in the SAMPLE (highest peaks)

use strict;
use warnings;
use diagnostics;

my $usage = qq{
    $0

    -----------------------------
    mandatory parameters:
    
    -f CHiP_file               	
    -c control_file
    -o output file
    -----------------------------
    optional parameters:

    -n number of files to create
    
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

my $nBootstrap = 1;


## optional arguments

## parse command line arguments

while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;
    if ( $this_arg eq '-h') {print "$usage\n"; exit; }

    elsif ( $this_arg eq '-f') {$filename = shift @ARGV;}
    elsif ( $this_arg eq '-c') {$controlFilename = shift @ARGV;}
    elsif ( $this_arg eq '-o') {$output_fname = shift @ARGV;}
    elsif ( $this_arg eq '-n') {$nBootstrap = shift @ARGV;}

    elsif ( $this_arg =~ m/^-/ ) { print "unknown flag: $this_arg\n";}
}

if ( $filename eq ""){
    die "you should specify chip file\n";
}
if( $controlFilename eq ""){
    die "you should specify control file\n";
}

if( $output_fname eq ""){
    die "you should specify output filename\n";
}


print "\n-----------------\n\n";

my %hash;
my $chipCount = 0;
my @header;
open FILE, "< $filename " || die "$filename : $!\n";    
while(<FILE>){
    $chipCount++;	    
}
close FILE;
#print "ChIP: $chipCount\n";


my $controlCount = 0;

open FILE, "< $controlFilename " || die "$controlFilename : $!\n";    
while(<FILE>){
    next if (/track/);
    my $entry = $_;
    my @fields = split(/\t/,$_);
    $hash{$entry} = $fields[4]; 
    $controlCount++;
}
#print "controlCount : $controlCount\n";

close FILE;
open OUT, "> $output_fname" || die "$output_fname: $!\n";
my $count = 0;
if ($controlCount>$chipCount) {
    my $prob = $chipCount/$controlCount*1.1;
    for my $entry (sort {$hash{$b}<=>$hash{$a}} keys %hash) {	
	my $yes = rand();
	if ($yes<=$prob) {$yes=1;}else {$yes=0;}	
	if ($yes) {
		print OUT $entry ;	
		$count++;	
	}
	if ($count >=$chipCount) {
	    last;
	}
    }
   
} else {    
    for my $entry (keys %hash) {	
	print OUT $entry;		
    }
}

close OUT;


for my $try (2..$nBootstrap) {

	open OUT, "> $output_fname$try" || die "$output_fname$try: $!\n";
	my $count = 0;
	if ($controlCount>$chipCount) {
	    my $prob = $chipCount/$controlCount*1.1;
	    for my $entry (sort {$hash{$b}<=>$hash{$a}} keys %hash) {	
		my $yes = rand();
		if ($yes<=$prob) {$yes=1;}else {$yes=0;}		
		if ($yes) {
			print OUT $entry ;	
			$count++;	
		}
		
		if ($count >=$chipCount ) {
		    last;
		}
	    }
	    if ($count <$chipCount) { #do it again!

		for my $entry (sort {$hash{$b}<=>$hash{$a}} keys %hash) {	
			my $yes = rand();
			if ($yes<=$prob) {$yes=1;}else {$yes=0;}		
			if ($yes) {
				print OUT $entry ;	
				$count++;	
			}
		
			if ($count >=$chipCount ) {
			    last;
			}
	        }

	    }
	   
	} else {    
	    for my $entry (keys %hash) {	
		print OUT $entry;		
	    }
	}

	close OUT;

}