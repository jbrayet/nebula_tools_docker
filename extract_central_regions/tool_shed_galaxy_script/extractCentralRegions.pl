#!/usr/bin/perl -w
use strict;

#creats a BED file with central area of peaks


my $usage = qq{
    $0

    -----------------------------
    mandatory parameters:
    
    -f   filename 	file with sites in BED format
    -l   value	 	length of the cental regions
    
    -----------------------------
    optional parameters:
    -v				for verbose   
    -head	1/0		if there is a header
    -o		filename	output file
};

if(scalar(@ARGV) <2){
    print $usage;
    exit(0);
}

my $flank = 0;
my $ResFilename = "";
my $file = "";
my $header = 0;
my $verbose = 0;

while(scalar(@ARGV) > 0){
    my $this_arg = shift @ARGV;
    if ( $this_arg eq '-h') {print "$usage\n"; exit; }
    elsif ( $this_arg eq '-f') {$file = shift @ARGV;}
    elsif ( $this_arg eq '-v') {$verbose = 1;}  
    elsif ( $this_arg eq '-head') {$header = shift @ARGV;}
    elsif ( $this_arg eq '-o') {$ResFilename = shift @ARGV;}
    elsif ( $this_arg eq '-l') {$flank = shift @ARGV;$flank /=2;}   
    elsif ( $this_arg =~ m/^-/ ) { print "unknown flag: $this_arg\n";}
}




my $count = 0;
my %hash;

open (FILE, "<$file") or die "Cannot open file $file!!!!: $!";
open (OUT, ">$ResFilename") or die "Cannot open file $ResFilename!!!!: $!";

if ($header) {
<FILE>;
}

while (<FILE>) {		
	s/\R//g; 
	next if (/^#/);
	next if (/track/);
	next if (/summit/);
        my @a = split /\s/;
        my $chr = $a[0];
        my $maxPos =  $a[3];
        my $score = $a[4];
	if ($maxPos=~/\D/) {
		$maxPos = int(($a[1]+$a[2])/2);
	} elsif ($maxPos < $a[1]){ #MACS intervals
		$maxPos = $a[1]+$a[4];
		$score = $a[5];
	}
        my $firstPos = $maxPos-$flank;
        my $lastPos = $maxPos+$flank;

        my $ID=$chr.":".$firstPos."_".$lastPos."_".$score ;
        unless (exists($hash{$ID})) {
            $hash{$ID}=1;
            $count++;            
            print OUT "$chr\t$firstPos\t$lastPos\t$score\n";
        }
}
    
#print "$file\t$count\n";
close FILE;
close OUT;

