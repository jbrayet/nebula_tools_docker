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
    for my $entry (sort {$hash{$b}<=>$hash{$a}} keys %hash) {	
	print OUT $entry;
	$count++;
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
