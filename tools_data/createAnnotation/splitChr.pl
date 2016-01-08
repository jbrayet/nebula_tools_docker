#!/usr/bin/perl

$f = $ARGV[0]; #get the file name
$p = $ARGV[1];
$o = $ARGV[2];

open (INFILE, "<$p/$f")
or die "Can't open: $p/$f $!";

while (<INFILE>) {
$line = $_;
chomp $line;
if ($line =~ /\>/) { #if has fasta >
close OUTFILE;
$new_file = substr($line,1);
@new_file = split(/ /,$new_file);
$new_file = $new_file[0];
$new_file .= ".fa";
open (OUTFILE, ">$o/$new_file")
or die "Can't open: $o/$new_file $!";
}
print OUTFILE "$line\n";
}
close OUTFILE;

