#!/usr/bin/perl 
use strict;

my ($infile,$annofile,$outfile)=@ARGV;

my %tx2gene;
open(IN,$annofile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	$tx2gene{$array[0]}=$array[1];
}
close IN;

open(IN,$infile) || die $!;
open(OUT,">$outfile") || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	
	unless($_=~/^ENST/) {
		print OUT $_,"\n";
	}
	else {
		print OUT join("\t",$tx2gene{$array[0]},@array[1..$#array]),"\n";
	}
}
close IN;
close OUT;
