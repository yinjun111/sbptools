#!/usr/bin/perl -w
use strict;

#read DM peaks, convert to bed format

my ($infile,$outfile)=@ARGV;

open(IN,$infile) || die $!;
open(OUT,">$outfile") || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	if($array[5] eq "1" || $array[5] eq "-1") {
		if($array[0]=~/(\w+)\:(\d+)\-(\d+)/) {
			print OUT join("\t",$1,$2,$3,$array[0]),"\n";
		}
	}
}
close IN;
close OUT;

