#!/usr/bin/perl -w
use strict;



########
#Prerequisites
########



########
#Interface
########


my $version="0.1";


my $usage="

rename_bed
version: $version
Usage: perl rename_bed.pl file.bed file_rename.bed

Description: select peaks merged from at least two different samples

Parameters:	
	
";


unless (@ARGV) {
	print STDERR $usage;
	exit;
}

my $params=join(" ",@ARGV);
#then call different scripts


########
#Parameters
########





########
#Process
########

my ($infile,$outfile,$minsamplenum)=@ARGV;


if(!(defined $minsamplenum)) {
	$minsamplenum=2;
}


open(IN,$infile) || die $!;
open(OUT,">$outfile") || die $!;

while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	
	my %samples;
	
	#peaks are separated by ","
	my @peaks=split(",",$array[3]);
	
	foreach my $peak (@peaks) {
		my $sample=(split("\\\|",$peak))[1];
		$samples{$sample}++;
	}
	
	if(keys %samples>=$minsamplenum) {
		#rename the name of the peak to coord
		print OUT join("\t",$array[0],$array[1],$array[2],$array[0].":".($array[1]+1)."-".$array[2]),"\n";
	}
}

close IN;
close OUT;
		