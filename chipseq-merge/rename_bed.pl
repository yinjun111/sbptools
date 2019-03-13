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

Description: Prepare bed files for merging

Parameters:	
	
";


unless (@ARGV) {
	print STDERR $usage;
	exit;
}

my $params=join(" ",@ARGV);
#then call different scripts


########
#Process
########

my ($infile,$outfile)=@ARGV;

my $samplename;

#bed files are in standard file names
if($infile=~/([^\/]+)_Peaks_sorted.bed/) {
	$samplename=$1;
}
else {
	print STDERR "ERROR:$infile name is wrong.\n";
}

open(IN,$infile) || die $!;
open(OUT,">$outfile") || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	print OUT join("\t",@array[0..2],$array[3]."|".$samplename,@array[4..$#array]),"\n";
}
close IN;
close OUT
	
