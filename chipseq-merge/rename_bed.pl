#!/usr/bin/perl -w
use strict;
use Getopt::Long;


########
#Prerequisites
########



########
#Interface
########


my $version="0.2";


my $usage="

rename_bed
version: $version
Usage: perl rename_bed.pl -i file.bed -o file_rename.bed -s samplename

Description: Prepare bed files for merging

Parameters:	
    --in|-i           Input file
    --out|-o          Output file
    --bysample|-s     Rename by sample name
    --bycoord|-c      Rename by coordinate
	
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

my ($infile,$outfile);
my $bysample;
my $bycoord;


GetOptions(
	"in|i=s" => \$infile,
	"out|o=s" => \$outfile,
	"bysample|s=s" => \$bysample,
	"bycoord|c" => \$bycoord,	
);


unless((defined $bysample && length($bysample)>0) || $bycoord) {
	print STDERR "ERROR:either --bysample or --bycoord needs to be defined.\n\n";
	exit;
}

########
#Process
#######

open(IN,$infile) || die $!;
open(OUT,">$outfile") || die $!;
while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^#/;
	my @array=split/\t/;
	if(defined $bysample && length($bysample)>0) {
		if(@array==4) {
			print OUT join("\t",@array[0..2],$array[3]."|".$bysample),"\n";
		}
		elsif(@array<4) {
			print OUT join("\t",$array[0],$array[1],$array[2],$array[0].":".($array[1]+1)."-".$array[2]."|".$bysample),"\n";
		}
		else {
			print OUT join("\t",@array[0..2],$array[3]."|".$bysample,@array[4..$#array]),"\n";
		}
	}
	elsif($bycoord) {
		if(@array<=4) {
			print OUT join("\t",$array[0],$array[1],$array[2],$array[0].":".($array[1]+1)."-".$array[2]),"\n";
		}
		else {
			print OUT join("\t",$array[0],$array[1],$array[2],$array[0].":".($array[1]+1)."-".$array[2],@array[4..$#array]),"\n";
		}
	}
}
close IN;
close OUT
	
