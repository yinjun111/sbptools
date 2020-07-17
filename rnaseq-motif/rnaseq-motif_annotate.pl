#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);

#annotate gene promoter binding using pre-computed calculations

########
#Updates
########



########
#Prerequisites
########

#none


########
#Interface
########


my $version="0.1";


my $usage="

rnaseq-motif_annotate
version: $version
Usage: 

Description: 

Parameters:
    --in|-i           Input file, a bed file with name as symbol|transcript
    --out|-o          Output file, an annotated file
    --anno|-a         Annotation file
	
	
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

my $infile;
my $outfile;
my $annofile;

my $verbose=1;


GetOptions(
	"in|i=s" => \$infile,
	"out|o=s" => \$outfile,
	"anno|a=s" => \$annofile,

	"verbose" => \$verbose,
	"help|h" => sub {print STDERR $usage;exit;}
	
);


#my $logfile="rnaseq-var_filter_run.log";


######
#Process input file
######

print STDERR "\nrnaseq-motif_anno version $version running...\n";


#precomputed annotation
my %tx2anno;
my $title;

open(IN,$annofile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	
	if($_=~/^Gene/) {
		$title=$_;
	}
	else {
		$tx2anno{$array[1]."|".$array[2]}=$_;
	}
}
close IN;


#
open(IN,$infile) || die $!;
open(OUT,">$outfile") || die $!;

print OUT $title,"\n";

while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
		
	if(defined $tx2anno{$array[3]}) {
		print OUT $array[3],"\t",$tx2anno{$array[3]},"\n";
	}
	else {
		print STDERR $array[3]," not defined.\n";
	}
}
close IN;
close OUT;
