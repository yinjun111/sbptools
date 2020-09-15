#!/usr/bin/perl -w
use strict;
use Getopt::Long;


my $version=0.2;

#separate annotation columns
#rename sample columns

#v0.2 add config file to change sample order ....!!!



my $usage="

reformat_peak_count
version: $version
Usage: reformat_peak_count.pl [parameters]

Description: 


Mandatory Parameters:
    --in|-i           Input file
    --out|-o          Output file
    --config|-c       Configuration file to match the sample order
    --anno|-a         Annotation file for peaks

";



unless (@ARGV) {
	print STDERR $usage;
	exit;
}

my $params=join(" ",@ARGV);
#then call different scripts



#####
#parameters
#####

my ($infile,$outfile,$annofile,$configfile);

GetOptions(
	"in|i=s" => \$infile,
	"config|c=s" => \$configfile,
	"out|o=s" => \$outfile,
	"anno|a=s" => \$annofile,
);


#####
# Program running
#####

my %sample2order;

my $linenum=0;
open(IN,$configfile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	
	unless ($_=~/^Sample/) {
		$sample2order{$array[0]}=$linenum;
	}
	
	$linenum++;
}
close IN;


open(IN,$infile) || die $!;
open(OUT1,">$outfile") || die $!;
open(OUT2,">$annofile") || die $!;

$linenum=0;
my @samples;
my $sample;
my @sampleorders;
my %sample2info;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;

	if($linenum==0) {
		#/data/jyin/XuLab/Processed/ATAC-Seq/Aligned/A528TAb_1/A528TAb_1_TagDir Tag Count in given bp (63065796.5 Total, normalization factor = 1, effective total = 10000000)
		foreach my $colname (@array[19..$#array]) {
			if($colname=~/([^\/]+)_TagDir /) {
				push @samples,$1;
				
			}
		}
		
		@sampleorders=cal_order(map {$sample2order{$_} }@samples);
		
		#data starts from column 20		
		#print sample according to config order
		print OUT1 join("\t","PeakID",@samples[@sampleorders]),"\n";
		print OUT2 join("\t","PeakID",@array[1..18]),"\n";
		
	}
	else {
		print OUT1 join("\t",@array[0,(19..$#array)[@sampleorders]]),"\n";
		
		print OUT2 join("\t",@array[0..18]),"\n";
	}
	$linenum++;
}
close IN;
close OUT1;
close OUT2;


####
#Functions
####

#return order of samples
sub cal_order {
	my @nums=@_;
	
	#give order of numbers
	my %num2pos;
	
	for(my $n=0;$n<@nums;$n++) {
		$num2pos{$nums[$n]}=$n;
	}
	
	return map {$num2pos{$_}} sort {$a<=>$b} keys %num2pos;
}

