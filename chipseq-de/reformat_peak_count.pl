#!/usr/bin/perl -w
use strict;

my $version=0.1;

#separate annotation columns
#rename sample columns

#need to add config file to change sample order ....!!!

my ($infile,$outfile,$annofile)=@ARGV;

open(IN,$infile) || die $!;
open(OUT1,">$outfile") || die $!;
open(OUT2,">$annofile") || die $!;

my $linenum=0;
my @samples;
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
		
		#data starts from column 20		
		print OUT1 join("\t","PeakID",@samples),"\n";
		print OUT2 join("\t","PeakID",@array[1..18]),"\n";
		
	}
	else {
		print OUT1 join("\t",@array[0,19..$#array]),"\n";
		print OUT2 join("\t",@array[0..18]),"\n";
	}
	$linenum++;
}
close IN;
close OUT1;
close OUT2;

	