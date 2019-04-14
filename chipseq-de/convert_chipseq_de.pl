#!/usr/bin/perl -w
use strict;

#convert chipseq peak de to gene centric result

#retrieve only promoter for comparison


my ($infile,$outfile) = @ARGV;

my $cate="promoter-TSS";

my %gene2info;
my $title;
my $linenum=0;

open(IN,$infile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	
	if($linenum==0) {
		$title="Gene\t".join("\t",@array[1..6]);
	}
	else {
		#25, 27 #cate, gene
		if($array[25] eq $cate) {
			$gene2info{$array[27]}{join("\t",@array[1..6])}++;
		}
	}
	$linenum++;
}
close IN;


open(OUT,">$outfile") || die $!;
print OUT $title,"\n";
foreach my $gene (sort keys %gene2info) {
	print OUT $gene,"\t";
	
	my $minp=1.1;
	my $best;
	foreach my $info (sort keys %{$gene2info{$gene}}) {
		my @nums=split("\t",$info);
		if($nums[4] eq "NA") {
			$nums[4]=1;
		}
		if($nums[4]<$minp) {
			$best=$info;
		}
	}
	
	print OUT $best,"\n";
}
close OUT;
	

