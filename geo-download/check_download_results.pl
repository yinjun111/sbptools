#!/usr/bin/perl -w
use strict;


#check fastq files by name, and match to download list, generate new sh for files not downloaded
my ($downloadfolder,$downloadscript, $newscript)=@ARGV;

my @fastqfiles=glob("$downloadfolder/*.fastq.gz");

my %srrs;
foreach my $fastqfile (@fastqfiles) {
	if($fastqfile=~/(SRR\d+)/) {
		$srrs{$1}++;
	}
}

open(IN,$downloadscript) || die $!;
open(OUT,">$newscript") || die $!;
while(<IN>) {
	tr/\r\n//d;
	if($_=~/(SRR\d+)/) {
		unless(defined $srrs{$1}) {
			print OUT $_,"\n";
		}
	}
}
close IN;
close OUT;
