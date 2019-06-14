#!/usr/bin/perl -w
use strict;

my ($infile,$outfile) = @ARGV;


my %peak2motif;


open(IN,$infile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	$peak2motif{$array[3]}{$array[7]}++;
}
close IN;

open(OUT,">$outfile") || die $!;
print OUT "Peak\tHomer Motifs\n";

foreach my $peak (sort keys %peak2motif) {
	print OUT $peak,"\t",join(",",sort keys %{$peak2motif{$peak}}),"\n";
}
close OUT;


