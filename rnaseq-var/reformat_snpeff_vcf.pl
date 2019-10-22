#!/usr/bin/perl -w
use strict;

my ($infile,$outfile) = @ARGV;


my $header="F";


#infile, snpeff vcf
#outfile, split by ann

open(IN,$infile) || die $!;
open(OUT,">$outfile") || die $!;

print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\t";
#snpeff annotations
print OUT "Allele	Annotation	Putative_impact	Gene Name	Gene ID	Feature type	Feature ID	Transcript biotype	Rank / total	HGVS.c	HGVS.p	cDNA_position / (cDNA_len optional)	CDS_position / (CDS_len optional):	Protein_position / (Protein_len optional)	Distance to feature	Errors, Warnings or Information messages\n";

while(<IN>) {
	tr/\r\n//d;
	if ($_=~/^\#/) {
		if($header eq "T") {
			print OUT $_,"\n";
		}
	}
	else {	
		my @array=split/\t/;
		if($array[7]=~/ANN=([^;]+)/) {
			my $annos=$1;
			my $info=$`.$';
			
			foreach my $anno (split(",",$annos)) {
				print OUT join("\t",@array[0..6],$info,@array[8,9], split("\\|",$anno)),"\n";
			}
		}
		else {
			print OUT $_,"\t";
			print OUT join("\t",(" ") x 16),"\n";
		}
	}
}
close IN;
close OUT;

