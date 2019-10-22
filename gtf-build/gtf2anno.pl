#!/usr/bin/perl -w
use strict;

#Ensembl GTF, txt anno
my ($infile,$outfile) = @ARGV;

#attr used from gtf
my @attrs=qw(gene_name gene_biotype transcript_id strand chromosome start end);

my %gene2info;

#read file
open(IN,$infile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^#/;
	my @array=split/\t/;
	if($array[2] eq "gene") {
		my %terms=process_info($array[8]);
		
		my $geneid=$terms{"gene_id"};
		
		foreach  my $attr (sort keys %terms) {
			$gene2info{$geneid}{$attr}{$terms{$attr}}++;
		}
		
		#str
		$gene2info{$geneid}{"strand"}=$array[6];
		#chr
		$gene2info{$geneid}{"chromosome"}=$array[0];
		#start, end
		$gene2info{$geneid}{"start"}=$array[3];
		$gene2info{$geneid}{"end"}=$array[4];
	}
	
	
	if($array[2] eq "transcript") {
		my %terms=process_info($array[8]);
		
		my $geneid=$terms{"gene_id"};
		
		foreach  my $attr (sort keys %terms) {
			$gene2info{$geneid}{$attr}{$terms{$attr}}++;
		}
	}
}


open(OUT,">$outfile") || die $!;
print OUT "Gene\t",join("\t",@attrs),"\n";
foreach my $geneid (sort keys %gene2info) {
	print OUT $geneid,"\t";
	my @contents;
	
	foreach my $attr (@attrs) {
		if(defined $gene2info{$geneid}{$attr}) {
			if(ref($gene2info{$geneid}{$attr}) eq "HASH") {
				push @contents,join(",",sort keys %{$gene2info{$geneid}{$attr}});
			}
			else {
				push @contents,$gene2info{$geneid}{$attr};
			}
		}
		else {
			push @contents," ";
		}
	}
	
	print OUT join("\t",@contents),"\n";
}

close OUT;
		

sub process_info {
	my $info=shift @_;
	
	my %infos;
	
	while($info=~/(\w+) ([^;]+);/g) {
		my $attr=$1;
		my $value=$2;
		$value=~tr/"//d;
		$infos{$attr}=$value;
	}
	
	return %infos;
}