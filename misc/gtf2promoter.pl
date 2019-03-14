#!/usr/bin/perl -w
use strict;

#use GTF file to output promoter bed file

#File1, all txs
#File2, longest tx per gene

use strict;
use Getopt::Long;

########
#Prerequisites
########

#none


########
#Interface
########


my $version="0.1";


my $usage="

gtf2promoter
version: $version
Usage: perl gtf2promoter.pl --gtf human.gtf -u 5000 -d 100 -o human_promoter_5kup100down

Description: extract promoter sequences for all transcripts and longest transcripts per gene from gtf

Parameters:

    --gtf             GTF file
    --up              Upstream length [2000]
    --down            Downstream length [0]
    --out|-o          Output file suffix
                        out_anno.txt
                        out_alltxs.bed
                        out_longesttxs.bed
	
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

my $gtf;
my $upstream=2000;
my $downstream=0;
my $size;
my $out;

my $verbose;

GetOptions(
	"gtf=s" => \$gtf,
	"up|u=s" => \$upstream,
	"down|d=s" => \$downstream,
	"size=s" => \$size,
	"out|o=s" => \$out,
	"verbose|v" => \$verbose,
);



########
#Process
########

my $txannofile="$out\_anno.txt";
my $txallbed="$out\_alltxs.bed";
my $txlongbed="$out\_longesttxs.bed";

my %gene2tx;
my %tx2info; #chr, start, end, str

open(IN,$gtf) || die $!;

while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^#/;
	my @array=split/\t/;
	
	if($array[2] eq "transcript") {
		my %terms=process_info($array[8]);
		
		$gene2tx{$terms{"gene_id"}}{$terms{"transcript_id"}}=$array[4]-$array[3]+1;
		$tx2info{$terms{"transcript_id"}}=[@array[0,3,4,6]];#chr,start,end,strand
	}
}
close IN;

#read whole chr size
my %chr2size;
open(IN,$size) || die $!;
while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^#/;
	my @array=split/\t/;
	$chr2size{$array[0]}=$array[1];
}
close IN;


#outputs
open(OUT1,">$txannofile") || die $!;
open(OUT2,">$txallbed") || die $!;
open(OUT3,">$txlongbed") || die $!;

print OUT1 "Tx\tGene\tChr\tStart\tEnd\tStr\tPromoterStart_Up${upstream}Down$downstream\tPromoterEnd_Up${upstream}Down$downstream\tLengthRank\n";
foreach my $gene (sort keys %gene2tx) {
	my @txs=sort {$gene2tx{$gene}{$b}<=>$gene2tx{$gene}{$a}} keys %{$gene2tx{$gene}};
	
	for(my $num=0;$num<@txs;$num++) {
		my $tx=$txs[$num];
		my $upcoord;
		my $downcoord;
		
		
		unless($tx2info{$tx}[3] eq "-") {
			#+ str
			$upcoord=$tx2info{$tx}[1]>$upstream?$tx2info{$tx}[1]-$upstream-1:0;  #bed coord, chr 5'
			$downcoord=$tx2info{$tx}[1]+$downstream-1; #chr 3'
		}
		else {
			#- str
			$upcoord=$tx2info{$tx}[2]>$downstream?$tx2info{$tx}[2]-$downstream:0;  #bed coord, chr 5'
			$downcoord=$tx2info{$tx}[2]+$upstream<$chr2size{$tx2info{$tx}[0]}?$tx2info{$tx}[2]+$upstream:$chr2size{$tx2info{$tx}[0]}; #chr 3'
		}
		
		if($num==0) {
			print OUT3 $tx2info{$tx}[0],"\t$upcoord\t$downcoord\t$tx\t.\t",$tx2info{$tx}[3],"\n";
		}
		
		print OUT2 $tx2info{$tx}[0],"\t$upcoord\t$downcoord\t$tx\t.\t",$tx2info{$tx}[3],"\n";
		
		print OUT1 $tx,"\t",$gene,"\t",join("\t",@{$tx2info{$tx}}),"\t",$upcoord,"\t",$downcoord,"\t",$num+1,"\n";
	}
			
}
close OUT1;
close OUT2;
close OUT3;



########
#Functions
########


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
