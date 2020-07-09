#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);

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

rnaseq-var_filter
version: $version
Usage: 

Description: 

Parameters:
    --in|-i           Input file, a vcf file
    --out|-o          Output file, a filtered vcf file
    --com|-c          Common SNP vcf
    --dp|-d           DP filter
	
	
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
my $commonsnp;
my $dp;

my $runmode=0;
my $verbose=1;


GetOptions(
	"in|i=s" => \$infile,
	"out|o=s" => \$outfile,
	"com|c=s" => \$commonsnp,
	"dp|d=s" => \$dp,

	"verbose" => \$verbose,
	"help|h" => sub {print STDERR $usage;exit;}
	
);


my $logfile="rnaseq-var_filter_run.log";


######
#Process input file
######

print STDERR "\nrnaseq-var_filter version $version running...\n";

#Read common snp inforamtion
my %commonsnps;
if(defined $commonsnp && length($commonsnp)>0) {
	open(IN,$commonsnp) || die $!;
	while(<IN>) {
		tr/\r\n//d;
		next if ($_=~/^\#/);

		my @array=split/\t/;
		
		my $snpname=join(":",@array[0,1,3,4]); #index by chr:coord:ref:alt

		$commonsnps{$snpname}++;
	}
	close IN;
}


open(IN,$infile) || die $!;
open(OUT,">$outfile") || die $!;

my $anum=0;#all snps
my $cnum=0;#common snps
my $dnum=0;#dp snps

while(<IN>) {
	tr/\r\n//d;
	if ($_=~/^\#/) {
		print OUT $_,"\n";
	}
	else {	
		$anum++;
		
		my @array=split/\t/;
		
		my $snpname=join(":",@array[0,1,3,4]); #index by chr:coord:ref:alt
		
		#filter common snp
		if(defined $commonsnps{$snpname}) {
			$cnum++;
			next;
		}
		
		#filter existing DP
		if( (!defined $dp || length($dp)==0) && $array[6] =~/DP/) {
			$dnum++;
			next;
		}
		
		#use new dp filter
		if(defined $dp && length($dp)>0) {
			if($array[7]=~/DP=(\d+)/) {
				if($1<=$dp) {
					$dnum++;
					next;
				}
			}
		}
		
		print OUT $_,"\n";
		
	}
}
close IN;
close OUT;

print STDERR "$anum SNPs are processed using rnaseq-var_filter.\n$cnum SNPs are filtered by matching common SNPs.\n$dnum SNPs are filtered by matching DP filter.\n\n";