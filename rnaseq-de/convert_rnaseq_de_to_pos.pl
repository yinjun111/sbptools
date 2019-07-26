#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);


########
#Prerequisites
########


#my $mergefiles="perl /home/jyin/Projects/Pipeline/sbptools/mergefiles/mergefiles_caller.pl";

########
#Interface
########


my $version="0.1";

my $usage="

convert_rnaseq_de_to_pos
version: $version
Usage: perl  [parameters]

Description: convert rnaseq de peaks to pos for tfbs analysis


Mandatory Parameters:
    --in|-i           Input file from chipseq-de
    --out|--o         Output file, pos files for up,down,both
		
    --fccutoff        Log2 FC cutoff [1]
    --qcutoff         Corrected P cutoff [0.05]

    --runmode|-r      Where to run the scripts, local, server or none [none]
    --verbose|-v      Verbose
	
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

my ($infile,$outfile);

my $fccutoff;
my $qcutoff;
	

my $verbose=1;
my $runmode="none";

GetOptions(
	"in|i=s" => \$infile,
	"out|o=s" => \$outfile,
	
	"fccutoff=s" => \$fccutoff,
	"qcutoff=s" => \$qcutoff,

	"runmode|r=s" => \$runmode,		
	"verbose|v" => \$verbose,
);


#convert chipseq de file to pos file for following TFBS analysis

#infile is ChIP-Seq DE
#outfile is DE pos

#pos format
#1. peak name (should be unique)
#2. chromsome
#3. starting position [integer] (1-indexed)
#4. end position [integer]
#5. strand [either 0/1 or +/-] (in HOMER strand of 0 is +, 1 is -)

my $upoutfile=$outfile;
$upoutfile=~s/(\.\w+)$/_up$1/;

my $downoutfile=$outfile;
$downoutfile=~s/(\.\w+)$/_down$1/;



open(IN,$infile) || die $!;
open(OUT,">$outfile") || die $!;
open(OUT1,">$upoutfile") || die $!;
open(OUT2,">$downoutfile") || die $!;

while(<IN>) {
	tr/\r\n//d;
	next unless $_=~/^\w/;
	my @array=split/\t/;
	#array[0] should be real coordinate
	
	if($array[0]=~/([^:]+)\:(\d+)-(\d+)/) {
		unless(defined $fccutoff && length($fccutoff)>0 && defined $qcutoff && length($qcutoff)>0) {
			if($array[5] eq "1" || $array[5] eq "-1") {	
				if($array[5] eq "1") {
					print OUT1 $array[0],"\t",$1,"\t",$2,"\t",$3,"\t",".","\n";
				}
				else {
					print OUT2 $array[0],"\t",$1,"\t",$2,"\t",$3,"\t",".","\n";
				}
				print OUT $array[0],"\t",$1,"\t",$2,"\t",$3,"\t",".","\n";
			}
		}
		else {
			if($array[1] >= $fccutoff && $array[4] ne "NA" && $array[4] <$qcutoff) {
				print OUT1 $array[0],"\t",$1,"\t",$2,"\t",$3,"\t",".","\n";
				print OUT $array[0],"\t",$1,"\t",$2,"\t",$3,"\t",".","\n";
			}
			elsif($array[1] <= $fccutoff && $array[4] ne "NA" && $array[4] <$qcutoff) {
				print OUT2 $array[0],"\t",$1,"\t",$2,"\t",$3,"\t",".","\n";
				print OUT $array[0],"\t",$1,"\t",$2,"\t",$3,"\t",".","\n";	
			}
		}
	}
}
close IN;
close OUT;
close OUT1;
close OUT2;

	
	