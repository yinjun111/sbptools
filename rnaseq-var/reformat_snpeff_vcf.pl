#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);


########
#Interface
########


my $version="0.2";

#v0.2, add summary for TMB calculation

my $usage="

reformat_snpeff_vcf
version: $version
Usage: perl reformat_snpeff_vcf.pl -i infile -o outfile

Description: reformat vcf file to split different functional annotations into different rows


Parameters:
    --input|-i        Input file
    --output|-o       Output file
								   
";

#--verbose|-v      Verbose

unless (@ARGV) {
	print STDERR $usage;
	exit;
}

#my $params=join(" ",@ARGV);



########
#Parameters
########

my $infile;
my $outfile;
my $dev=0; #developmental version

my $verbose=1;
my $header="F";


GetOptions(
	"input|i=s" => \$infile,
	"output|o=s" => \$outfile,
	"verbose|v" => \$verbose,
	"dev" => \$dev,		
);


########
#Prerequisites
########



my $sbptoolsfolder="/apps/sbptools/";

#adding --dev switch for better development process
if($dev) {
	$sbptoolsfolder="/home/jyin/Projects/Pipeline/sbptools/";
}
else {
	#the tools called will be within the same folder of the script
	$sbptoolsfolder=get_parent_folder(abs_path(dirname($0)));
}


my $tmbconfig="$sbptoolsfolder/rnaseq-var/Snpeff_anno_rev_V2.txt";



########
#Program starts
########

my $summaryfile=$outfile;

$summaryfile=~s/\.txt$/\_summary.txt/; #summary for TMBs


#read TMB config
my $linenum=0;
my %effect2tmb;
my @tmbs;
my @effects;
my %effects;

open(IN,$tmbconfig) || die "ERROR:$tmbconfig can't be read.$!\n";
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	if($linenum==0) {
		@tmbs=@array[4,5,6];
	}
	else {
		for(my $num=4;$num<=6;$num++) {
			if($array[$num] eq "Y") {
				$effect2tmb{$array[0]}{$tmbs[$num-4]}++;
			}
		}
		
		push @effects,$array[0];
		
	}
	$linenum++;
}
close IN;

%effects=map {$_,1} @effects;


#infile, snpeff vcf
#outfile, split by ann

my %tmb2vars;
my %effect2vars;

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
				my @anno_results=split("\\|",$anno);
				print OUT join("\t",@array[0..6],$info,@array[8,9], @anno_results),"\n";
				
				foreach my $effect (split("&",$anno_results[1])) {
					if(defined $effect2tmb{$effect}) {
						foreach my $tmb (sort keys %{$effect2tmb{$effect}}) {
							$tmb2vars{$tmb}{join("-",@array[0,1,3,4])}++; #chr,pos,ref,alt
						}
					}
					else {
						unless(defined $effects{$effect}) {					
							print STDERR "ERROR:$effect not defined.\n";
						}
					}
					
					$effect2vars{$effect}{join("-",@array[0,1,3,4])}++;
				}
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

#produce a report for missense mutations, all exonic mutations, all mutations 
open(OUT,">$summaryfile") || die $!;

foreach my $tmb (@tmbs) {
	if(defined $tmb2vars{$tmb}) {
		print OUT $tmb,"\t",scalar(keys %{$tmb2vars{$tmb}}),"\n";
	}
	else {
		print OUT $tmb,"\t0\n";
	}
}
print OUT "\n";

foreach my $effect (@effects) {
	if(defined $effect2vars{$effect}) {
		print OUT $effect,"\t",scalar(keys %{$effect2vars{$effect}}),"\n";
	}
	else {
		print OUT $effect,"\t0\n";
	}
}

close OUT;



#######
#Functions
#######


sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}

sub getsysoutput {
	my $command=shift @_;
	my $output=`$command`;
	$output=~tr/\r\n//d;
	return $output;
}


sub build_timestamp {
	my ($now,$opt)=@_;
	
	if($opt eq "long") {
		$now=~tr/ /_/;
		$now=~tr/://d;
	}
	else {
		$now=substr($now,0,10);
	}
	
	return $now;
}


sub find_program {
	my $fullprogram=shift @_;
	
	#use defined program as default, otherwise search for this program in PATH
	
	my $program;
	if($fullprogram=~/([^\/]+)$/) {
		$program=$1;
	}
	
	if(-e $fullprogram) {
		return $fullprogram;
	}
	else {
		my $sysout=`$program`;
		if($sysout) {
			my $location=`which $program`;
			return $location;
		}
		else {
			print STDERR "ERROR:$fullprogram or $program not found in your system.\n\n";
			exit;
		}
	}
}


sub get_parent_folder {
	my $dir=shift @_;
	
	if($dir=~/^(.+\/)[^\/]+\/?/) {
		return $1;
	}
}



