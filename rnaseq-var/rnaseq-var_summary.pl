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

rnaseq-var_summary
version: $version
Usage: 

Description: 

Parameters:
    --in|-i           Input folder, containing rnaseq-var results
    --out|-o          Output file
    --cate|-c         Category, variants_impact_HIGH (default), variants_impact_LOW, variants_impact_MODERATE or variants_impact_MODIFIER

	
	
";

#In Falco, BSR Linux server, use screen+parallel to parallelly running jobs in background.


unless (@ARGV) {
	print STDERR $usage;
	exit;
}

my $params=join(" ",@ARGV);
#then call different scripts


########
#Parameters
########

my $infolder;
my $outfile;
my $cate="variants_impact_HIGH";

my $runmode=0;
my $verbose=1;


GetOptions(
	"in|i=s" => \$infolder,
	"out|o=s" => \$outfile,
	"cate|c=s" => \$cate,

	"verbose" => \$verbose,
	"help|h" => sub {print STDERR $usage;exit;}
	
);


my $logfile="rnaseq-var_summary_run.log";


######
#Process input file
######

#
my %cate2col=qw(variants_impact_HIGH 4 variants_impact_LOW 5 variants_impact_MODERATE 6 variants_impact_MODIFIER 7);


#get all summary.genes.txt
my @infiles=glob("$infolder/*/snpanno/*snpEff_summary.genes.txt");

my %sample2snp;
my %samples;
my %genes;

foreach my $file (@infiles) {
	
	print STDERR "$file\n";
	
	my $filename=basename($file);
	my $samplename;
	
	if($filename=~/(.+)snpEff_summary.genes.txt/) {
		$samplename=$1;
		$samples{$samplename}++;
	}
	
	open(IN,$file) || die "ERROR:Can't read $file.$!\n";
	while(<IN>) {
		tr/\r\n//d;
		next if $_=~/^#/;
		
		my @array=split/\t/;
		
		#variants_impact_HIGH
		$sample2snp{$samplename}{$array[0]}=$array[$cate2col{$cate}];
		$genes{$array[0]}++;
	}
	close IN;
}


open(OUT,">$outfile") || die $!;
print OUT "Gene\t",join("\t",sort keys %samples),"\n";

foreach my $gene (sort keys %genes) {
	my @marks;
	
	foreach my $sample (sort keys %samples) {
		if(defined $sample2snp{$sample}) {
			if(defined $sample2snp{$sample}{$gene}) {
				push @marks,$sample2snp{$sample}{$gene};
			}
			else {
				push @marks,0;
			}
		}
		else {
			push @marks,0;
		}
	}
	
	print OUT $gene,"\t",join("\t",@marks),"\n";
}
close OUT;

	

	
	






















