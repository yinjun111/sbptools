#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);

#to be implemented !!!!


########
#Prerequisites
########



########
#Interface
########


my $version="0.1";

my $usage="

gtf-build
version: $version
Usage: sbptools gtf-build [parameters]

Description: Using an Ensembl gtf file to produce necessary files for NGS analysis, including:

	1) Convert Ensembl format to UCSC format
	2) Get composite_exons, 1kb_promoter_longesttx, 1kb_promoter_alltxs
	3) Generate gene/tx annotation to be used in RNA-Seq based on GTF
	4) Build library for STAR and RSEM


Parameters:

    --in|-i           Input folder(s)
	


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

my $samples;
my $inputfolders;

my $configfile;
my $outputfolder;
my $verbose;
my $tx;
my $runmode="none";

GetOptions(
	"in|i=s" => \$inputfolders,
	"samples|s=s" => \$samples,	
	"config|c=s" => \$configfile,
	"output|o=s" => \$outputfolder,
	"tx|t=s" => \$tx,	
	"runmode|r=s" => \$runmode,		
	"verbose|v" => \$verbose,
);

