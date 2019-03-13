#!/usr/bin/perl -w
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

#v0.1, now for RNA-Seq using Cybase


my $usage="

prepare_trackhub
version: $version
Usage: sbptools 

Description: Prepare trackhub files for selected samples

Parameters:

    --in|-i           input file
    --out|-o          out folder for trackhub
    --config|-c       sample config

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

my $infile;
my $outfolder;
my $type;
my $verbose;

GetOptions(
	"in|i=s" => \$infile,
	"out|o=s" => \$outfolder,
	"config|c=s" => \$config,

	"verbose|v" => \$verbose,
);


my $logfile=$outfile;
$logfile=~s/\.\w+$/_prepare_trackhub.log/;


#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";
print LOG "\n";



########
#Process
########






