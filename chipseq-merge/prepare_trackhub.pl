#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);

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

    --in|-i           input folder
    --out|-o          out folder for trackhub
    --config|-c       sample config (optional)
    --tx|-t           Transcriptome
                        Current support Human.B38.Ensembl84, Mouse.B38.Ensembl84	
    --project|-p      Project name

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

my $infolder;
my $outfolder;
my $config;
my $type;
my $verbose;
my $tx;
my $project;


GetOptions(
	"in|i=s" => \$infolder,
	"out|o=s" => \$outfolder,
	"config|c=s" => \$config,
	"project|p=s" => \$project,	
	"tx=s" =>\$tx,
	"verbose|v" => \$verbose,
);





if(!-e $outfolder) {
	mkdir($outfolder);
}

$outfolder = abs_path($outfolder);



#my $logfile="prepare_trackhub.log";


#write log file
#open(LOG, ">$outfolder/$logfile") || die "Error writing into $logfile. $!";

#my $now=current_time();

#print LOG "perl $0 $params\n\n";
#print LOG "Start time: $now\n\n";
#print LOG "Current version: $version\n\n";
#print LOG "\n";


########
#Process
########

#read config

my %attr2key;

open(IN,$config) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	$attr2key{$array[0]}=$array[1];
}
close IN;



####
#make hub.txt
###
open(OUT,">$outfolder/hub.txt") || die $!;

print OUT <<EOF;
hub $project
shortLabel $project
longLabel $project
genomesFile $attr2key{"genomes.txt"}
email jyin\@sbpdiscovery.org
descriptionUrl
EOF
close OUT;

####
#make genome .txt
####

open(OUT,">$outfolder/genomes.txt") || die $!;

print OUT <<EOF;
genome mm10
trackDb $attr2key{"trackdb.txt"}
EOF
close OUT;


####
#support bb,bw now
####

opendir(DIR,$infolder) || die $!;
my %file2type;

my @files=readdir(DIR);

close DIR;


####
#tracks
####

my $genome;
if($tx=~/Human.B38/) {
	$genome="hg38";
}
elsif($tx=~/Mouse.B38/) {
	$genome="mm10";
}

unless(-e "$outfolder/$genome") {
	mkdir("$outfolder/$genome");
}

open(OUT,">$outfolder/$genome/trackdb.txt") || die $!;

foreach my $file (@files) {
	my $filename=basename($file);

	my $samplename=$filename;
	$samplename=~s/\.\w+$//;

	if($file=~/.bw$/) {
		print OUT <<EOF;
track $samplename-Signal
bigDataUrl $attr2key{"trackdb"}/$filename
shortLabel $samplename-Signal
longLabel $samplename-Signal
type bigWig
maxHeightPixels 50:50:0
viewLimitsMax 5

EOF
	}
	
	if($file=~/.bb$/) {
		print OUT <<EOF;
track $samplename-Peaks
bigDataUrl $attr2key{"trackdb"}/$filename
shortLabel $samplename-Peaks
longLabel $samplename-Peaks
type bigBed

EOF
	}
}

close OUT;


########
#Functions
########

sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}
