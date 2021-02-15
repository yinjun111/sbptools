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

#v0.1, summarize 



my $usage="

dnaseq-summary
version: $version
Usage: Summarize TMBs and variants for different TMB categories

Description: 

Parameters:
    --in|-i           Input folder(s), containing dnaseq results
    --out|-o          Output folder

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

my $inputfolder;
my $outputfolder;

my $runmode=0;
my $verbose=1;


GetOptions(
	"in|i=s" => \$inputfolder,
	"out|o=s" => \$outputfolder,

	"verbose" => \$verbose,
	"help|h" => sub {print STDERR $usage;exit;}
	
);


########
#Require file
########


my $snpeffanno="/apps/sbptools/rnaseq-var/Snpeff_anno_rev_V2.txt";



########
#Parameters
########

if(!-e $outputfolder) {
	mkdir($outputfolder);
}

$outputfolder = abs_path($outputfolder);

my $logfile="$outputfolder/dnaseq_summary_run.log";

#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();
my $timestamp=build_timestamp($now,"long");

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";



#output files

my $tmbsummary="$outputfolder/dnaseq_tmb_summary.txt";

#genotype
#to be implemented


######
#Process input file
######


print STDERR "\nsbptools dnaseq-summary $version running ...\n\n" if $verbose;
print LOG "\nsbptools dnaseq-summary $version running ...\n\n";


#Write TMB summary file
#--------------------
print STDERR "Summarizing TMB counts:\n";
print LOG "Summarizing TMB counts:\n";

my %samplefolders; #remember the folders for next step
my %samplenames;
my %sample2tmb;

foreach my $folder (split(",",$inputfolder)) {
	my @subfolders=glob("$folder/*"); #get all folders inside each infolder
	
	foreach my $samplefolder (@subfolders) {
		
		if(-e "$samplefolder/snpanno/") {
			#this is a var calling folder
			
			my $samplename=basename($samplefolder);
			
			if(defined $samplenames{$samplename}) {
				$samplenames{$samplename}++;
				$samplename=$samplename."_".$samplenames{$samplename};
				$samplefolders{$samplename}=abs_path($samplefolder);
			}
			else {
				$samplenames{$samplename}++;
				$samplefolders{$samplename}=abs_path($samplefolder);
			}
			
			#for next dnaseq-process, it needs to remove sample name from ouptutfiles
			my $varsummaryfile=(glob("$samplefolder/snpanno/*.filtered-cleaned.snpeff.sift.annotated.edited_summary.txt"))[0];
			
			if(-e $varsummaryfile) {
				print STDERR "Processing $samplename\n";
				print LOG "Processing $samplename: $varsummaryfile\n";
				
				open(IN,"$varsummaryfile") || die $!;
				while(<IN>) {
					tr/\r\n//d;
					my @array=split/\t/;
					if($_=~/^TMB1/) {
						$sample2tmb{$samplename}{"TMB1"}=$array[1];
					}
					if($_=~/^TMB2/) {
						$sample2tmb{$samplename}{"TMB2"}=$array[1];
					}				
					if($_=~/^TMB3/) {
						$sample2tmb{$samplename}{"TMB3"}=$array[1];
					}
				}
				close IN;
			}
			else {
				print STDERR "ERROR:$samplename in $samplefolder doesn't have summary file $varsummaryfile.\n";
			}
		}
	}
}


#write TMB summary
open(OUT,">$tmbsummary") || die $!;
print OUT "Sample\tTMB1(Foundation)\tTMB2(Missense)\tTMB3(Protein Changing)\n";

foreach my $sample (sort keys %sample2tmb) {
	print OUT $sample,"\t";
	print OUT $sample2tmb{$sample}{"TMB1"},"\t";
	print OUT $sample2tmb{$sample}{"TMB2"},"\t";
	print OUT $sample2tmb{$sample}{"TMB3"},"\n";
}

close OUT;



close LOG;


####
#Functions
####



sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
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






