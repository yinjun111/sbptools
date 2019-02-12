#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);

#CutAdapt+FASTQC+RSEM+STAR


########
#Prerequisites
########

my $cutadapt="/home/jyin/.local/bin/cutadapt";
my $fastqc="/apps/FastQC/fastqc";
my $rsem="/home/jyin/Programs/RSEM-1.3.1/rsem-calculate-expression";
my $star="/home/jyin/Programs/STAR-2.6.1c/bin/Linux_x86_64/STAR";



########
#Interface
########


my $version="0.1";

my $usage="

rnaseq-process
version: $version
Usage: sbptools rnaseq-process [parameters]


Parameters:

    -c/--config       Configuration file
    -o/--output       Output folder


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


my $configfile;
my $outputfolder;
my $verbose;

GetOptions(
	"config|c=s" => \$configfile,
	"output|o=s" => \$outputfolder,
	"verbose|v" => \$verbose,
);

$outputfolder = abs_path($outputfolder);

if(!-e $outputfolder) {
	mkdir($outputfolder);
}


my $scriptfolder="$outputfolder/scripts";
my $logfile="$outputfolder/rnaseq-process_run.log";
my $newconfigfile="$outputfolder/rnaseq-process_config.txt";


#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";
print LOG "$cutadapt version:", getsysoutput("$cutadapt --version"),"\n";
print LOG "$fastqc version:", getsysoutput("$fastqc --version"),"\n";
print LOG "$rsem version:", getsysoutput("$rsem --version"),"\n";
print LOG "$star version:", getsysoutput("$star --version"),"\n";



########
#Process
########

print STDERR "\nsbptools rnaseq-process running ...\n\n" if $verbose;
print LOG "\nsbptools rnaseq-process running ...\n\n";

#open config files to find fastqs

my %sample2fastq;
my %sample2indexname;
my %configattrs;

my $fileline=0;
my $newconfigfiletitle;
open(IN,$configfile) || die "Error opening $configfile. $!";
open(OUT,">$newconfigfile") || die "Error opening $newconfigfile. $!";

while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	if($fileline==0) {
		for(my $num=0;$num<@array;$num++) {
			#case sensitive match to key words, Sample, FASTQ1, FASTQ2, Index (for foldernames)
			$configattrs{uc $array[$num]}=$num;
		}
		
		#check title
		unless(defined $configattrs{"SAMPLE"} && defined $configattrs{"FASTQ1"}) {
			#SAMPLE and FASTQ1 have to be defined. The others are optional
			print STDERR "ERROR: SAMPLE and FASTQ1 need to be defined in $configfile. Input config includes:",join(",",map {uc $_} @array),"\n";
			print LOG "ERROR: SAMPLE and FASTQ1 need to be defined in $configfile. Input config includes:",join(",",map {uc $_} @array),"\n";
			exit;
		}
		else {
			print STDERR "Input config $configfile includes:",join(",",map {uc $_} @array),"\n\n" if $verbose;
			print LOG "Input config $configfile includes:",join(",",map {uc $_} @array),"\n\n";
		}
		
		if(defined $configattrs{"INDEX"}) {
			$newconfigfiletitle=join("\t","INDEX",map {uc $_} (splice @array,$configattrs{"INDEX"},1));
		}
		else {
			$newconfigfiletitle=join("\t","INDEX",map {uc $_} @array);
		}
		
		print OUT $newconfigfiletitle,"\n";
		
		print STDERR "New config $newconfigfile includes:",join(",",split("\t",$newconfigfiletitle)),"\n\n" if $verbose;
		print LOG "New config $newconfigfile includes:",join(",",split("\t",$newconfigfiletitle)),"\n\n";		
	
	}
	else {
		
		my $indexname;
		my $newconfigline;
		
		#need to change non word chars to words
		if(defined $configattrs{"INDEX"}) {
			$indexname=$array[$configattrs{"INDEX"}];
		}
		else {
			$indexname=$array[$configattrs{"SAMPLE"}];
			$indexname=~s/[^\w-]/_/g;
		}
		
		$sample2indexname{$array[$configattrs{"SAMPLE"}]}=$indexname;
		
		#fastq1
		push @{$sample2fastq{$indexname}},$array[$configattrs{"FASTQ1"}];
		
		#fastq2 (Read2)
		if(defined $configattrs{"FASTQ2"}) {
			push @{$sample2fastq{$indexname}},$array[$configattrs{"FASTQ2"}];
		}
		
		if(defined $configattrs{"INDEX"}) {
			$newconfigline=join("\t",$indexname,(splice @array,$configattrs{"INDEX"},1));
		}
		else {
			$newconfigline=join("\t",$indexname,@array);
		}
		
		print OUT $newconfigline,"\n";		
		
	}
	
	$fileline++;
}

close IN;
close OUT;

			


#prepare fastq file for cutadapt
#assume truseq for now









close LOG;
########
#Functions
########

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












