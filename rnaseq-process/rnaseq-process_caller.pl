#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);

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

    --config|-c       Configuration file
    --output|-o       Output folder

    --tx|-t           Transcriptome
                        Current support Human.B38.Ensembl84, Mouse.B38.Ensembl84

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


my $configfile;
my $outputfolder;
my $verbose;
my $tx;
my $runmode="none";

GetOptions(
	"config|c=s" => \$configfile,
	"output|o=s" => \$outputfolder,
	"tx|t=s" => \$tx,	
	"runmode|r=s" => \$runmode,		
	"verbose|v" => \$verbose,
);

$outputfolder = abs_path($outputfolder);

if(!-e $outputfolder) {
	mkdir($outputfolder);
}


my $scriptfolder="$outputfolder/scripts";

if(!-e $scriptfolder) {
	mkdir($scriptfolder);
}


my $logfile="$outputfolder/rnaseq-process_run.log";
my $newconfigfile="$outputfolder/rnaseq-process_config.txt";

my $scriptfile1="$scriptfolder/rnaseq-process_run.sh";

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
print LOG "\n";


#test tx option

my %tx2ref=(
	"Human.B38.Ensembl84"=>"/home/jyin/Projects/Databases/Ensembl/v84/Human_STAR/Human_RSEM",
	"Mouse.B38.Ensembl84"=>"/home/jyin/Projects/Databases/Ensembl/v84/Mouse_STAR/Mouse_RSEM",
);


if(defined $tx2ref{$tx}) {
	print STDERR "Starting analysis using $tx.\n\n" if $verbose;
	print LOG "Starting analysis using $tx.\n\n";
}
else {
	print STDERR "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n" if $verbose;
	print LOG "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n";
}


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
open(IN,$configfile) || die "Error reading $configfile. $!";
open(OUT,">$newconfigfile") || die "Error reading $newconfigfile. $!";

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


#----------------
#create folder
print STDERR scalar(keys %sample2fastq)," samples identified from $configfile, including:\nSAMPLE\tINDEX\n",join("\n",map {$_."\t".$sample2indexname{$_}} sort keys %sample2indexname),"\n\n" if $verbose;
print LOG scalar(keys %sample2fastq)," samples identified from $configfile, including:\nSAMPLE\tINDEX\n",join("\n",map {$_."\t".$sample2indexname{$_}} sort keys %sample2indexname),"\n\n" if $verbose;

print STDERR "Make folders for different samples.\n\n" if $verbose;
print LOG "Make folders for different samples.\n\n";

foreach my $sample (sort keys %sample2fastq) {
	my $samplefolder="$outputfolder/$sample";
	if(!-e	$samplefolder) {
		mkdir $samplefolder;
	}
}

#----------------
#prepare fastq file for cutadapt
#assume truseq for now




my %sample2workflow;

#FASTQ files should all be .fastq.gz

if(defined $configattrs{"FASTQ2"}) {
	#PE
	print STDERR "Printing cutadapt PE script.\n\n" if $verbose;
	print LOG "Printing cutadapt PE script.\n\n";
	
	foreach my $sample (sort keys %sample2fastq) {
		my $samplefolder="$outputfolder/$sample";
		
		my $fastq1=$sample2fastq{$sample}[0];
		my $fastq2=$sample2fastq{$sample}[0];
		
		my $fastq1trim="$samplefolder/".basename($fastq1);
		$fastq1trim=~s/\.fastq\.gz/_trimmed.fastq.gz/;
		push @{$sample2fastq{$sample}},$fastq1trim;

		my $fastq2trim="$samplefolder/".basename($fastq2);
		$fastq2trim=~s/\.fastq\.gz/_trimmed.fastq.gz/;
		push @{$sample2fastq{$sample}},$fastq2trim;

		$sample2workflow{$sample}.="$cutadapt -j 4 -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o $fastq1trim -p $fastq2trim $fastq1 $fastq2;";
	}
}
else {
	#SE
	print STDERR "Printing cutadapt SE script.\n\n" if $verbose;
	print LOG "Printing cutadapt SE script.\n\n";
	
	foreach my $sample (sort keys %sample2fastq) {
		my $samplefolder="$outputfolder/$sample";
	
		my $fastq1=$sample2fastq{$sample}[0];
		
		my $fastq1trim="$samplefolder/".basename($fastq1);
		$fastq1trim=~s/\.fastq\.gz/_trimmed.fastq.gz/;
		push @{$sample2fastq{$sample}},$fastq1trim;

		$sample2workflow{$sample}.="$cutadapt -j 4 -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o $fastq1trim $fastq1;";
	}
}
		

#----------------
#FASTQC for trimmed reads


if(defined $configattrs{"FASTQ2"}) {
	#PE
	print STDERR "Printing FASTQC PE script.\n\n" if $verbose;
	print LOG "Printing FASTQC PE script.\n\n";
	
	foreach my $sample (sort keys %sample2fastq) {
		my $samplefolder="$outputfolder/$sample";
		
		$sample2workflow{$sample}.="$fastqc --nogroup -o $samplefolder -f fastq ".$sample2fastq{$sample}[2]." ".$sample2fastq{$sample}[3].";";
	}
}
else {
	#SE
	print STDERR "Printing FASTQC SE script.\n\n" if $verbose;
	print LOG "Printing FASTQC SE script.\n\n";
	
	foreach my $sample (sort keys %sample2fastq) {
		my $samplefolder="$outputfolder/$sample";
		
		$sample2workflow{$sample}.="$fastqc --nogroup -o $samplefolder -f fastq ".$sample2fastq{$sample}[1].";";
	}
}	
		
	
#----------------
#RSEM for trimmed reads

if(defined $configattrs{"FASTQ2"}) {
	#PE
	print STDERR "Printing RSEM PE script.\n\n" if $verbose;
	print LOG "Printing RSEM PE script.\n\n";
	
	foreach my $sample (sort keys %sample2fastq) {
		my $samplefolder="$outputfolder/$sample";
		
		$sample2workflow{$sample}.="$rsem -p 4 --output-genome-bam --sort-bam-by-coordinate --star-gzipped-read-file --star --paired-end ".$sample2fastq{$sample}[2]." ".$sample2fastq{$sample}[3]." ".$tx2ref{$tx}." $samplefolder/$sample;";
	}
}
else {
	#SE
	print STDERR "Printing RSEM SE script.\n\n" if $verbose;
	print LOG "Printing RSEM SE script.\n\n";
	
	foreach my $sample (sort keys %sample2fastq) {
		my $samplefolder="$outputfolder/$sample";
		
		$sample2workflow{$sample}.="$rsem -p 4 --output-genome-bam --sort-bam-by-coordinate --star-gzipped-read-file --star ".$sample2fastq{$sample}[1]." ".$tx2ref{$tx}." $samplefolder/$sample;";
	}
}



########
#Print out commands, for local and server run
########


open(S1,">$scriptfile1") || die "Error writing $scriptfile1. $!";

foreach my $sample (sort keys %sample2workflow) {
	print S1 $sample2workflow{$sample},"\n";
}

close S1;


#local mode
print STDERR "To run locally, in shell type: sh $scriptfile1\n\n";
print LOG "To run locally, in shell type: sh $scriptfile1\n\n";



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












