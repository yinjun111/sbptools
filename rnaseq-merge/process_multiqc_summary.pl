#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);

########
#Interface
########


my $version="0.1";


my $usage="

process_multiqc_summary
version: $version
Usage: perl process_multiqc_summary [parameters]

Description: Process multiqc mltiqc_general_stats.txt


Parameters:

    --in|-i           Input file, multiqc_general_stats.txt
    --output|-o       Output file, multiqc_general_stats_rev.txt

    --config|-c       Configuration file
	
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
my $configfile;

my $verbose=1;
my $runmode="local";

my $dev=0; #developmental version

GetOptions(
	"in|i=s" => \$infile,
	"output|o=s" => \$outfile,
	"config|c=s" => \$configfile,

	"runmode|r=s" => \$runmode,		
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


my $mergefiles="$sbptoolsfolder/mergefiles/mergefiles_caller.pl";
my $text2excel="$sbptoolsfolder/text2excel/text2excel.pl";


########
#Code begins
########

my %configattrs;
#two output file, sample / R1+R2

my %samples;
my %fastq2sample;

#read config file to find fastq name
open(IN,$configfile) || die $!;
my $fileline=0;

while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	if($fileline==0) {
		#record title into hash
		for(my $num=0;$num<@array;$num++) {
			#case sensitive match to key words, Sample, FASTQ1, FASTQ2, Index (for foldernames)
			$configattrs{uc $array[$num]}=$num;
		}
		
		#should always use first column as sample name and index to maintain continuous workflow
		
		#first column as sample
		
		#check title
		unless(defined $configattrs{"FASTQ1"}) {
			#SAMPLE and FASTQ1 have to be defined. The others are optional
			print STDERR "ERROR: FASTQ1 need to be defined in $configfile. Input config includes:",join(",",map {uc $_} @array),"\n";
			exit;
		}
		else {
			print STDERR "Input config $configfile includes:",join(",",map {uc $_} @array),"\n\n" if $verbose;
		}
			
	}
	else {
		
		$samples{$array[$configattrs{"SAMPLE"}]}++;
		
		if(defined $configattrs{"FASTQ1"}) {
			my $fastq1name=basename($array[$configattrs{"FASTQ1"}]);
			
			if($fastq1name=~/(.+).fastq.gz/) {
				$fastq2sample{$1}=$array[$configattrs{"SAMPLE"}];
			}
		}


		if(defined $configattrs{"FASTQ2"}) {
			my $fastq2name=basename($array[$configattrs{"FASTQ2"}]);
			
			if($fastq2name=~/(.+).fastq.gz/) {
				$fastq2sample{$1}=$array[$configattrs{"SAMPLE"}];
			}
		}
	}
	
	$fileline++;
}

close IN;

my @tags=("TagDir","cutadapt");


#if fastq and sample can't merge, produce two files. If fastq and samples can merge, one should be enough

my $linenum=0;
my $title;

my %fastq2info;
my %sample2info;

open(IN,$infile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	
	if($linenum==0) {
		$title=$_;
	}
	else {
		if(defined $fastq2sample{$array[0]}) {
			$fastq2info{$array[0]}=[$fastq2sample{$array[0]},@array];
		}
		
		#rename sample name, remove tag
		my $samplename=$array[0];
		
		foreach my $tag (@tags) {
			if($array[0]=~/(.+)_$tag/) {
				$samplename=$1;
			}
		}
		
		#check sample info
		
		if(defined $samples{$samplename}) {
			unless(defined $sample2info{$samplename}) {
				$sample2info{$samplename}=[@array];
			}
			else {
				$sample2info{$samplename}=merge_info($sample2info{$samplename},[@array]);
			}
		}
	}
	$linenum++;
}
close IN;

my $outfile1=$outfile;
$outfile1=~s/.txt/_fastq.txt/;

my $outfile2=$outfile;	
$outfile2=~s/.txt/_sample.txt/;

open(OUT,">$outfile1") || die $!;
print OUT "Sample\t$title\n";
foreach my $fastq (sort keys %fastq2info) {
	print OUT join("\t",@{$fastq2info{$fastq}}),"\n";
}
close OUT;


open(OUT,">$outfile2") || die $!;
print OUT "$title\n";
foreach my $sample (sort keys %sample2info) {
	print OUT join("\t",@{$sample2info{$sample}}),"\n";
}
close OUT;	




#####
#Function
#####

sub merge_info {
	my ($array1,$array2)=@_;
	
	my $max_col=scalar(split("\t",$title));
	
	my @combined;
	
	for(my $num=0;$num<$max_col;$num++) {
		if(defined $array1->[$num]) {
			push @combined,$array1->[$num];
		}
		elsif(defined $array2->[$num]) {
			push @combined,$array2->[$num];
		}
		else {
			push @combined," ";
		}
	}
	
	return [@combined];
}

	
	


sub get_parent_folder {
	my $dir=shift @_;
	
	if($dir=~/^(.+\/)[^\/]+\/?/) {
		return $1;
	}
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


	
	


