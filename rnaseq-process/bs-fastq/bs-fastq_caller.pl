#!/usr/bin/perl -w
use strict;
use File::Basename qw(basename);
use Getopt::Long;
use Cwd qw(abs_path);

#find and download fastq files from Basespace basemount directories

########
#Prerequisites
########

#basemount 
my $basemount="/usr/local/bin/basemount";


########
#Interface
########


my $version="0.2";

#v0.2 adding function to merge fastq files by R1, R2

my $usage="

bs-fastq
version: $version
Usage: sbptools bs-fastq -p projectname -f basemount_folder -o outputfolder

Description: Copy Fastq files from basemount project folder

Parameters:

    --folder|-f       Basemount folder name
    --project|-p      Project name
    --sample|-s       Sample name(s) (optional)
    --out|-o          output file

    --merge|-m        Whether to merge fastq files from each folder by R1/R2 [T]
    --config|-c       Config file, tab delimited file Folder,Project,Sample (optional)

	
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

my $folder;
my $project;
my $sample;
my $config;
my $outputfolder;
my $merge="T";
my $runmode="none";

my $verbose;

GetOptions(
	"folder|f=s" => \$folder,
	"project|p=s" => \$project,
	"sample|s=s" => \$sample,	
	"out|o=s" => \$outputfolder,
	"config|c=s" => \$config,
	"merge|m=s" => \$merge,
	"runmode|r" => \$runmode,
	
	"verbose|v" => \$verbose,
);



########
#Output folder & log
########

$outputfolder = abs_path($outputfolder);

if(!-e $outputfolder) {
	mkdir($outputfolder);
}

my $mergedfolder="$outputfolder/Merged";
if($merge eq "T") {
	mkdir($mergedfolder);
}


my $logfile="$outputfolder/bs-fastq_run.log";
my $scriptfile="$outputfolder/bs-fastq_run.sh";

#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";
print LOG "$basemount version:", getsysoutput("$basemount -V"),"\n";


########
#Progress
########

my %samplelist;


if(defined $config && length($config)>0) {
	#use config file to retrieve information
	my %attr2col;
	my $cline=0;
	open(IN,$config) || die "Error openning $config. $!";
	while(<IN>) {
		tr/\r\n//d;
		my @array=split/\t/;
		if($cline==0) {
			for(my $num=0;$num<@array;$num++) {
				$attr2col{uc $array[$num]}=$num;
			}
			
			unless(defined $attr2col{"PROJECT"} && defined $attr2col{"FOLDER"} && defined $attr2col{"SAMPLE"}) {
				print STDERR "ERROR:PROJECT, FOLDER, SAMPLE need to be defined in $config.\n";
				print LOG "ERROR:PROJECT, FOLDER, SAMPLE need to be defined in $config.\n";
				exit;
			}
			
		}
		else {
			if(-e $array[$attr2col{"FOLDER"}]) {
				my $infolder=abs_path($array[$attr2col{"FOLDER"}]);
				
				my $samplefolder="$infolder/Projects/".$array[$attr2col{"PROJECT"}]."/Samples/".$array[$attr2col{"SAMPLE"}];
				
				if(-e $samplefolder) {
					$samplelist{$infolder}{$array[$attr2col{"PROJECT"}]}{$array[$attr2col{"SAMPLE"}]}++;
				}
				else {
					print STDERR "ERROR:$samplefolder doesn't exist.\n";
					print LOG "ERROR:$samplefolder doesn't exist.\n";
					exit;
				}
			}
			else {
				print STDERR "ERROR:",$array[$attr2col{"FOLDER"}]," doesn't exist in $_\n";
				print LOG "ERROR:",$array[$attr2col{"FOLDER"}]," doesn't exist in $_\n";
				exit;
			}
		}
		$cline++;
	}
	close IN;
	
}
else {
	
	if(-e $folder) {
		my $infolder=abs_path($folder);
		
		#search for samples
		if(defined $sample && length($sample)>0) {
			foreach my $s (split(",",$sample)) {
				$samplelist{$infolder}{$project}{$s}++;
			}
		}
		else {
			#~/Data/Basespace/Projects/XU_XZ_ATAC_1_14_19/Samples
			
			my $inputfolder="$folder/Projects/$project/Samples/";
			
			opendir my($dh), $inputfolder || die "ERROR:Couldn't open dir $inputfolder. $!";
			my @files = readdir $dh;
			closedir $dh;
			
			#find samples
			foreach my $file (@files) {
				if($file=~/^[^\.]/) {
					#dont' copy hidden files
					my $s=basename($file);
					$samplelist{$infolder}{$project}{$s}++;
				}
			}	
		}
	}
	else {
		print STDERR "ERROR:$folder doesn't exist.\n";
		print LOG "ERROR:$folder doesn't exist.\n";
		exit;		
	}
}

#message
foreach my $f (sort keys %samplelist) {
	foreach my $p (sort keys %{$samplelist{$f}}) {
		print STDERR "In Folder: $f, Project: $p, ",scalar(keys %{$samplelist{$f}{$p}})," samples were identified.\n" if $verbose;
		print STDERR "Samples: ",join(",", sort keys %{$samplelist{$f}{$p}}),"\n\n" if $verbose;

		print LOG "In Folder: $f, Project: $p, ",scalar(keys %{$samplelist{$f}{$p}})," samples were identified.\n";
		print LOG "Samples: ",join(",", sort keys %{$samplelist{$f}{$p}}),"\n\n";
	}
}
#output command


open(OUT,">$scriptfile") || die $!;

foreach my $f (sort keys %samplelist) {
	foreach my $p (sort keys %{$samplelist{$f}}) {
		if(!-e "$outputfolder/$p") {
			mkdir "$outputfolder/$p";
		}
		
		foreach my $s (sort keys %{$samplelist{$f}{$p}}) {
		
			if(!-e "$outputfolder/$p/$s") {
				mkdir "$outputfolder/$p/$s";
			}
			
			my $sfolder="$f/Projects/$p/Samples/$s/Files";
			
			#print STDERR $sfolder,"\n";
			my @fastqfiles=glob("'${sfolder}'"."/*fastq.gz"); #deal with special chars in the variable for glob
			
			print STDERR scalar(@fastqfiles)," FASTQ files found for Project $p, Sample $s.\n" if $verbose;
			print LOG scalar(@fastqfiles)," FASTQ files found for Project $p, Sample $s.\n";
			
			my @newfastqfiles_r1;
			my @newfastqfiles_r2;
			foreach my $file (@fastqfiles) {
				print OUT "cp \'$file\' \'$outputfolder/$p/$s\';";
				
				if($merge eq "T") {
					if($file=~/(R\d)_\d+.fastq.gz/) {
						my $filename=basename($file);
						if($1 eq "R1") {
							push @newfastqfiles_r1,"\'$outputfolder/$p/$s/$filename\'";
						}
						elsif($1 eq "R2") {
							push @newfastqfiles_r2,"\'$outputfolder/$p/$s/$filename\'";
						}
						else {
							print STDERR "ERROR:$file doesn't conform naming convention ",'(R\d)_\d+.fastq.gz',".\n";
						}
					}
					else {
						print STDERR "ERROR:$file doesn't conform naming convention ",'(R\d)_\d+.fastq.gz',".\n";
					}
				}
			}
			
			if($merge eq "T") {
				#provide option to merge fastq files
				print OUT "cat ",join(" ",@newfastqfiles_r1)," > \'$mergedfolder/$s\_R1.fastq.gz\';";
				if(@newfastqfiles_r2) {
					print OUT "cat ",join(" ",@newfastqfiles_r2)," > \'$mergedfolder/$s\_R2.fastq.gz\';";
				}
			}		
			print OUT "\n";
			
		}
	}
}
close OUT;


########
#Execuate
########

#local mode
print STDERR "\nTo run locally, in shell type: sh $scriptfile\n\n";
print LOG "\nTo run locally, in shell type: sh $scriptfile\n\n";

if($runmode eq "local") {
	print STDERR "Running locally...\n" if $verbose;
	print LOG "Running locally...\n";
	system("sh $scriptfile");
	print STDERR "Done\n\n" if $verbose;
	print LOG "Done\n\n";
}

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

