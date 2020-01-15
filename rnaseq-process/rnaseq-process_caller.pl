#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);

#CutAdapt+FASTQC+RSEM+STAR


########
#Interface
########


my $version="0.41";

#0.2b change ensembl to UCSC format
#0.2c add bw generation
#0.2d correct bug for PE
#0.3 add runmode, always show -v
#v0.31, solves screen envinroment problem
#v0.4 add find_program, and queuejob, --dev switch, add 30gb requirement
#v0.41 option to turn off bamcoverage due to long processing time. changed cutadapt logging

my $usage="

rnaseq-process
version: $version
Usage: sbptools rnaseq-process [parameters]

Description: 


Parameters:

    --config|-c       Configuration file
                           first column as sample name.
                           Controled headers include:
                               fastq1,fastq2 columns for fastqs
                               index for sample index used for folder names
    --output|-o       Output folder

    --tx|-t           Transcriptome
                        Current support Human.B38.Ensembl84, Mouse.B38.Ensembl84

    --task            Number of tasks to be paralleled. By default 4 tasks for local mode, 8 tasks for cluster mode.
    --ncpus           No. of cpus for each task [4]
    --mem|-m          Memory usage for each process, e.g. 100mb, 100gb [40gb]

    --nobamcoverage   Use --nobamcoverage turn off producing bw file for bam files


    --runmode|-r      Where to run the scripts, local, cluster or none [none]
                            local is to run locally using \"parallel\", recommended for Falco
                            cluster is to submit jobs to PBS queue in the HPC, recommended for Firefly
                            none is to generate scripts only, after that,
                                   you can use \"sh rnaseq-merge_local_submission.sh\" to run locally, or
                                   you can use \"sh rnaseq-merge_cluster_submission.sh\" to submit job to PBS

								   
";

#--verbose|-v      Verbose

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
my $verbose=1;
my $tx;
my $task;
my $ncpus=4;
my $nobamcoverage=0;
my $mem="40gb";
my $runmode="none";

my $dev=0; #developmental version


GetOptions(
	"config|c=s" => \$configfile,
	"output|o=s" => \$outputfolder,
	"tx|t=s" => \$tx,
	"task=s" => \$task,
	"ncpus=s" => \$ncpus,
	"mem=s" => \$mem,
	"nobamcoverage" => \$nobamcoverage,	
	"runmode|r=s" => \$runmode,		
	"verbose|v" => \$verbose,
	"dev" => \$dev,		
);


#tasks, local 4, cluster 8
unless(defined $task && length($task)>0) {
	if($runmode eq "cluster") {
		$task=8;
	}
	else {
		$task=4;
	}
}



########
#Prerequisites
########



my $sbptoolsfolder="/apps/sbptools/";

#adding --dev switch for better development process
if($dev) {
	$sbptoolsfolder="/home/jyin/Projects/Pipeline/sbptools/";
}

my $mergefiles="$sbptoolsfolder/mergefiles/mergefiles_caller.pl";
my $parallel_job="$sbptoolsfolder/parallel-job/parallel-job_caller.pl";
my $rnaseqmergefilter="$sbptoolsfolder/rnaseq-merge/rnaseq-merge_filter.R";
my $count2cpm="$sbptoolsfolder/rnaseq-merge/count2cpm.R";


my $cutadapt=find_program("/apps/python-3.5.2/bin/cutadapt");
my $fastqc=find_program("/apps/FastQC/fastqc");
my $rsem=find_program("/apps/RSEM-1.3.1/rsem-calculate-expression");
my $star=find_program("/apps/STAR-master/bin/Linux_x86_64/STAR");
my $bamcoverage=find_program("/apps/python-3.5.2/bin/bamCoverage");



#######
#Output folder
#######

unless(defined $outputfolder && length($outputfolder)>0 ) {
	print STDERR "\nERROR: -o outputfolder needs to be defined without default value.\n\n";
	exit;
}

if(!-e $outputfolder) {
	mkdir($outputfolder);
}

$outputfolder = abs_path($outputfolder);

my $scriptfolder="$outputfolder/scripts";

if(!-e $scriptfolder) {
	mkdir($scriptfolder);
}


my $logfile="$outputfolder/rnaseq-process_run.log";
my $newconfigfile="$outputfolder/rnaseq-process_config.txt";

my $scriptfile1="$scriptfolder/rnaseq-process_run.sh";


my $scriptlocalrun="$outputfolder/rnaseq-process_local_submission.sh";
my $scriptclusterrun="$outputfolder/rnaseq-process_cluster_submission.sh";



#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();
my $timestamp=build_timestamp($now,"long");

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
	"Human.B38.Ensembl84"=>"/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl84_STAR/Human_RSEM",
	"Mouse.B38.Ensembl84"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl84_STAR/Mouse_RSEM",
);


########
#Process
########

print STDERR "\nsbptools rnaseq-process $version running ...\n\n" if $verbose;
print LOG "\nsbptools rnaseq-process $version running ...\n\n";


if(defined $tx2ref{$tx}) {
	print STDERR "Starting analysis using $tx.\n\n" if $verbose;
	print LOG "Starting analysis using $tx.\n\n";
}
else {
	print STDERR "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n" if $verbose;
	print LOG "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n";
	exit;
}


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
		
		#should always use first column as sample name and index to maintain continuous workflow
		
		#first column as sample
		
		#check title
		unless(defined $configattrs{"FASTQ1"}) {
			#SAMPLE and FASTQ1 have to be defined. The others are optional
			print STDERR "ERROR: FASTQ1 need to be defined in $configfile. Input config includes:",join(",",map {uc $_} @array),"\n";
			print LOG "ERROR: FASTQ1 need to be defined in $configfile. Input config includes:",join(",",map {uc $_} @array),"\n";
			exit;
		}
		else {
			print STDERR "Input config $configfile includes:",join(",",map {uc $_} @array),"\n\n" if $verbose;
			print LOG "Input config $configfile includes:",join(",",map {uc $_} @array),"\n\n";
		}
		
		#version 0.2, removed self defined INDEX
		
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
		#if defined INDEX, use index, otherwise, use first column
		if(defined $configattrs{"INDEX"}) {
			$indexname=$array[$configattrs{"INDEX"}];
		}
		else {
			$indexname=$array[0]; #change first column as sample column in v0.2
			$indexname=~s/[^\w-]/_/g;
		}
		
		
		$sample2indexname{$array[0]}=$indexname;

		
		#fastq1
		if(-e $array[$configattrs{"FASTQ1"}]) {
			push @{$sample2fastq{$indexname}},$array[$configattrs{"FASTQ1"}];
		}
		else {
			print STDERR "\n\nERROR:",$array[$configattrs{"FASTQ1"}]," doesn't exist. Please check $configfile setting.\n\n";
			print LOG "\n\nERROR:",$array[$configattrs{"FASTQ1"}]," doesn't exist. Please check $configfile setting.\n\n";
		}
		
		#fastq2 (Read2)
		if(defined $configattrs{"FASTQ2"}) {
			if(-e $array[$configattrs{"FASTQ2"}]) {
				push @{$sample2fastq{$indexname}},$array[$configattrs{"FASTQ2"}];
			}
			else {
				print STDERR "\n\nERROR:",$array[$configattrs{"FASTQ2"}]," doesn't exist. Please check $configfile setting.\n\n";
				print LOG "\n\nERROR:",$array[$configattrs{"FASTQ2"}]," doesn't exist. Please check $configfile setting.\n\n";	
			}
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
print STDERR scalar(keys %sample2indexname)," samples identified from $configfile, including:\nINDEX\tSAMPLE\n",join("\n",map {$sample2indexname{$_}."\t".$_} sort keys %sample2indexname),"\n\n" if $verbose;
print LOG scalar(keys %sample2indexname)," samples identified from $configfile, including:\nINDEX\tSAMPLE\n",join("\n",map {$sample2indexname{$_}."\t".$_} sort keys %sample2indexname),"\n\n" if $verbose;

#print STDERR scalar(keys %sample2fastq)," samples identified from $configfile, including:",join("\n",sort keys %sample2fastq),"\n\n" if $verbose;
#print LOG scalar(keys %sample2fastq)," samples identified from $configfile, including:",join("\n",sort keys %sample2fastq),"\n\n";


print STDERR "Make folders for different samples.\n\n" if $verbose;
print LOG "Make folders for different samples.\n\n";

foreach my $sample (sort keys %sample2fastq) {
	my $samplefolder="$outputfolder/$sample"; #edited to enable space in folder name
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
		my $fastq2=$sample2fastq{$sample}[1]; #
		
		my $fastq1trim="$samplefolder/".basename($fastq1);
		$fastq1trim=~s/\.fastq\.gz/_trimmed.fastq.gz/;
		push @{$sample2fastq{$sample}},$fastq1trim;

		my $fastq2trim="$samplefolder/".basename($fastq2);
		$fastq2trim=~s/\.fastq\.gz/_trimmed.fastq.gz/;
		push @{$sample2fastq{$sample}},$fastq2trim;
		
		my $cutadaptlog="$samplefolder/$sample\_cutadapt.log";

		#$sample2workflow{$sample}.="$cutadapt -j 4 -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o $fastq1trim -p $fastq2trim $fastq1 $fastq2 > $cutadaptlog 2>&1;";

		#$sample2workflow{$sample}.="$cutadapt -j 4 -m 20 --interleaved -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT $fastq1 $fastq2 | $cutadapt --interleaved -j 4 -m 20 -a \"A{100}\" -A \"A{100}\" - | $cutadapt --interleaved -j 4 -m 20 -a \"T{100}\" -A \"T{100}\" - -o $fastq1trim -p $fastq2trim >> $cutadaptlog 2>&1;";

		
		#add poly-A/T trimming, used shorter -A adapter trimming
		#$sample2workflow{$sample}.="$cutadapt -j 4 -m 20 --interleaved -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT $fastq1 $fastq2 | $cutadapt --interleaved -j 4 -m 20 -a \"A{100}\" -A \"A{100}\" - | $cutadapt --interleaved -j 4 -m 20 -a \"T{100}\" -A \"T{100}\" - -o $fastq1trim -p $fastq2trim >> $cutadaptlog 2>&1;";
		
		$sample2workflow{$sample}.="$cutadapt -j 4 -m 20 --interleaved -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT $fastq1 $fastq2 2>>$cutadaptlog | $cutadapt --interleaved -j 4 -m 20 -a \"A{100}\" -A \"A{100}\" - 2>>$cutadaptlog | $cutadapt --interleaved -j 4 -m 20 -a \"T{100}\" -A \"T{100}\" - -o $fastq1trim -p $fastq2trim 1>>$cutadaptlog;";
		
		
	}
}
else {
	#SE
	print STDERR "Printing cutadapt SE script.\n\n" if $verbose;
	print LOG "Printing cutadapt SE script.\n\n";
	
	
	
	foreach my $sample (sort keys %sample2fastq) {
		my $samplefolder="$outputfolder/$sample";
		
		my $cutadaptlog="$samplefolder/$sample\_cutadapt.log";
		
		my $fastq1=$sample2fastq{$sample}[0];
		
		my $fastq1trim="$samplefolder/".basename($fastq1);
		$fastq1trim=~s/\.fastq\.gz/_trimmed.fastq.gz/;
		push @{$sample2fastq{$sample}},$fastq1trim;

		#$sample2workflow{$sample}.="$cutadapt -j 4 -m 20 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o $fastq1trim $fastq1 > $cutadaptlog 2>&1;";
		
		#add poly-A/T trimming
		#$sample2workflow{$sample}.="$cutadapt -j 4 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC $fastq1 | $cutadapt -j 4 -a \"A{100}\" - | $cutadapt -j 4 -m 20 -a \"T{100}\" - -o $fastq1trim >> $cutadaptlog 2>&1;";
		
		$sample2workflow{$sample}.="$cutadapt -j 4 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC $fastq1 2>>$cutadaptlog | $cutadapt -j 4 -a \"A{100}\" - 2>>$cutadaptlog | $cutadapt -j 4 -m 20 -a \"T{100}\" - -o $fastq1trim 1>>$cutadaptlog;";

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

	#     rsem-calculate-expression [options] --paired-end upstream_read_file(s) downstream_read_file(s) reference_name sample_name 
	
	foreach my $sample (sort keys %sample2fastq) {
		my $samplefolder="$outputfolder/$sample";
		
		my $rsemlog="$samplefolder/$sample\_rsem.log";
		
		$sample2workflow{$sample}.="$rsem -p 4 --output-genome-bam --sort-bam-by-coordinate --star-gzipped-read-file --star --paired-end ".$sample2fastq{$sample}[2]." ".$sample2fastq{$sample}[3]." ".$tx2ref{$tx}." $samplefolder/$sample > $rsemlog 2>&1;";
		
		if($nobamcoverage==0) {
			$sample2workflow{$sample}.="$bamcoverage --numberOfProcessors 4 --bam $samplefolder/$sample.genome.sorted.bam --normalizeUsing CPM --binSize 5 -o $samplefolder/$sample.genome.sorted.bw;";
		}
	}
}
else {
	#SE
	print STDERR "Printing RSEM SE script.\n\n" if $verbose;
	print LOG "Printing RSEM SE script.\n\n";

	#     rsem-calculate-expression [options] upstream_read_file(s) reference_name sample_name 
	
	foreach my $sample (sort keys %sample2fastq) {
		my $samplefolder="$outputfolder/$sample";
		my $rsemlog="$samplefolder/$sample\_rsem.log";
		
		$sample2workflow{$sample}.="$rsem -p 4 --output-genome-bam --sort-bam-by-coordinate --star-gzipped-read-file --star ".$sample2fastq{$sample}[1]." ".$tx2ref{$tx}." $samplefolder/$sample > $rsemlog 2>&1;";
		
		if($nobamcoverage==0) {
			$sample2workflow{$sample}.="$bamcoverage --numberOfProcessors 4 --bam $samplefolder/$sample.genome.sorted.bam --normalizeUsing CPM --binSize 5 -o $samplefolder/$sample.genome.sorted.bw;";
		}
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


#######
#Run mode
#######

open(LOUT,">$scriptlocalrun") || die "ERROR:can't write to $scriptlocalrun. $!";
open(SOUT,">$scriptclusterrun") || die "ERROR:can't write to $scriptclusterrun. $!";


my @scripts_all=($scriptfile1);


#print out command for local and cluster parallel runs
my $jobnumber=0;
my $jobname="rnaseq-process-$timestamp";

if($task eq "auto") {
	$jobnumber=0;
}
else {
	$jobnumber=$task;
}

my @local_runs;
my @script_names;

foreach my $script (@scripts_all) {
	push @local_runs,"cat $script | parallel -j $jobnumber";

	if($script=~/([^\/]+)\.\w+$/) {
		push @script_names,$1;
	}
}

my $localcommand="screen -S $jobname -dm bash -c \"source ~/.bashrc;".join(";",@local_runs).";\"";

print LOUT $localcommand,"\n";
close LOUT;

#print out command for cluster parallel runs

my $clustercommand="perl $parallel_job -i ".join(",", @scripts_all)." -o $scriptfolder -n ".join(",",@script_names)." --tandem -t $task --ncpus $ncpus --env"; #changed here for none version

if(defined $mem && length($mem)>0) {
	$clustercommand.=" -m $mem";	
}

print SOUT "sh $outputfolder/scripts/parallel-job_submit.sh\n"; #submit step
close SOUT;


system("$clustercommand");
print LOG "$clustercommand;\n\n";



if($runmode eq "none") {
	print STDERR "\nTo run locally, in shell type: sh $scriptlocalrun\n";
	print STDERR "To run in cluster, in shell type: sh $scriptclusterrun\n";
	
	print LOG "\nTo run locally, in shell type: sh $scriptlocalrun\n";
	print LOG "To run in cluster, in shell type: sh $scriptclusterrun\n";
}
elsif($runmode eq "local") {
	#local mode
	#implemented for Falco
	
	system("sh $scriptlocalrun");
	print LOG "sh $scriptlocalrun;\n\n";

	print STDERR "Starting local paralleled processing using $jobnumber tasks. To monitor process, use \"screen -r $jobname\".\n\n";
	print LOG "Starting local paralleled processing using $jobnumber tasks. To monitor process, use \"screen -r $jobname\".\n\n";
	
}
elsif($runmode eq "cluster") {
	#cluster mode
	#implement for Firefly
	
	system("sh $scriptclusterrun");
	print LOG "sh $scriptclusterrun;\n\n";

	print STDERR "Starting cluster paralleled processing using $jobnumber tasks. To monitor process, use \"qstat\".\n\n";

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






