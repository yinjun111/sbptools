#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);

#CutAdapt+FASTQC+RSEM+STAR


########
#Interface
########


my $version="0.4";

#v0.2, add filter function to get files for PCA
#v0.3, removed -v, add -r implementation for local
#v0.31, solves screen envinroment problem
#v0.32, add CPM
#v0.4, major updates, to add parallel job, and --dev option

my $usage="

rnaseq-merge
version: $version
Usage: sbptools rnaseq-merge [parameters]

Description: Merge rnaseq-process folder to get summarized QC, count, TPM etc.

Parameters:

    --in|-i           Input folder(s)
	
	#two ways of retrieving samples, either by names or by config files
    --samples|-s      Samples
    --config|-c       Configuration file

    --output|-o       Output folder


    --tx|-t           Transcriptome
                        Current support Human.B38.Ensembl84, Mouse.B38.Ensembl84
    --anno|-a         Add annotation

    --filter          Signal filter [auto]
                         automatically defined signal cutoff as
                           Count >= No. of samples * 5
                         or can define a count number
	
    --runmode|-r      Where to run the scripts, local, cluster or none [none]
	
	
	#Parallel computing controls
    --task            Number of tasks to be paralleled. By default 7 tasks. [7]
    --ncpus           No. of cpus for each task for tasks can't use multiple nodes [2]
    --mem|-m          Memory usage for each process

";

#    --verbose|-v      Verbose

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
my $verbose=1;
my $tx;
my $task=7;
my $ncpus=2;
my $mem;
my $filter="auto";
my $runmode="none";

my $dev=0; #developmental version

GetOptions(
	"in|i=s" => \$inputfolders,
	"samples|s=s" => \$samples,	
	"config|c=s" => \$configfile,
	"output|o=s" => \$outputfolder,
	"tx|t=s" => \$tx,	
	"filter=s" => \$filter,
	"task=s" => \$task,
	"ncpus=s" => \$ncpus,
	"mem=s" => \$mem,	

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


my $multiqc="/apps/python-3.5.2/bin/multiqc";

my $mergefiles="$sbptoolsfolder/mergefiles/mergefiles_caller.pl";
my $parallel_job="$sbptoolsfolder/parallel-job/parallel-job_caller.pl";
my $rnaseqmergefilter="$sbptoolsfolder/rnaseq-merge/rnaseq-merge_filter.R";
my $count2cpm="$sbptoolsfolder/rnaseq-merge/count2cpm.R";



#used programs

my $Rscript=find_program("/apps/R-3.5.2/bin/Rscript");



########
#default ouputs
########


my $genecountmerged="gene.results.merged.count.txt";
my $genetpmmerged="gene.results.merged.tpm.txt";
my $genefpkmmerged="gene.results.merged.fpkm.txt";
my $genecpmmerged="gene.results.merged.cpm.txt";

my $txcountmerged="tx.results.merged.count.txt";
my $txtpmmerged="tx.results.merged.tpm.txt";
my $txfpkmmerged="tx.results.merged.fpkm.txt";
my $txcpmmerged="tx.results.merged.cpm.txt";

#Create folders


if(!-e $outputfolder) {
	mkdir($outputfolder);
}

$outputfolder = abs_path($outputfolder);


my $scriptfolder="$outputfolder/scripts";
my $tempfolder="$outputfolder/temp";

if(!-e $scriptfolder) {
	mkdir($scriptfolder);
}

if(!-e $tempfolder) {
	mkdir($tempfolder);
}

my $logfile="$outputfolder/rnaseq-merge_run.log";

my $scriptfile1="$scriptfolder/rnaseq-merge_run1.sh";
my $scriptfile2="$scriptfolder/rnaseq-merge_run2.sh";

my $scriptlocalrun="$outputfolder/rnaseq-merge_local_submission.sh";
my $scriptclusterrun="$outputfolder/rnaseq-merge_cluster_submission.sh";


#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();
my $timestamp=build_timestamp($now,"long");


print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";
print LOG "$multiqc version:", getsysoutput("$multiqc --version"),"\n";
print LOG "\n";


#test tx option

#my %tx2ref=(
#	"Human.B38.Ensembl84"=>"/home/jyin/Projects/Databases/Ensembl/v84/Human_STAR/Human_RSEM",
#	"Mouse.B38.Ensembl84"=>"/home/jyin/Projects/Databases/Ensembl/v84/Mouse_STAR/Mouse_RSEM",
#);



my %tx2ref=(
	"Human.B38.Ensembl84"=> { 
		"star"=>"/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl84_STAR",
		"chrsize"=>"/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl84_STAR/chrNameLength.txt",
		"fasta"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.dna.primary_assembly_ucsc.fa",
		"gtf"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc.gtf",
		"homeranno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_homeranno.txt",
		"geneanno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_gene_annocombo_rev.txt",
		"txanno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_tx_annocombo.txt"},
	"Mouse.B38.Ensembl84"=>{ 
		"star"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl84_STAR",
		"chrsize"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl84_STAR/chrNameLength.txt",
		"fasta"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.dna.primary_assembly_ucsc.fa",
		"gtf"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc.gtf",
		"homeranno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_homeranno.txt",
		"geneanno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_gene_annocombo_rev.txt",
		"txanno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_tx_anno.txt"}
);



########
#Process
########

print STDERR "\nsbptools rnaseq-merge running ...\n\n" if $verbose;
print LOG "\nsbptools rnaseq-merge running ...\n\n";


if(defined $tx2ref{$tx}) {
	print STDERR "Starting analysis using $tx.\n\n" if $verbose;
	print LOG "Starting analysis using $tx.\n\n";
}
else {
	print STDERR "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n" if $verbose;
	print LOG "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n";
	exit;
}

#open config files to find samples

my @samples_array;
my %samples_hash; #check no dup samples


if(defined $samples && length($samples)>0) {
	foreach my $sample (split(",",$samples)) {
		#keep the order of samples
		push @samples_array, $sample;
		unless (defined $samples_hash{$sample}) {
			$samples_hash{$sample}++
		}
		else {
			print STDERR "ERROR:$sample is defined multiple times in $samples.\n";
			print LOG "ERROR:$sample is defined multiple times in $samples.\n";
			exit;
		}
	}
}
elsif(defined $configfile && length($configfile)>0) {
	#another way of getting samples
	
	open(IN,$configfile) || die "Error reading $configfile. $!";
	#first column should be sample, with a header

	my $fileline=0;
	while(<IN>) {
		tr/\r\n//d;
		my @array=split/\t/;
		
		if($fileline>0) {
			my $sample=$array[0];
			push @samples_array, $sample;
			unless (defined $samples_hash{$sample}) {
				$samples_hash{$sample}++
			}
			else {
				print STDERR "ERROR:$sample is defined multiple times in $configfile.\n";
				print LOG "ERROR:$sample is defined multiple times in $configfile.\n";
				exit;
			}
		}
		
		$fileline++;
	}
	close IN;
}
else {
	#adding error message for no samples defined
	print STDERR "ERROR:Either --config or --samples needs to be defined.\n";
	print LOG "ERROR:Either --config or --samples needs to be defined.\n";

	exit;

}

print STDERR scalar(@samples_array)," samples identified: ",join(",",@samples_array),"\n\n" if $verbose;
print LOG scalar(@samples_array)," samples identified: ",join(",",@samples_array),"\n\n";

#----------------
#Find files for merging
#Examples
#JY_315_Isl.genes.results
#JY_315_Isl.isoforms.results

print STDERR "\nReading sample folders.\n" if $verbose;
print LOG "\nReading sample folders.\n";

my %sample2gene;
my %sample2tx;
my %sample2folder;

foreach my $infolder (split(",",$inputfolders)) {
	#samples in different folders
	#if different folders have the same sample name, use the first found sample
	
	my @infiles=glob("$infolder/*");
	
	foreach my $file (@infiles) {
		if(-d $file) {
			my $samplename=basename($file);
			
			#see whether it is a defined sample
			if(defined $samples_hash{$samplename}) {
				print STDERR $samplename," is found in ",$file,".\n" if $verbose;
				print LOG $samplename," is found in ",$file,".\n" if $verbose;
				
				$sample2folder{$samplename}=abs_path($file);
				
				my @samplefiles=glob("$file/*");
				
				foreach  my $samplefile (@samplefiles) {
					#gene rsem
					if($samplefile=~/.genes.results/) {
						$sample2gene{$samplename}=abs_path($samplefile);
					}
					
					#tx rsem
					if($samplefile=~/.isoforms.results/) {
						$sample2tx{$samplename}=abs_path($samplefile);
					}
				}
			}
		}
	}
}



print STDERR scalar(keys %sample2gene)," samples identified with gene results.\n" if $verbose;
print LOG scalar(keys %sample2gene)," samples identified with gene results.\n";

print STDERR scalar(keys %sample2tx)," samples identified with isoform results.\n" if $verbose;
print LOG scalar(keys %sample2tx)," samples identified with isoform results.\n";


if( scalar(@samples_array) != scalar(keys %sample2gene) || scalar(@samples_array) != scalar(keys %sample2tx)) {
	print STDERR "ERROR:Not all samples have gene or isoform results.\n\n";
	print LOG "ERROR:Not all samples have gene or isoform results.\n\n";
	exit;
}


#print OUT model to combine gene and tx results

system("cut -f 1 ".$sample2gene{$samples_array[0]}." > $tempfolder/genes.txt");
system("cut -f 1 ".$sample2tx{$samples_array[0]}." > $tempfolder/txs.txt");

#print title for gene and txs
open(OUT,">$tempfolder/gene_title.txt") || die $!;
print OUT "Gene\t",join("\t",@samples_array),"\n";
close OUT;

open(OUT,">$tempfolder/tx_title.txt") || die $!;
print OUT "Tx\t",join("\t",@samples_array),"\n";
close OUT;


#get all files
my @genefiles;
my @txfiles;

foreach my $sample (@samples_array) {
	if(defined $sample2gene{$sample}) {
		push @genefiles,$sample2gene{$sample};
	}
	else {
		print STDERR "ERROR: $sample gene.results not defined in $inputfolders.\n\n";
		print LOG "ERROR: $sample gene.results not defined in $inputfolders.\n\n";
		exit;
	}

	if(defined $sample2tx{$sample}) {
		push @txfiles,$sample2tx{$sample};
	}
	else {
		print STDERR "ERROR: $sample isoform.results not defined in $inputfolders.\n\n";
		print LOG "ERROR: $sample isoform.results not defined in $inputfolders.\n\n";
		exit;
	}
}


########
#Print out commands, for local and cluster run
########

open(S1,">$scriptfile1") || die "Error writing $scriptfile1. $!";


#----------------
#Multiqc with specific folders


#file for multiqc folder search
open(OUT,">$tempfolder/samplefolders.txt") || die $!;
foreach my $sample (@samples_array) {
	print OUT $sample2folder{$sample},"\n";
}
close OUT;

print S1 "$multiqc -l $tempfolder/samplefolders.txt -o $outputfolder/multiqc;\n";

#------------------
#may need to implement annotation ...
#need to implement filter step for PCA ready file

#copy gene annotation
print S1 "cp ".$tx2ref{$tx}{"geneanno"}." $outputfolder/geneanno.txt;";
print S1 "cp ".$tx2ref{$tx}{"txanno"}." $outputfolder/txanno.txt;";


#Gene count
print S1 "$mergefiles -m $tempfolder/genes.txt -i ",join(",",@genefiles)," -l 5 -o $tempfolder/$genecountmerged\_wrongtitle;";
print S1 "tail -n +2 $tempfolder/$genecountmerged\_wrongtitle > $tempfolder/$genecountmerged\_notitle;";
print S1 "cat $tempfolder/gene_title.txt $tempfolder/$genecountmerged\_notitle > $outputfolder/$genecountmerged;";
print S1 "$Rscript $count2cpm --count $outputfolder/$genecountmerged --cpm $outputfolder/$genecpmmerged;";
print S1 "rm $tempfolder/$genecountmerged\_wrongtitle;rm $tempfolder/$genecountmerged\_notitle;";
print S1 "\n";

#Gene tpm
print S1 "$mergefiles -m $tempfolder/genes.txt -i ",join(",",@genefiles)," -l 6 -o $tempfolder/$genetpmmerged\_wrongtitle -v;";
print S1 "tail -n +2 $tempfolder/$genetpmmerged\_wrongtitle > $tempfolder/$genetpmmerged\_notitle;";
print S1 "cat $tempfolder/gene_title.txt $tempfolder/$genetpmmerged\_notitle > $outputfolder/$genetpmmerged;";
print S1 "rm $tempfolder/$genetpmmerged\_wrongtitle;rm $tempfolder/$genetpmmerged\_notitle;";
print S1 "\n";

#Gene fpkm
print S1 "$mergefiles -m $tempfolder/genes.txt -i ",join(",",@genefiles)," -l 7 -o $tempfolder/$genefpkmmerged\_wrongtitle -v;";
print S1 "tail -n +2 $tempfolder/$genefpkmmerged\_wrongtitle > $tempfolder/$genefpkmmerged\_notitle;";
print S1 "cat $tempfolder/gene_title.txt $tempfolder/$genefpkmmerged\_notitle > $outputfolder/$genefpkmmerged;";
print S1 "rm $tempfolder/$genefpkmmerged\_wrongtitle;rm $tempfolder/$genefpkmmerged\_notitle;";
print S1 "\n";

#tx count
print S1 "$mergefiles -m $tempfolder/txs.txt -i ",join(",",@txfiles)," -l 5 -o $tempfolder/$txcountmerged\_wrongtitle -v;";
print S1 "tail -n +2 $tempfolder/$txcountmerged\_wrongtitle > $tempfolder/$txcountmerged\_notitle;";
print S1 "cat $tempfolder/tx_title.txt $tempfolder/$txcountmerged\_notitle > $outputfolder/$txcountmerged;";
print S1 "$Rscript $count2cpm --count $outputfolder/$txcountmerged --cpm $outputfolder/$txcpmmerged;";
print S1 "rm $tempfolder/$txcountmerged\_wrongtitle;rm $tempfolder/$txcountmerged\_notitle;";
print S1 "\n";

#tx tpm
print S1 "$mergefiles -m $tempfolder/txs.txt -i ",join(",",@txfiles)," -l 6 -o $tempfolder/$txtpmmerged\_wrongtitle -v;";
print S1 "tail -n +2 $tempfolder/$txtpmmerged\_wrongtitle > $tempfolder/$txtpmmerged\_notitle;";
print S1 "cat $tempfolder/tx_title.txt $tempfolder/$txtpmmerged\_notitle > $outputfolder/$txtpmmerged;";
print S1 "rm $tempfolder/$txtpmmerged\_wrongtitle;rm $tempfolder/$txtpmmerged\_notitle;";
print S1 "\n";

#tx fpkm
print S1 "$mergefiles -m $tempfolder/txs.txt -i ",join(",",@txfiles)," -l 7 -o $tempfolder/$txfpkmmerged\_wrongtitle -v;";
print S1 "tail -n +2 $tempfolder/$txfpkmmerged\_wrongtitle > $tempfolder/$txfpkmmerged\_notitle;";
print S1 "cat $tempfolder/tx_title.txt $tempfolder/$txfpkmmerged\_notitle > $outputfolder/$txfpkmmerged;";
print S1 "rm $tempfolder/$txfpkmmerged\_wrongtitle;rm $tempfolder/$txfpkmmerged\_notitle;";
print S1 "\n";


close S1;


##########
#filter for PCA
##########

open(S2,">$scriptfile2") || die "Error writing $scriptfile2. $!";

#gene
print S2 "$Rscript $rnaseqmergefilter --count $outputfolder/$genecountmerged --fpkm $outputfolder/$genefpkmmerged --tpm $outputfolder/$genetpmmerged --filter $filter;";

my $genecountmerged_filtered="gene.results.merged.count.filtered.$filter.txt";
my $genecpmmerged_filtered="gene.results.merged.cpm.filtered.$filter.txt";

print S2 "$Rscript $count2cpm --count $outputfolder/$genecountmerged_filtered --cpm $outputfolder/$genecpmmerged_filtered;";
print S2 "\n";

#tx
print S2 "$Rscript $rnaseqmergefilter --count $outputfolder/$txcountmerged --fpkm $outputfolder/$txfpkmmerged --tpm $outputfolder/$txtpmmerged --filter $filter;\n";

my $txcountmerged_filtered="tx.results.merged.count.filtered.$filter.txt";
my $txcpmmerged_filtered="tx.results.merged.cpm.filtered.$filter.txt";

print S2 "$Rscript $count2cpm --count $outputfolder/$txcountmerged_filtered --cpm $outputfolder/$txcpmmerged_filtered;";

print S2 "\n";

close S2;




#######
#Run mode
#######

open(LOUT,">$scriptlocalrun") || die "ERROR:can't write to $scriptlocalrun. $!";
open(SOUT,">$scriptclusterrun") || die "ERROR:can't write to $scriptclusterrun. $!";


my @scripts_all=($scriptfile1,$scriptfile2);


#print out command for local and cluster parallel runs
my $jobnumber=0;
my $jobname="rnaseq-merge-$timestamp";

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
		#push @script_names,$1."_".basename_short($outputfolder);
		push @script_names,$1;
	}
}

my $localcommand="screen -S $jobname -dm bash -c \"source ~/.bashrc;".join(";",@local_runs).";\"";

print LOUT $localcommand,"\n";
close LOUT;

#print out command for cluster parallel runs

my $clustercommand="perl $parallel_job -i ".join(",", @scripts_all)." -o $scriptfolder -n ".join(",",@script_names)." --tandem -t $task --ncpus $ncpus -r ";

if(defined $mem && length($mem)>0) {
	$clustercommand.=" -m $mem";
}

print SOUT $clustercommand,"\n";
close SOUT;



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

sub getsysoutput {
	my $command=shift @_;
	my $output=`$command`;
	$output=~tr/\r\n//d;
	return $output;
}

sub basename_short {
	my $filename=shift @_;
	
	my $basename;
	
	if($filename=~/([^\/]+)\/?$/) {
		$basename=$1;
	}
	
	return $basename;

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


