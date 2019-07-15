#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);

#CutAdapt+FASTQC+RSEM+STAR


########
#Prerequisites
########

my $multiqc="/apps/python-3.5.2/bin/multiqc";
my $mergefiles="/apps/sbptools/mergefiles/mergefiles_caller.pl";

my $rnaseqmergefilter="/apps/sbptools/rnaseq-merge/rnaseq-merge_filter.R";


########
#Interface
########


my $version="0.2";

#v0.2, add filter function to get files for PCA


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
my $filter="auto";
my $runmode="none";

GetOptions(
	"in|i=s" => \$inputfolders,
	"samples|s=s" => \$samples,	
	"config|c=s" => \$configfile,
	"output|o=s" => \$outputfolder,
	"tx|t=s" => \$tx,	
	"filter=s" => \$filter,	
	"runmode|r=s" => \$runmode,		
	"verbose|v" => \$verbose,
);


#default ouputs

my $genecountmerged="gene.results.merged.count.txt";
my $genetpmmerged="gene.results.merged.tpm.txt";
my $genefpkmmerged="gene.results.merged.fpkm.txt";

my $txcountmerged="tx.results.merged.count.txt";
my $txtpmmerged="tx.results.merged.tpm.txt";
my $txfpkmmerged="tx.results.merged.fpkm.txt";


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

#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";
print LOG "$multiqc version:", getsysoutput("$multiqc --version"),"\n";
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


#Files for merging
#Examples
#JY_315_Isl.genes.results
#JY_315_Isl.isoforms.results


########
#Process
########

print STDERR "\nsbptools rnaseq-merge running ...\n\n" if $verbose;
print LOG "\nsbptools rnaseq-merge running ...\n\n";

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
#Print out commands, for local and server run
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


#Gene count
print S1 "$mergefiles -m $tempfolder/genes.txt -i ",join(",",@genefiles)," -l 5 -o $tempfolder/$genecountmerged\_wrongtitle;";
print S1 "tail -n +2 $tempfolder/$genecountmerged\_wrongtitle > $tempfolder/$genecountmerged\_notitle;";
print S1 "cat $tempfolder/gene_title.txt $tempfolder/$genecountmerged\_notitle > $outputfolder/$genecountmerged;";
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


#filter for PCA
open(S2,">$scriptfile2") || die "Error writing $scriptfile2. $!";

print S2 "Rscript $rnaseqmergefilter --count $outputfolder/$genecountmerged --fpkm $outputfolder/$genefpkmmerged --tpm $outputfolder/$genetpmmerged --filter $filter;\n";
print S2 "Rscript $rnaseqmergefilter --count $outputfolder/$txcountmerged --fpkm $outputfolder/$txfpkmmerged --tpm $outputfolder/$txtpmmerged --filter $filter;\n";

close S2;


#local mode
print STDERR "\nTo run locally, in shell type: sh $scriptfile1;sh $scriptfile2\n\n";
print LOG "\nTo run locally, in shell type: sh $scriptfile1;sh $scriptfile2\n\n";



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


