#!/usr/bin/perl -w
use strict;
use Getopt::Long;



my $version="0.5beta";

#v0.3, runmode implementations in rnaseq and chipseq
#v0.31, add rnaseq-motif
#v0.4, major updates to support Firefly
#v0.41 rnaseq-var is supported in Firefly
#v0.5, new procedure to expand rnaseq processing functions

my $usage="

sbptools
version: $version
Usage: sbptools [tool] [parameters]


Parameters:

    ########
    #Supported in both Firefly and Falco
    ########
    rnaseq-process    RNA-seq QC, Align, and RSEM for FASTQ files
    rnaseq-merge      Merge rnaseq-process results for downstream analyses
    rnaseq-de         Perform DE analysis using DESeq2
    rnaseq-summary    Summarize RNA-Seq DE results

    rnaseq-var        RNA-seq variant calling pipeline
    rnaseq-motif      RNA-seq TFBS motif finding pipeline
    motif-finder      Transcription factor binding motif prediction
	
    mergefiles        Use a model file to merge different files together
    text2excel        Build excel file using text file (by Andrew Hodges)

    ########
    #Supported only in Firefly
    ########
    parallel-job      Batch job submission in Firefly

    ########
    #Supported only in Falco
    ########	

    chipseq-process   ChIP-seq QC, Align, and Peak Calling
    chipseq-merge     Summarize ChIP-Seq results
    chipseq-de        Perform DE analysis for ChIP-Seq
    chipseq-summary   Summarize ChIP-seq DE results
	
    bs-fastq          Download and merge FASTQ files from Basespace	

    ensembl2ucsc      Convert Ensembl gtf/fasta/bed into UCSC format


";


unless (@ARGV) {
	print STDERR $usage;
	exit;
}


#functions

my ($command,@params)=@ARGV;

####
#check whether to use --dev version
####
my $params_list=join(" ",@params);

my $dev=0; #developmental version

if($params_list=~/--dev/) {
	$dev=1;
}


my $sbptoolsfolder="/apps/sbptools/";

#adding --dev switch for better development process
if($dev) {
	$sbptoolsfolder="/home/jyin/Projects/Pipeline/sbptools/";
}


#####
#then call different scripts
#####


my $bs_fastq="$sbptoolsfolder/bs-fastq/bs-fastq_caller.pl";

my $rnaseq_process="$sbptoolsfolder/rnaseq-process/rnaseq-process_caller.pl";
my $rnaseq_merge="$sbptoolsfolder/rnaseq-merge/rnaseq-merge_caller.pl";
my $rnaseq_de="$sbptoolsfolder/rnaseq-de/rnaseq-de_caller.pl";
my $rnaseq_summary="$sbptoolsfolder/rnaseq-summary/rnaseq-summary_caller.pl";

my $rnaseq_var="sh $sbptoolsfolder/rnaseq-var/gatk3_rnaseq_variant_v1.sh";
my $rnaseq_motif="$sbptoolsfolder/rnaseq-motif/rnaseq-motif_caller.pl";

my $chipseq_process="$sbptoolsfolder/chipseq-process/chipseq-process_caller.pl";
my $chipseq_merge="$sbptoolsfolder/chipseq-merge/chipseq-merge_caller.pl";
my $chipseq_de="$sbptoolsfolder/chipseq-de/chipseq-de_caller.pl";
my $chipseq_summary="$sbptoolsfolder/chipseq-summary/chipseq-summary_caller.pl";

my $motif_finder="$sbptoolsfolder/motif-finder/motif-finder_caller.pl";

my $mergefiles="$sbptoolsfolder/mergefiles/mergefiles_caller.pl";
my $text2excel ="perl $sbptoolsfolder/text2excel/text2excel.pl";

my $parallel_job ="perl $sbptoolsfolder/parallel-job/parallel-job_caller.pl";

my %commands2program=(
    "bs-fastq"=>$bs_fastq,
	
    "rnaseq-process"=>$rnaseq_process,
    "rnaseq-merge"=>$rnaseq_merge,
    "rnaseq-de"=>$rnaseq_de,
	"rnaseq-summary"=>$rnaseq_summary,

	"rnaseq-var"=>$rnaseq_var,
	"rnaseq-motif"=>$rnaseq_motif,
	
    "chipseq-process"=>$chipseq_process,
    "chipseq-merge"=>$chipseq_merge,
    "chipseq-de"=>$chipseq_de,
	"chipseq-summary"=>$chipseq_summary,	

	"motif-finder"=>$motif_finder,
    
	"mergefiles"=>$mergefiles,
    "text2excel"=>$text2excel,
	
	"parallel-job"=>$parallel_job,
);



if(defined $commands2program{$command}) {
	system($commands2program{$command}." ".join(" ",@params));
}
else {
	print STDERR "ERORR $command not found in sbptools.\n\n";
	system("sbptools");
}
