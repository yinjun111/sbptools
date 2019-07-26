#!/usr/bin/perl -w
use strict;
use Getopt::Long;



my $version="0.3";

#v0.3, runmode implementations in rnaseq and chipseq

my $usage="

sbptools
version: $version
Usage: sbptools [tool] [parameters]


Parameters:

    
    rnaseq-process    RNA-seq QC, Align, and RSEM for FASTQ files
    rnaseq-merge      Merge rnaseq-process results for downstream analyses
    rnaseq-de         Perform DE analysis using DESeq2
    rnaseq-summary    Summarize RNA-Seq DE results

    rnaseq-var        RNA-seq variant calling pipeline
    
    chipseq-process   ChIP-seq QC, Align, and Peak Calling
    chipseq-merge     Summarize ChIP-Seq results
    chipseq-de        Perform DE analysis for ChIP-Seq
    chipseq-summary   Summarize ChIP-seq DE results

    motif-finder      Transcription factor binding motif prediction
	
    bs-fastq          Download and merge FASTQ files from Basespace	

    ensembl2ucsc      Convert Ensembl gtf/fasta/bed into UCSC format

    mergefiles        Use a model file to merge different files together
    text2excel        Build excel file using text file (by Andrew Hodges)

";


unless (@ARGV) {
	print STDERR $usage;
	exit;
}

#then call different scripts



my $bs_fastq="/apps/sbptools/bs-fastq/bs-fastq_caller.pl";

my $rnaseq_process="/apps/sbptools/rnaseq-process/rnaseq-process_caller.pl";
my $rnaseq_merge="/apps/sbptools/rnaseq-merge/rnaseq-merge_caller.pl";
my $rnaseq_de="/apps/sbptools/rnaseq-de/rnaseq-de_caller.pl";
my $rnaseq_summary="/apps/sbptools/rnaseq-summary/rnaseq-summary_caller.pl";

my $rnaseq_var="sh /apps/sbptools/rnaseq-var/gatk3_rnaseq_variant_v1.sh";

my $chipseq_process="/apps/sbptools/chipseq-process/chipseq-process_caller.pl";
my $chipseq_merge="/apps/sbptools/chipseq-merge/chipseq-merge_caller.pl";
my $chipseq_de="/apps/sbptools/chipseq-de/chipseq-de_caller.pl";
my $chipseq_summary="/apps/sbptools/chipseq-summary/chipseq-summary_caller.pl";

my $motif_finder="/apps/sbptools/motif-finder/motif-finder_caller.pl";

my $mergefiles="/apps/sbptools/mergefiles/mergefiles_caller.pl";
my $text2excel ="perl /apps/sbptools/text2excel/text2excel.pl";


my %commands2program=(
    "bs-fastq"=>$bs_fastq,
	
    "rnaseq-process"=>$rnaseq_process,
    "rnaseq-merge"=>$rnaseq_merge,
    "rnaseq-de"=>$rnaseq_de,
	"rnaseq-summary"=>$rnaseq_summary,

	"rnaseq-var"=>$rnaseq_var,
	
    "chipseq-process"=>$chipseq_process,
    "chipseq-merge"=>$chipseq_merge,
    "chipseq-de"=>$chipseq_de,
	"chipseq-summary"=>$chipseq_summary,	

	"motif-finder"=>$motif_finder,
    
	"mergefiles"=>$mergefiles,
    "text2excel"=>$text2excel,
);

my ($command,@params)=@ARGV;


if(defined $commands2program{$command}) {
	system($commands2program{$command}." ".join(" ",@params));
}
else {
	print STDERR "ERORR $command not found in sbptools.\n\n";
	system("sbptools");
}
