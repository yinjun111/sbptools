#!/usr/bin/perl -w
use strict;
use Getopt::Long;



my $version="0.1";



my $usage="

sbptools
version: $version
Usage: sbptools [tool] [parameters]


Parameters:

    
    rnaseq-process    RNA-seq QC, Align, and RSEM for FASTQ files
    rnaseq-merge      Merge rnaseq-process results for downstream analyses
    rnaseq-de         Perform DE analysis using DESeq2
    
	bs-fastq          Download and merge FASTQ files from Basespace	

    ensembl2ucsc      Convert Ensembl gtf/fasta/bed into UCSC format

    mergefiles        Use a model file to merge different files together
    text2excel        Build excel file using text file [not implemented]

";


unless (@ARGV) {
	print STDERR $usage;
	exit;
}

#then call different scripts



my $bs_fastq="perl /apps/sbptools/bs-fastq/bs-fastq_caller.pl";
my $rnaseq_process="perl /apps/sbptools/rnaseq-process/rnaseq-process_caller.pl";
my $rnaseq_merge="perl /apps/sbptools/rnaseq-merge/rnaseq-merge_caller.pl";
my $rnaseq_de="perl /apps/sbptools/rnaseq-de/rnaseq-de_caller.pl";
my $mergefiles="perl /apps/sbptools/mergefiles/mergefiles_caller.pl";
my $text2excel ="perl /apps/sbptools/text2excel/text2excel_caller.pl";


my %commands2program=(
    "bs-fastq"=>$bs_fastq,
    "rnaseq-process"=>$rnaseq_process,
    "rnaseq-merge"=>$rnaseq_merge,
    "rnaseq-de"=>$rnaseq_de,
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
