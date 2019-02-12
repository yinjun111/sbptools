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
    mergefiles        Use a model file to merge different files together
    text2excel        Build excel file using text file

";


unless (@ARGV) {
	print STDERR $usage;
	exit;
}

#then call different scripts



