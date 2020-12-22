#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);

my $version="0.62";

#v0.3, runmode implementations in rnaseq and chipseq
#v0.31, add rnaseq-motif
#v0.4, major updates to support Firefly
#v0.41 rnaseq-var is supported in Firefly
#v0.5, new procedure to expand rnaseq processing functions
#v0.51, update rnaseq-var to v2
#v0.52, start to support different sbptools versions
#v0.6, major updates planned for R4.0, chipseq Firefly compatibility, Ensembl v100, de novo assembler, git compatibility
#v0.61, add geo-download
#v0.62, add dnaseq-process, rnaseq-var improved

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
    rnaseq-motif-summary  RNA-seq TFBS motif finding results summary	
    motif-finder      Transcription factor binding motif prediction
	
    chipseq-process   ChIP-seq QC, Align, and Peak Calling
    chipseq-merge     Summarize ChIP-Seq results
    chipseq-de        Perform DE analysis for ChIP-Seq
    chipseq-summary   Summarize ChIP-seq DE results

    dnaseq-process    DNA-Seq (Exome/Genome-Seq) processing based on GATK4

    mergefiles        Use a model file to merge different files together
    text2excel        Build excel file using text file (by Andrew Hodges)

    gsea-gen          Generate files ready for GSEA analysis

    ########
    #Supported only in Firefly
    ########
    parallel-job      Batch job submission in Firefly

    ########
    #Supported only in Falco
    ########	
	
    bs-fastq          Download and merge FASTQ files from Basespace	
    ensembl2ucsc      Convert Ensembl gtf/fasta/bed into UCSC format
    geo-download      Download raw FASTQ files from GEO

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
else {
	#the tools called will be within the same folder of the script
	$sbptoolsfolder=dirname(abs_path($0));
}



#####
#then call different scripts
#####


my $bs_fastq="$sbptoolsfolder/bs-fastq/bs-fastq_caller.pl";
my $geo_download="$sbptoolsfolder/geo-download/geo-download_caller.pl";

my $rnaseq_process="$sbptoolsfolder/rnaseq-process/rnaseq-process_caller.pl";
my $rnaseq_merge="$sbptoolsfolder/rnaseq-merge/rnaseq-merge_caller.pl";
my $rnaseq_de="$sbptoolsfolder/rnaseq-de/rnaseq-de_caller.pl";
my $rnaseq_summary="$sbptoolsfolder/rnaseq-summary/rnaseq-summary_caller.pl";

my $rnaseq_var="sh $sbptoolsfolder/rnaseq-var/rnaseq-var_caller.sh";
my $rnaseq_var_summary="$sbptoolsfolder/rnaseq-var/rnaseq-var_summary.pl";
my $rnaseq_motif="$sbptoolsfolder/rnaseq-motif/rnaseq-motif_caller.pl";
my $rnaseq_motif_summary="$sbptoolsfolder/rnaseq-motif-summary/rnaseq-motif-summary.pl";

my $gsea_gen="perl $sbptoolsfolder/gsea-gen/gsea-gen_caller.pl";

my $chipseq_process="$sbptoolsfolder/chipseq-process/chipseq-process_caller.pl";
my $chipseq_merge="$sbptoolsfolder/chipseq-merge/chipseq-merge_caller.pl";
my $chipseq_de="$sbptoolsfolder/chipseq-de/chipseq-de_caller.pl";
my $chipseq_summary="$sbptoolsfolder/chipseq-summary/chipseq-summary_caller.pl";

my $dnaseq_process="$sbptoolsfolder/dnaseq-process/dnaseq-process_caller.pl";

my $motif_finder="$sbptoolsfolder/motif-finder/motif-finder_caller.pl";

my $mergefiles="$sbptoolsfolder/mergefiles/mergefiles_caller.pl";
my $text2excel ="perl $sbptoolsfolder/text2excel/text2excel.pl";

my $parallel_job ="perl $sbptoolsfolder/parallel-job/parallel-job_caller.pl";

my %commands2program=(
    "bs-fastq"=>$bs_fastq,
	"geo-download"=>$geo_download,
	
    "rnaseq-process"=>$rnaseq_process,
    "rnaseq-merge"=>$rnaseq_merge,
    "rnaseq-de"=>$rnaseq_de,
	"rnaseq-summary"=>$rnaseq_summary,

	"rnaseq-var"=>$rnaseq_var,
	"rnaseq-var-summary"=>$rnaseq_var_summary,
	"rnaseq-motif"=>$rnaseq_motif,
	"rnaseq-motif-summary"=>$rnaseq_motif_summary,
	"gsea-gen"=>$gsea_gen,
	
    "chipseq-process"=>$chipseq_process,
    "chipseq-merge"=>$chipseq_merge,
    "chipseq-de"=>$chipseq_de,
	"chipseq-summary"=>$chipseq_summary,	

    "dnaseq-process"=>$dnaseq_process,

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
