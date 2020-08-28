#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);



########
#Interface
########


my $version="0.6";

#version 0.2a, add r version log
#v0.3 add runmode
#v0.31, solves screen envinroment problem
#v0.4, add server option, updating r script
#v0.5, use R 4.0.2, add comparison argument for multiple comparisons
#v0.51, versioning
#v0.6, major updates planned for R4.0, comparisons of multiple groups. turn off txde

my $usage="

rnaseq-de
version: $version
Usage: sbptools rnaseq-de [parameters]

Description: Differential Expression (DE) tests using DESeq2. This script works for most of the counting based data, e.g. RNA-Seq, ChIP-Seq, ATAC-Seq


Mandatory Parameters:
    --in|-i           Input folder from rnaseq-merge

    --output|-o       Output folder
                         Changes since v0.5. If your output folder in -o DE/,
                         the comparison results will be save in DE/treatment_vs_reference folder

    --config|-c       Configuration file match the samples in the rnaseq-merge folder
                           first column as sample name.

    --tx|-t           Transcriptome
                        Current support Human.B38.Ensembl84, Mouse.B38.Ensembl84

    --formula|-f      Formula for GLM, e.g. ~Group.
                          the last factor of the formula is used for comparison

    #if you have multiple comparisons to perform in a project
    --comparisons     Tab delimited file with first column as treatment groups
                        and second column as reference groups for multiple pairwise comparisons

    #if you only have one comparison
    --treatment       Treatment group name
    --reference       Reference group name

						
Optional Parameters:
    --pmethod         DESeq2 method, default as Wald test [Wald]
    --qmethod         Multiple testing correction method [BH]

    --useallsamples   Use all samples in config file for 
                           gene dispersion calcultation [F]

    --filter          Signal filter [auto]
                         automatically defined signal cutoff as
                           Count >= No. of samples * 5
                         or can define a count number

    --fccutoff        Log2 FC cutoff [1]
    --qcutoff         Corrected P cutoff [0.05]

    --txde            Run DE tests for Tx [F]

    --runmode|-r      Where to run the scripts, local, cluster or none [none]

	#Parallel computing controls	
    --task            Number of tasks to be paralleled. By default 5 tasks. [5]
    --ncpus           No. of cpus for each task [2]
    --mem|-m          Memory usage for each process, e.g. 100mb, 100gb	


";

#    --verbose|-v      Verbose

#add --comparison argument to read comparisons.txt for pairwise comparisons...
#may need to turn off -o, and use comparison name as output folder name to be compatible with gsea-gen


#R parameters
my $rparams="
  -i, --in IN			Expr input file
  -a, --anno ANNO			Annotation file
  -o, --out OUT			Output file
  -f, --formula FORMULA			DESeq formula
  -t, --treat TREAT			treatment name
  -r, --ref REF			reference name
  --fccutoff FCCUTOFF			Log2 FC cutoff [default: 1]
  -q, --qcutoff QCUTOFF			qcutoff [default: 0.05]
  -p, --pmethod PMETHOD			Method used in DESeq2 [default: Wald]
  --qmethod QMETHOD			FDR Method [default: BH]
  --filter FILTER			Count filter for all the exps [default: 10]
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
	
my $inputfolder;	
my $configfile;
my $outputfolder;

my $formula;
my $comparisons;
my $treatment;
my $reference;

my $pmethod="Wald";
my $qmethod="BH";
my $useallsamples="F";
my $filter="auto";

my $txde="F";
my $fccutoff=1;
my $qcutoff=0.05;

my $verbose=1;
my $task=5;
my $ncpus=2;
my $mem;

my $dev=0; #developmental version 


my $tx;
my $runmode="none";

GetOptions(
	"in|i=s" => \$inputfolder,
	"config|c=s" => \$configfile,
	"out|o=s" => \$outputfolder,
	"formula|f=s" => \$formula,
	
	"comparisons=s" => \$comparisons,
	"treatment=s" => \$treatment,
	"reference=s" => \$reference,
	
	"pmethod=s" => \$pmethod,
	"qmethod=s" => \$qmethod,
	"useallsamples=s" => \$useallsamples,
	"filter=s" => \$filter,
	
	"fccutoff=s" => \$fccutoff,
	"qcutoff=s" => \$qcutoff,
	
	"txde=s" => \$txde,
	
	"tx|t=s" => \$tx,	
	"runmode|r=s" => \$runmode,		
	"task=s" => \$task,	
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
	$sbptoolsfolder=get_parent_folder(dirname(abs_path($0)));
}


#sbptools
my $descript="$sbptoolsfolder/rnaseq-de/de_test_caller.R";
my $mergefiles="$sbptoolsfolder/mergefiles/mergefiles_caller.pl";
my $parallel_job="$sbptoolsfolder/parallel-job/parallel-job_caller.pl";


my $r=find_program("/apps/R-4.0.2/bin/R");
my $rscript=find_program("/apps/R-4.0.2/bin/Rscript");

#my $r=find_program("/apps/R-3.4.1/bin/R");
#my $rscript=find_program("/apps/R-3.4.1/bin/Rscript");


#######
#Input/Output
#######

$inputfolder = abs_path($inputfolder);

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


my $logfile="$outputfolder/rnaseq-de_run.log";

my $rlogfile="$outputfolder/rnaseq-de_r_env.log";

my $scriptfile1="$scriptfolder/rnaseq-de_run.sh";

my $scriptlocalrun="$outputfolder/rnaseq-de_local_submission.sh";
my $scriptclusterrun="$outputfolder/rnaseq-de_cluster_submission.sh";



#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();
my $timestamp=build_timestamp($now,"long");


print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";

print LOG "\n";

print STDERR "\nsbptools rnaseq-de $version running ...\n\n" if $verbose;
print LOG "\nsbptools rnaseq-de $version running ...\n\n";



#Report R package version here !!!

open(RLOG,"|$r --no-restore --no-save --slave") || die $!;
select RLOG;
print << "CODE";

rinfo<-c()

Rversion<-getRversion()

rinfo<-rbind(rinfo,c("R",R.Version()\$version.string))
rinfo<-rbind(rinfo,c("Rscript","$rscript"))
rinfo<-rbind(rinfo,c("R library",paste(.libPaths(), collapse=",")))

for (package in c("DESeq2","argparser","ggplot2","EnhancedVolcano")) {
	rinfo<-rbind(rinfo,c(package,packageDescription(package,fields="Version")))
}

colnames(rinfo)<-c("Package","Version")

write.table(file="$rlogfile",rinfo,sep="\t",quote=F,row.names=F)

q()
CODE
close RLOG;


##test tx option
#may need to change for different annotation versions

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


if(defined $tx2ref{$tx}) {
	print STDERR "Starting analysis using $tx.\n\n" if $verbose;
	print LOG "Starting analysis using $tx.\n\n";
}
else {
	print STDERR "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n" if $verbose;
	print LOG "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n";
	exit;
}


#RNA-Seq

#my $genecountmerged="gene.results.merged.count.txt";
#my $txcountmerged="tx.results.merged.count.txt";
#my $genederesult="gene.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt";
#my $txderesult="tx.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt";


my %rnaseq2files=(
	"gene"=> { 
		"count"=> "gene.results.merged.count.txt",
		"selected"=> "gene.results.merged.count.selected.txt",
		"result"=> "gene.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt",
		"resultanno"=> "gene.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.anno.txt",
	},
	"tx"=> { 
		"count"=> "tx.results.merged.count.txt",
		"selected"=> "tx.results.merged.count.selected.txt",
		"result"=> "tx.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt",
		"resultanno"=> "tx.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.anno.txt",
	}
);


#ChIP/ATAC-Seq

my %chip2files=(
	"gene"=> { 
		"count"=> "gene.results.merged.count.txt",
		"selected"=> "gene.results.merged.count.selected.txt",
		"result"=> "gene.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt"
	},
	"tx"=> { 
		"count"=> "tx.results.merged.count.txt",
		"selected"=> "tx.results.merged.count.selected.txt",
		"result"=> "tx.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt"
	}
);


########
#Process
########


#open config files to find fastqs
#----------------
#read design

my %factors;
my @factors_array;

if($formula=~/\~(.+)/) {
	my $cont=$1;
	#a bit loose on factor name
	while($cont=~/([^\+]+)/g) {
		
		#remove trailing white spaces
		my $factor=$1;
		$factor=~s/^\s+//;
		$factor=~s/\s+$//;
		
		$factors{$factor}++;
		push @factors_array,$factor;
	}
}

print STDERR join(",",@factors_array)," factors are identified from -f $formula\n\n" if $verbose;
print LOG join(",",@factors_array)," factors are identified from -f $formula\n\n";


#---------------
#read comparisons

my %comparisons_all;

my @trts;
my @refs;

if(defined $comparisons && length($comparisons)>0) {
	open(IN,$comparisons) || die $!;
	while(<IN>) {
		tr/\r\n//d;
		#no header
		
		my @array=split/\t/;
		
		push @trts,$array[0];
		push @refs,$array[1];
	}
}
else {
	push @trts,$treatment;
	push @refs,$reference;
}


		
#----------------
open(S1,">$scriptfile1") || die "Error writing $scriptfile1. $!";

for(my $compnum=0;$compnum<@trts;$compnum++) {

	my $trt=$trts[$compnum];
	my $ref=$refs[$compnum];


	my $outputfolder_de="$outputfolder/$trt\_vs\_$ref";

	if(!-e $outputfolder_de) {
		mkdir($outputfolder_de);
	}

	
	my $newconfigfile="$outputfolder_de/rnaseq-de_config.txt";
	
	#read config file
	my %sample2fastq;
	my %sample2indexname;
	my %configattrs;
	my %attr2name;

	my @configsamples;

	my $fileline=0;
	my @attrselcols;
	my @sampleselrows;

	open(IN,$configfile) || die "Error reading $configfile. $!";
	open(OUT,">$newconfigfile") || die "Error reading $newconfigfile. $!";

	while(<IN>) {
		tr/\r\n//d;
		my @array=split/\t/;
		if($fileline==0) {
			for(my $num=0;$num<@array;$num++) {
				#case insensitive match to key words
				$configattrs{uc $array[$num]}=$num;
			}
			
			foreach my $factor (@factors_array) {
				if(defined $configattrs{uc $factor}) {
					print STDERR "$factor is identified at the ",$configattrs{uc $factor}+1,"(th) column of $configfile.\n" if $verbose;
					print LOG "$factor is identified at the ",$configattrs{uc $factor}+1,"(th) column of $configfile.\n";
					push @attrselcols,$configattrs{uc $factor};
				}
				else {
					print STDERR "ERROR:$factor is not defined in $configfile.\n";
					print LOG "ERROR:$factor is not defined in $configfile.\n";
					exit;
				}
			}
			
			#print out new config for DE
			print OUT "Sample\t",join("\t",@factors_array),"\n";
		
		}
		else {
			
			push @configsamples,$array[0];
			
			#print out config file for DE
			if($useallsamples eq "T") {
				print OUT join("\t",@array[0,@attrselcols]),"\n";
			}
			else {
				#only print out used samples
				if($array[$configattrs{uc $factors_array[$#factors_array]}] eq $trt || $array[$configattrs{uc $factors_array[$#factors_array]}] eq $ref) {
					print OUT join("\t",@array[0,@attrselcols]),"\n";
					push @sampleselrows,$fileline+1; #sample row number in config is the same with sample col number in count file
				}
			}
			
			foreach my $factor (@factors_array) {
				$attr2name{$factor}{$array[$configattrs{uc $factor}]}++;
			}
		}
		$fileline++;
	}

	close IN;
	close OUT;

	#----------------
	#test treat and ref
	if(defined $attr2name{$factors_array[$#factors_array]}{$trt}) {
		print STDERR $attr2name{$factors_array[$#factors_array]}{$trt}," samples identified for $trt in $configfile.\n" if $verbose;
		print LOG $attr2name{$factors_array[$#factors_array]}{$trt}," samples identified for $trt in $configfile.\n";
	}
	else {
		print STDERR "ERROR:$trt not defined in $configfile.\n";
		print LOG "ERROR:$trt not defined in $configfile.\n";	
		exit;
	}

	if(defined $attr2name{$factors_array[$#factors_array]}{$ref}) {
		print STDERR $attr2name{$factors_array[$#factors_array]}{$ref}," samples identified for $ref in $configfile.\n" if $verbose;
		print LOG $attr2name{$factors_array[$#factors_array]}{$ref}," samples identified for $ref in $configfile.\n";
	}
	else {
		print STDERR "ERROR:$ref not defined in $configfile.\n";
		print LOG "ERROR:$ref not defined in $configfile.\n";	
		exit;
	}


	#----------------
	#check input folder

	if(-e "$inputfolder/".$rnaseq2files{"gene"}{"count"}) {
		open(IN,"$inputfolder/".$rnaseq2files{"gene"}{"count"}) || die $!;
		while(<IN>) {
			tr/\r\n//d;
			my @array=split/\t/;
			my @mergesamples=@array[1..$#array];
			
			if(join(",",@configsamples) ne join(",",@mergesamples)) {
				print STDERR "ERROR:Sample order different.\n";
				print STDERR "ERROR:Configure file $configfile sample order:",join(",",@configsamples),"\n";
				print STDERR "ERROR:Merged file $inputfolder/",$rnaseq2files{"gene"}{"count"}," sample order:",join(",",@mergesamples),"\n";
				
				print LOG "ERROR:Sample order different.\n";
				print LOG "ERROR:Configure file $configfile sample order:",join(",",@configsamples),"\n";
				print LOG "ERROR:Merged file $inputfolder/",$rnaseq2files{"gene"}{"count"}," sample order:",join(",",@mergesamples),"\n";

				exit;
			}
			else {
				print STDERR "Sample order matched.\n" if $verbose;
				print STDERR "Configure file $configfile sample order:",join(",",@configsamples),"\n" if $verbose;
				print STDERR "Merged file $inputfolder/",$rnaseq2files{"gene"}{"count"}," sample order:",join(",",@mergesamples),"\n\n" if $verbose;
				
				print LOG "Sample order matched.\n";
				print LOG "Configure file $configfile sample order:",join(",",@configsamples),"\n";
				print LOG "Merged file $inputfolder/",$rnaseq2files{"gene"}{"count"}," sample order:",join(",",@mergesamples),"\n\n";			
			}
			last;
		}
		close IN;
	}
	else {
		print STDERR "ERROR:$inputfolder/",$rnaseq2files{"gene"}{"count"}," doesn't exist. You need to provide a rnaseq-merge folder.\n";
		print LOG "ERROR:$inputfolder/",$rnaseq2files{"gene"}{"count"}," doesn't exist. You need to provide a rnaseq-merge folder.\n";
		exit;
	}

	########
	#Filter samples
	########

	if($useallsamples eq "T") {
		#copy counting files to de folder
		print STDERR "--useallsamples T defined. All samples are used for DE test.\n\n" if $verbose;
		print LOG "--useallsamples T defined. All samples are used for DE test.\n\n";
		
		system("cp $inputfolder/".$rnaseq2files{"gene"}{"count"}." $outputfolder_de/".$rnaseq2files{"gene"}{"selected"});
		
		if($txde eq "T") {
			system("cp $inputfolder/".$rnaseq2files{"tx"}{"count"}." $outputfolder_de/".$rnaseq2files{"tx"}{"selected"});
		}
	}
	else {
		#copy counting files to de folder, using selected samples
		print STDERR "--useallsamples F defined. Only selected samples are used for DE test.\n" if $verbose;
		print STDERR "Sample columns ".join(",",@sampleselrows)." are used.\n\n" if $verbose;
		
		print LOG "--useallsamples F defined. Only selected samples are used for DE test.\n";
		print LOG "Sample columns ".join(",",@sampleselrows)." are used.\n\n";

		system("cut -f 1,".join(",",@sampleselrows)." $inputfolder/".$rnaseq2files{"gene"}{"count"}." > $outputfolder_de/".$rnaseq2files{"gene"}{"selected"});
		
		if($txde eq "T") {
			system("cut -f 1,".join(",",@sampleselrows)." $inputfolder/".$rnaseq2files{"tx"}{"count"}." > $outputfolder_de/".$rnaseq2files{"tx"}{"selected"});
		}
	}



	#Print out script
	
	#foreach my $sample (sort keys %sample2workflow) {
	#	print S1 $sample2workflow{$sample},"\n";
	#}

	#Gene #changed after v0.6
	print S1 "$rscript $descript -i $outputfolder_de/",$rnaseq2files{"gene"}{"selected"}," -c $newconfigfile -o $outputfolder_de/",$rnaseq2files{"gene"}{"result"}," -f \"$formula\" -t $trt -r $ref --fccutoff $fccutoff --qcutoff $qcutoff --qmethod $qmethod --pmethod $pmethod --filter $filter -a ",$tx2ref{$tx}{"geneanno"}," > $outputfolder_de/gene_de_test_run.log 2>&1;"; #need to check here #add r script running log

	#Gene anno
	print S1 "$mergefiles -m $outputfolder_de/",$rnaseq2files{"gene"}{"result"}," -i ".$tx2ref{$tx}{"geneanno"}." -o $outputfolder_de/",$rnaseq2files{"gene"}{"resultanno"},";\n";


	if($txde eq "T") {
		#Tx
		print S1 "$rscript $descript -i $outputfolder_de/",$rnaseq2files{"tx"}{"selected"}," -c $newconfigfile -o $outputfolder_de/",$rnaseq2files{"tx"}{"result"}," -f \"$formula\" -t $trt -r $ref --fccutoff $fccutoff --qcutoff $qcutoff --qmethod $qmethod --pmethod $pmethod --filter $filter"," > $outputfolder_de/tx_de_test_run.log 2>&1;";

		#tx anno
		print S1 "$mergefiles -m $outputfolder_de/",$rnaseq2files{"tx"}{"result"}," -i ".$tx2ref{$tx}{"txanno"}." -o $outputfolder_de/",$rnaseq2files{"tx"}{"resultanno"},";\n";
	}

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
my $jobname="rnaseq-de-$timestamp";

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

my $clustercommand="perl $parallel_job -i ".join(",", @scripts_all)." -o $scriptfolder -n ".join(",",@script_names)." --tandem -t $task --ncpus $ncpus --env -r ";

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

sub get_parent_folder {
	my $dir=shift @_;
	
	if($dir=~/^(.+\/)[^\/]+\/?/) {
		return $1;
	}
}



