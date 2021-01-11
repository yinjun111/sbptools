#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);

########
#Interface
########


my $version="0.41";

#version 0.2, add de summary script
#v0.3, add de summary by peak calling
#v0.31, solves screen envinroment problem
#v0.4, Firefly compatible. Versioning. R 4.0
#v0.41, skip input for ChIP-Seq. New inputfile for DE

my $usage="

chipseq-de
version: $version
Usage: sbptools chipseq-de [parameters]

Description: Differential Expression (DE) tests using DESeq2 for ChIP-Seq and ATAC-Seq


Mandatory Parameters:
    --in|-i           Input folder from chipseq-merge
    --output|-o       Output folder
    --config|-c       Configuration file match the samples in the chipseq-merge folder
                           first column as sample name.

    --formula|-f      Formula for GLM, e.g. ~Group.
                          the last factor of the formula is used for comparison
    --treatment       Treatment group name
    --reference       Reference group name

    --tx|-t           Transcriptome
                        Current support Human.B38.Ensembl84, Mouse.B38.Ensembl84
						
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

    --runmode|-r      Where to run the scripts, cluster, local, or none [none]

    Parallel computating parameters
    --task            Number of tasks to be paralleled. By default 4 tasks for local mode, 8 tasks for cluster mode.
    --ncpus           No. of cpus for each task [4]
    --mem|-m          Memory usage for each process, e.g. 100mb, 100gb [40gb]

	
";


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
my $treatment;
my $reference;

my $pmethod="Wald";
my $qmethod="BH";
my $useallsamples="F";
my $filter="auto";

my $fccutoff=1;
my $qcutoff=0.05;

my $verbose=1;
my $tx;
my $runmode="none";

my $task;
my $ncpus=4;
my $mem="40gb";
my $dev=0; #developmental version


GetOptions(
	"in|i=s" => \$inputfolder,
	"config|c=s" => \$configfile,
	"out|o=s" => \$outputfolder,
	"formula|f=s" => \$formula,
	
	"treatment=s" => \$treatment,
	"reference=s" => \$reference,
	
	"pmethod=s" => \$pmethod,
	"qmethod=s" => \$qmethod,
	"useallsamples=s" => \$useallsamples,
	"filter=s" => \$filter,
	
	"fccutoff=s" => \$fccutoff,
	"qcutoff=s" => \$qcutoff,
	
	"tx|t=s" => \$tx,	
	"runmode|r=s" => \$runmode,		

	"task=s" => \$task,
	"ncpus=s" => \$ncpus,
	"mem=s" => \$mem,

	"dev" => \$dev,			

	"verbose|v" => \$verbose,
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
else {
	#the tools called will be within the same folder of the script
	$sbptoolsfolder=get_parent_folder(abs_path(dirname($0)));
}


my $parallel_job="$sbptoolsfolder/parallel-job/parallel-job_caller.pl";
my $mergefiles="$sbptoolsfolder/mergefiles/mergefiles_caller.pl";
my $processmergebed="$sbptoolsfolder/chipseq-merge/process_mergebed.pl";
my $selectbed="$sbptoolsfolder/chipseq-merge/select_mergebed.pl";

my $descript="$sbptoolsfolder/rnaseq-de/de_test_caller.R";
my $reformatpeakcount="$sbptoolsfolder/chipseq-de/reformat_peak_count.pl";
my $desummary="$sbptoolsfolder/chipseq-de/summarize_dm_peaks.pl";
my $desummarybycalling="$sbptoolsfolder/chipseq-de/summarize_dm_peaks_bycalling.pl";
my $convertdetopos="$sbptoolsfolder/chipseq-de/convert_chipseq_de_to_pos.pl";


#homer
my $homer="/apps/homer/bin/";
my $findmotifsgenome="$homer/findMotifsGenome.pl";


my $r=find_program("/apps/R-4.0.2/bin/R");
my $rscript=find_program("/apps/R-4.0.2/bin/Rscript");



########
#I/O
########


$inputfolder = abs_path($inputfolder);

if(!-e $outputfolder) {
	mkdir($outputfolder);
}

$outputfolder = abs_path($outputfolder);


my $scriptfolder="$outputfolder/scripts";

if(!-e $scriptfolder) {
	mkdir($scriptfolder);
}


my $logfile="$outputfolder/chipseq-de_run.log";
my $newconfigfile="$outputfolder/chipseq-de_config.txt";

my $scriptfile1="$scriptfolder/chipseq-de_run.sh";

my $scriptlocalrun="$outputfolder/chipseq-de_local_submission.sh";
my $scriptclusterrun="$outputfolder/chipseq-de_cluster_submission.sh";


#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();
my $timestamp=build_timestamp($now,"long");


print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";

#Report R package version here !!!



print LOG "\n";

print STDERR "\nsbptools chipseq-de $version running ...\n\n" if $verbose;
print LOG "\nsbptools chipseq-de $version running ...\n\n";


#test tx option

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
}

my $genomeversion;

if($tx=~/Human.B38/) {
	$genomeversion="hg38";
}
elsif($tx=~/Mouse.B38/) {
	$genomeversion="mm10";
}

#RNA-Seq

#my $genecountmerged="gene.results.merged.count.txt";
#my $txcountmerged="tx.results.merged.count.txt";
#my $genederesult="gene.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt";
#my $txderesult="tx.results.merged.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt";


my %chipseq2files=(
	"all"=> { 
		"raw"=> "all.reprod.peak.merged.raw.count.txt",
		"count"=> "all.reprod.peak.merged.raw.countonly.txt",
		"anno"=> "all.reprod.peak.merged.raw.anno.txt",	
		"selected"=> "all.reprod.peak.merged.raw.countonly.selected.txt",
		"result"=> "all.reprod.peak.merged.raw.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt",

		"depos"=> "all.reprod.peak.merged.raw.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.pos",
		"deposup"=> "all.reprod.peak.merged.raw.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff\_up.pos",
		"deposdown"=> "all.reprod.peak.merged.raw.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff\_down.pos",

		"resultanno"=> "all.reprod.peak.merged.raw.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.anno.txt",
		"summary"=> "all.reprod.peak.merged.raw.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.summary.txt",
		"summaryanno"=> "all.reprod.peak.merged.raw.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.summary.anno.txt",

		"resultbycalling"=> "all.reprod.peak.merged.summary.txt",
		"resultbycallinganno"=> "all.reprod.peak.merged.dm.bycalling.anno.txt",
		"summarybycalling"=> "all.reprod.peak.merged.dm.bycalling.summary.txt",
		"summarybycallinganno"=> "all.reprod.peak.merged.dm.bycalling.summary.anno.txt",		
	},
	"1000u0d_longest"=> { 
		"raw"=> "1000u0d_longest_promoter.merged.raw.count.txt",
		"count"=> "1000u0d_longest_promoter.merged.raw.countonly.txt",
		"anno"=> "1000u0d_longest_promoter.merged.raw.anno.txt",	
		"selected"=> "1000u0d_longest_promoter.merged.raw.countonly.selected.txt",
		"result"=> "1000u0d_longest_promoter.merged.raw.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt",
		"resultanno"=> "1000u0d_longest_promoter.merged.raw.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.anno.txt",
	},
	"1000u0d_all"=> { 
		"raw"=> "1000u0d_all_promoter.merged.raw.count.txt",
		"count"=> "1000u0d_all_promoter.merged.raw.countonly.txt",
		"anno"=> "1000u0d_all_promoter.merged.raw.anno.txt",	
		"selected"=> "1000u0d_all_promoter.merged.raw.countonly.selected.txt",
		"result"=> "1000u0d_all_promoter.merged.raw.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.txt",
		"resultanno"=> "1000u0d_all_promoter.merged.raw.count.DESeq2.$pmethod.FC$fccutoff.$qmethod.P$qcutoff.anno.txt",
	},	

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
		
#----------------
#read config file, generate new config for samples subset

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
				print STDERR "$factor is identified at the ",$configattrs{uc $factor},"(th) column of $configfile.\n" if $verbose;
				print LOG "$factor is identified at the ",$configattrs{uc $factor},"(th) column of $configfile.\n";
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
		
		#skip input rows		

		#find input sample
		#there should not be any INPUT in the config file
		if(defined $array[$configattrs{"INPUT"}] && length($array[$configattrs{"INPUT"}])>0 && $array[$configattrs{"INPUT"}]=~/Y/i) {
			next;
		}
		else {			
			push @configsamples,$array[0];
			
			#print out config file for DE
			if($useallsamples eq "T") {
				print OUT join("\t",@array[0,@attrselcols]),"\n";
			}
			else {
				#only print out used samples
				if($array[$configattrs{uc $factors_array[$#factors_array]}] eq $treatment || $array[$configattrs{uc $factors_array[$#factors_array]}] eq $reference) {
					print OUT join("\t",@array[0,@attrselcols]),"\n";
					push @sampleselrows,$fileline+1; #sample row number in config is the same with sample col number in count file
				}
			}
			
			foreach my $factor (@factors_array) {
				$attr2name{$factor}{$array[$configattrs{uc $factor}]}++;
			}
		}
		
	}
	$fileline++;
}

close IN;
close OUT;

#----------------
#test treat and ref
if(defined $attr2name{$factors_array[$#factors_array]}{$treatment}) {
	print STDERR $attr2name{$factors_array[$#factors_array]}{$treatment}," samples identified for $treatment in $configfile.\n" if $verbose;
	print LOG $attr2name{$factors_array[$#factors_array]}{$treatment}," samples identified for $treatment in $configfile.\n";
}
else {
	print STDERR "ERROR:$treatment not defined in $configfile.\n";
	print LOG "ERROR:$treatment not defined in $configfile.\n";	
	exit;
}

if(defined $attr2name{$factors_array[$#factors_array]}{$reference}) {
	print STDERR $attr2name{$factors_array[$#factors_array]}{$reference}," samples identified for $reference in $configfile.\n" if $verbose;
	print LOG $attr2name{$factors_array[$#factors_array]}{$reference}," samples identified for $reference in $configfile.\n";
}
else {
	print STDERR "ERROR:$reference not defined in $configfile.\n";
	print LOG "ERROR:$reference not defined in $configfile.\n";	
	exit;
}


#double check order of treatment and ref




#----------------
#check input folder

#count and newconfig are in the same sample order after reformatpeakcount

#convert files
foreach my $attr ("all","1000u0d_longest","1000u0d_all") {
	print STDERR "Converting $inputfolder/",$chipseq2files{$attr}{"raw"}," for DE test.\n";
	print LOG "Converting $inputfolder/",$chipseq2files{$attr}{"raw"}," for DE test.\n";
	
	system("$reformatpeakcount -i $inputfolder/".$chipseq2files{$attr}{"raw"}." -o $outputfolder/".$chipseq2files{$attr}{"count"}." -a $outputfolder/".$chipseq2files{$attr}{"anno"}." -c $newconfigfile");
	print LOG "$reformatpeakcount -i $inputfolder/".$chipseq2files{$attr}{"raw"}." -o $outputfolder/".$chipseq2files{$attr}{"count"}." -a $outputfolder/".$chipseq2files{$attr}{"anno"}." -c $newconfigfile\n";
}

########
#Filter samples
########

#this step is unneccesary. 
#$outputfolder/".$chipseq2files{$attr}{"count"} should be the file needed for DE


#if($useallsamples eq "T") {
#	#copy counting files to de folder
#	print STDERR "--useallsamples T defined. All samples are used for DE test.\n\n" if $verbose;
#	print LOG "--useallsamples T defined. All samples are used for DE test.\n\n";
	
#	foreach my $attr ("all","1000u0d_longest","1000u0d_all") {
#		system("cp $outputfolder/".$chipseq2files{$attr}{"count"}." $outputfolder/".$chipseq2files{$attr}{"selected"});
#		print LOG "cp $outputfolder/".$chipseq2files{$attr}{"count"}." $outputfolder/".$chipseq2files{$attr}{"selected"},"\n";
#	}
#}
#else {
#	#copy counting files to de folder, using selected samples
#	print STDERR "--useallsamples F defined. Only selected samples are used for DE test.\n" if $verbose;
#	print STDERR "Sample columns ".join(",",@sampleselrows)." are used.\n\n" if $verbose;
	
#	print LOG "--useallsamples F defined. Only selected samples are used for DE test.\n";
#	print LOG "Sample columns ".join(",",@sampleselrows)." are used.\n\n";

#	foreach my $attr ("all","1000u0d_longest","1000u0d_all") {
#		system("cut -f 1,".join(",",@sampleselrows)." $outputfolder/".$chipseq2files{$attr}{"count"}." > $outputfolder/".$chipseq2files{$attr}{"selected"});
#		print LOG "cut -f 1,".join(",",@sampleselrows)." $outputfolder/".$chipseq2files{$attr}{"count"}." > $outputfolder/".$chipseq2files{$attr}{"selected"},"\n";
#	}
#}

########
#Print out commands, for local and server run
########


open(S1,">$scriptfile1") || die "Error writing $scriptfile1. $!";

#foreach my $sample (sort keys %sample2workflow) {
#	print S1 $sample2workflow{$sample},"\n";
#}


foreach my $attr ("all","1000u0d_longest","1000u0d_all") {
	
	#DE by Signal ...
	print S1 "$rscript $descript -i $outputfolder/",$chipseq2files{$attr}{"count"}," -c $newconfigfile -o $outputfolder/",$chipseq2files{$attr}{"result"}," -f \"$formula\" -t $treatment -r $reference --fccutoff $fccutoff --qcutoff $qcutoff --qmethod $qmethod --pmethod $pmethod --filter $filter;";


	#summary and annotation 
	
	if($attr eq "all") {
	
		#DE by Signal, annotation
		#Gene anno
		print S1 "$mergefiles -m $outputfolder/",$chipseq2files{$attr}{"result"}," -i $outputfolder/",$chipseq2files{$attr}{"anno"}," -o $outputfolder/",$chipseq2files{$attr}{"resultanno"},";";
		
		#summarize data
		print S1 "$desummary -i $outputfolder/",$chipseq2files{$attr}{"resultanno"}," --tx $tx --out $outputfolder/",$chipseq2files{$attr}{"summary"},";";
		
		#anno
		print S1 "$mergefiles -m $outputfolder/",$chipseq2files{$attr}{"summary"}," -i ",$tx2ref{$tx}{"geneanno"}," -o $outputfolder/",$chipseq2files{$attr}{"summaryanno"},";";
		
		print S1 "\n";
		
		
		#DE by peak calling ...
		print S1 "$desummarybycalling -i $inputfolder/",$chipseq2files{$attr}{"resultbycalling"}," --treatment $treatment --reference $reference --tx $tx --o1 $outputfolder/",$chipseq2files{$attr}{"resultbycallinganno"}," --o2 $outputfolder/",$chipseq2files{$attr}{"summarybycalling"}," -a $inputfolder/",$chipseq2files{$attr}{"raw"}," >> $logfile 2>&1   ;";		#need to add ref and treat
		
		print S1 "$mergefiles -m $outputfolder/",$chipseq2files{$attr}{"summarybycalling"}," -i ",$tx2ref{$tx}{"geneanno"}," -o $outputfolder/",$chipseq2files{$attr}{"summarybycallinganno"},";";
		
		
		#TFBS for DM by Signal #didn't work ....
		#print S1 "$convertdetopos -i $outputfolder/",$chipseq2files{$attr}{"result"}," -o $outputfolder/",$chipseq2files{$attr}{"depos"},";";
		#run the three processes in bg
		#print S1 "$findmotifsgenome $outputfolder/",$chipseq2files{$attr}{"depos"}," $genomeversion $outputfolder/deboth_tfbs -size given &";
		#print S1 "$findmotifsgenome $outputfolder/",$chipseq2files{$attr}{"deposup"}," $genomeversion $outputfolder/deup_tfbs -size given &";
		#print S1 "$findmotifsgenome $outputfolder/",$chipseq2files{$attr}{"deposdown"}," $genomeversion $outputfolder/dedown_tfbs -size given &";
		
		print S1 "\n";	
	}
	else {
		#txs
		#Gene anno
		print S1 "$mergefiles -m $outputfolder/",$chipseq2files{$attr}{"result"}," -i $outputfolder/",$chipseq2files{$attr}{"anno"},",",$tx2ref{$tx}{"txanno"}," -o $outputfolder/",$chipseq2files{$attr}{"resultanno"},";\n";
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
my $jobname="chipseq-de-$timestamp";

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
	
	#use defined program as default, otherwise search for this program in PATH
	
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





