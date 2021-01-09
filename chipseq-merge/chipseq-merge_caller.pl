#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);


########
#Interface
########


my $version="0.5";

#v0.2 now skips input samples for merging
#v0.3, removed -v, add -r implementation for local
#v0.31, solves screen envinroment problem
#v0.4, firefly support
#v0.5, add expr-qc

my $usage="

chipseq-merge
version: $version
Usage: sbptools chipseq-merge [parameters]

Description: Merge chipseq-process folder to get summarized QC, counting for promoter and merged peak, and so on

Parameters:

    --in|-i           Input chipseq-process folder(s)
	
    --config|-c       Configuration file #at this stage, configuration only, to make merged peaks
    --group|-g        Column name for sample groups to merge peaks

    --groupsubset|--gss     Subset of groups selected for the merged folder, optional
    --samplesubset|--sss    Subset of samples selected for the merged folder, optional

    --output|-o       Output folder


    --tx|-t           Transcriptome
                        Current support Human.B38.Ensembl84, Mouse.B38.Ensembl84
    --anno|-a         Add annotation

    --runmode|-r      Where to run the scripts, cluster, local, or none [none]

    Parallel computating parameters
    --task            Number of tasks to be paralleled. By default 4 tasks for local mode, 8 tasks for cluster mode.
    --ncpus           No. of cpus for each task [4]
    --mem|-m          Memory usage for each process, e.g. 100mb, 100gb [40gb]


";

#    --verbose|-v      Verbose, use -v 0 to turn off verbose [1]

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

my $group;
my $groupsubset;
my $samplesubset;

my $configfile;
my $outputfolder;
my $verbose=1;
my $tx;
my $runmode="none";
#my $jobs=5;

my $task;
my $ncpus=4;
my $mem="40gb";
my $dev=0; #developmental version

GetOptions(
	"in|i=s" => \$inputfolders,
	#"samples|s=s" => \$samples,	
	"config|c=s" => \$configfile,
	"group|g=s" => \$group,
	"groupsubset|gss=s" => \$groupsubset,
	"samplesubset|sss=s" => \$samplesubset,		
	"output|o=s" => \$outputfolder,
	"tx|t=s" => \$tx,
	#"jobs|j=s" => \$jobs,
	
	"task=s" => \$task,
	"ncpus=s" => \$ncpus,
	"mem=s" => \$mem,

	"dev" => \$dev,			

	"runmode|r=s" => \$runmode,		
	"verbose|v=s" => \$verbose,
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
my $expr_qc="$sbptoolsfolder/expr-qc/expr-qc.R";
my $reformatpeakcount="$sbptoolsfolder/chipseq-de/reformat_peak_count.pl";

my $multiqc=find_program("/apps/python-3.5.2/bin/multiqc");
my $cutadapt=find_program("/apps/python-3.5.2/bin/cutadapt");
my $fastqc=find_program("/apps/FastQC/fastqc");
my $rsem=find_program("/apps/RSEM-1.3.1/rsem-calculate-expression");
my $star=find_program("/apps/STAR-master/bin/Linux_x86_64/STAR");
my $bamcoverage=find_program("/apps/python-3.5.2/bin/bamCoverage");
my $samtools=find_program("/apps/samtools-1.3.1/bin/samtools");
my $mergebed=find_program("/apps/bedtools2-2.26.0/bin/bedtools")." merge";
my $bedtobigbed=find_program("/apps/ucsc/bedToBigBed");

my $r=find_program("/apps/R-4.0.2/bin/R");
my $rscript=find_program("/apps/R-4.0.2/bin/Rscript");

#homer
my $homer="/apps/homer/bin/";
#my $homer="/home/jyin/Programs/Homer/bin";
my $maketagdirectory="$homer/makeTagDirectory";
my $annotatepeaks="$homer/annotatePeaks.pl";
my $findpeaks="$homer/findPeaks";
my $makeucscfile="$homer/makeUCSCfile";
my $pos2bed="$homer/pos2bed.pl";
my $bed2pos="$homer/bed2pos.pl";


#########

#default ouputs

my $promotermergedcount_norm="promoter.merged.norm.count.txt"; #count ,annotated count, 
my $promotermergedcount_raw="promoter.merged.raw.count.txt"; #count ,annotated count, 
#my $samplemergedcount="promoter.merged.count.txt"; #peaks for samples from the same group are merged 
my $allmergedcount_norm="all.reprod.peak.merged.norm.count.txt"; #peaks for all samples are merged, based on group merged
my $allmergedcount_raw="all.reprod.peak.merged.raw.count.txt"; #peaks for all samples are merged, based on group merged

my $allmergedcount_norm_reformated="all.reprod.peak.merged.norm.count.reformated.txt"; #peaks for all samples are merged, based on group merged
my $allmergedcount_norm_reformated_anno="all.reprod.peak.merged.norm.count.reformated.anno.txt";

my $promotermergedcount_norm_reformated="promoter.merged.norm.count.reformated.txt"; #count ,annotated count, 
my $promotermergedcount_norm_reformated_anno="promoter.merged.norm.count.reformated.anno.txt"; 


my $allmergedbed="all.reprod.peak.merged.bed";
my $allmergedbed_sorted="all.reprod.peak.merged_sorted.bed";
my $allmergedpos="all.reprod.peak.merged.pos";
my $allmergedbb="all.reprod.peak.merged_sorted.bb";


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

my $logfile="$outputfolder/chipseq-merge_run.log";

my $scriptfile1="$scriptfolder/chipseq-merge_run1.sh";
my $scriptfile2="$scriptfolder/chipseq-merge_run2.sh";


my $scriptlocalrun="$outputfolder/chipseq-merge_local_submission.sh";
my $scriptclusterrun="$outputfolder/chipseq-merge_cluster_submission.sh";


#new config
my $newconfigfile="$outputfolder/chipseq-merge_config.txt";


#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();
my $timestamp=build_timestamp($now,"long");

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";
print LOG "$multiqc version:", getsysoutput("$multiqc --version"),"\n";
print LOG "\n";


print STDERR "\nsbptools chipseq-merge version $version running ...\n\n" if $verbose;
print LOG "\nsbptools chipseq-merge version $version running ...\n\n";

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



#to genome version
my %tx2gv=(
	"Human.B38.Ensembl84"=>"hg38",
	"Mouse.B38.Ensembl84"=>"mm10",
);

my %tx2promoter=(
	"Human.B38.Ensembl84"=> {
		"1000u0d_longest"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_promoter_1kup0down_longesttxs.pos",
		"1000u0d_all"=> "/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_promoter_1kup0down_alltxs.pos"
	},
	"Mouse.B38.Ensembl84"=> {
		"1000u0d_longest"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_promoter_1kup0down_longesttxs.pos",
		"1000u0d_all"=> "/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_promoter_1kup0down_alltxs.pos"
	}
);

my %tx2gtf=(
	"Human.B38.Ensembl84"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc.gtf",
	"Mouse.B38.Ensembl84"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc.gtf",
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


#open config files to find samples

my @samples_array;
my %samples_hash; #check no dup samples
my %groups_hash;

my %samples_hash_sel; #check no dup samples
my %groups_hash_sel;


my @inputsamples_array;
my $groupcol;
my %group2samples;
my %configattrs;

if(defined $samplesubset && length($samplesubset)>0 && defined $groupsubset && length($groupsubset)>0) {
	print STDERR "ERROR either --samplesubset or --groupsubset can be used, not both.\n\n";
	print LOG "ERROR either --samplesubset or --groupsubset can be used, not both.\n\n";
	exit;
}


#selected samples
if(defined $samplesubset && length($samplesubset)>0) {
	foreach my $sample (split(",",$samplesubset)) {
		#keep the order of samples
		push @samples_array, $sample;
		unless (defined $samples_hash_sel{$sample}) {
			$samples_hash_sel{$sample}++
		}
		else {
			print STDERR "ERROR:$sample is defined multiple times in $samples.\n";
			print LOG "ERROR:$sample is defined multiple times in $samples.\n";
			exit;
		}
	}
}

#selected groups
if(defined $groupsubset && length($groupsubset)>0) {
	foreach my $group (split(",",$groupsubset)) {
		#keep the order of samples
		#push @samples_array, $sample;
		unless (defined $groups_hash_sel{$group}) {
			$groups_hash_sel{$group}++
		}
		else {
			print STDERR "ERROR:$group is defined multiple times in $groupsubset.\n";
			print LOG "ERROR:$group is defined multiple times in $groupsubset.\n";
			exit;
		}
	}
}


open(OUT,">$newconfigfile") || die "Error writing $newconfigfile. $!";  #need to be implemented!!!!!

#mandatory, a config file
if(defined $configfile && length($configfile)>0) {
	#another way of getting samples
	
	open(IN,$configfile) || die "Error reading $configfile. $!";
	#first column should be sample, with a header

	my $fileline=0;
	while(<IN>) {
		tr/\r\n//d;
		my @array=split/\t/;
		
		if($fileline==0) {
			#record title into hash
			for(my $num=0;$num<@array;$num++) {
				#case sensitive match to key words, Sample, FASTQ1, FASTQ2, Index (for foldernames)
				$configattrs{uc $array[$num]}=$num;
			}
			
			if(defined $configattrs{uc $group}) {
				$groupcol=$configattrs{uc $group};
				print STDERR "--group $group is identifed at column ",$groupcol+1," in $configfile.\n\n" if $verbose;
				print LOG "--group $group is identifed at column ",$groupcol+1," in $configfile.\n\n";
			}
			else {
				print STDERR "ERROR:--group $group is not found in $configfile.\n\n" if $verbose;
				print LOG "ERROR:--group $group is not found in $configfile.\n\n";
				exit;
			}
			
			if(defined $configattrs{"INPUT"}) {
				print STDERR "Input samples is defined at column ",$configattrs{"INPUT"}+1," in $configfile.\n\n" if $verbose;
				print LOG "Input samples is defined at column ",$configattrs{"INPUT"}+1," in $configfile.\n\n";
			}

			print OUT $_,"\n";
		}
		else {
			my $sample=$array[0];
			
			unless(defined $configattrs{"INPUT"} && $array[$configattrs{"INPUT"}]=~/Y/i ) {
				#don't need input sample in the count file 

				#by sample selection
				if(defined $samplesubset && length($samplesubset)>0) {
					if(defined $samples_hash_sel{$sample}) {
						#push @samples_array, $sample;
					
						unless (defined $samples_hash{$sample}) {
							$samples_hash{$sample}++;
							$group2samples{$array[$groupcol]}{$sample}++;
							print OUT $_,"\n";							
						}
						else {
							print STDERR "ERROR:$sample is defined multiple times in $configfile.\n";
							print LOG "ERROR:$sample is defined multiple times in $configfile.\n";
							exit;
						}
					}
				}
				
				#By group selection 

				if(defined $groupsubset && length($groupsubset)>0) {
					
					if(defined $groups_hash_sel{$array[$groupcol]}) {
						#push @samples_array, $sample;
					
						unless (defined $samples_hash{$sample}) {
							push @samples_array, $sample;
							$samples_hash{$sample}++;
							$group2samples{$array[$groupcol]}{$sample}++;
							print OUT $_,"\n";
						}
						else {
							print STDERR "ERROR:$sample is defined multiple times in $configfile.\n";
							print LOG "ERROR:$sample is defined multiple times in $configfile.\n";
							exit;
						}
					}
				}

				#No subset selection
				unless((defined $samplesubset && length($samplesubset)>0) || (defined $groupsubset && length($groupsubset)>0)) {
					push @samples_array, $sample;
					$samples_hash{$sample}++;
					$group2samples{$array[$groupcol]}{$sample}++;
					print OUT $_,"\n";
				}
				
			}
			else {
				push @inputsamples_array, $sample;
			
			}
			
		}
		
		$fileline++;
	}
	close IN;
}

close OUT;


if(defined $samplesubset && length($samplesubset)>0 && keys %samples_hash != keys %samples_hash_sel) {
	print STDERR "ERROR: samples selected by --samplesubset ",join(",",sort keys %samples_hash_sel)," is different from --config $configfile ",join(",",sort keys %samples_hash),"\n";
	print LOG "ERROR: samples selected by --samplesubset ",join(",",sort keys %samples_hash_sel)," is different from --config $configfile ",join(",",sort keys %samples_hash),"\n";
	exit;
}

if(defined $groupsubset && length($groupsubset)>0 && keys %group2samples != keys %groups_hash_sel) {
	print STDERR "ERROR: groups selected by --groupsubset ",join(",",sort keys %groups_hash_sel)," is different from --config $configfile ",join(",",sort keys %group2samples),"\n";
	print LOG "ERROR: groups selected by --groupsubset ",join(",",sort keys %groups_hash_sel)," is different from --config $configfile ",join(",",sort keys %group2samples),"\n";
	exit;
}


print STDERR scalar(@samples_array)," samples identified: ",join(",",@samples_array),"\n\n" if $verbose;
print LOG scalar(@samples_array)," samples identified: ",join(",",@samples_array),"\n\n";

print STDERR scalar(@inputsamples_array)," input samples identified: ",join(",",@inputsamples_array),"\n\n" if $verbose;
print LOG scalar(@inputsamples_array)," input samples identified: ",join(",",@inputsamples_array),"\n\n";

print STDERR scalar(keys %group2samples)," groups of samples identified: ",join(",",sort keys %group2samples),"\n\n" if $verbose;
print LOG scalar(keys %group2samples)," groups of samples identified: ",join(",",sort keys %group2samples),"\n\n";


#----------------
#Find files for merging
#Examples
#JY_315_Isl.genes.results
#JY_315_Isl.isoforms.results

print STDERR "\nReading sample folders.\n" if $verbose;
print LOG "\nReading sample folders.\n";

#my %sample2gene;
#my %sample2tx;
my %sample2folder;

my %sample2tagdir;

foreach my $infolder (split(",",$inputfolders)) {
	#samples in different folders
	#if different folders have the same sample name, use the first found sample
	
	my @infiles=glob("$infolder/*");
	
	foreach my $file (@infiles) {
		if(-d $file) {
			my $samplename=basename($file);
			
			$sample2folder{$samplename}=abs_path($file);
			
			#see whether it is a defined sample
			if(defined $samples_hash{$samplename}) {
				print STDERR $samplename," is found in ",$file,".\n" if $verbose;
				print LOG $samplename," is found in ",$file,".\n" if $verbose;

				my @samplefiles=glob("$file/*");
				
				foreach  my $samplefile (@samplefiles) {
					#homer TagDir
					if($samplefile=~/_TagDir$/) {
						$sample2tagdir{$samplename}=abs_path($samplefile);
					}
				}
			}
		}
	}
}



print STDERR scalar(keys %sample2tagdir)," samples identified with Homer TagDir.\n\n" if $verbose;
print LOG scalar(keys %sample2tagdir)," samples identified with Homer TagDir.\n\n";


if( scalar(@samples_array) != scalar(keys %sample2tagdir) ) {
	print STDERR "ERROR:Not all samples have Homer TagDir.\n\n";
	print LOG "ERROR:Not all samples have Homer TagDir.\n\n";
	exit;
}



########
#Print out commands, for local and server run
########





open(S1,">$scriptfile1") || die "Error writing $scriptfile1. $!";


#----------------
#Multiqc with specific folders


#file for multiqc folder search
open(OUT,">$tempfolder/samplefolders.txt") || die $!;

#chipped and input samples are both used in the multiqc
foreach my $sample (@samples_array,@inputsamples_array) {
	print OUT $sample2folder{$sample},"\n";
}
close OUT;


print S1 "$multiqc -l $tempfolder/samplefolders.txt -o $outputfolder/multiqc;\n";

#------------------


#Get peaks for promoters

#all tag dirs
my $tagdirs=join(" ",map {$sample2tagdir{$_}} sort keys %sample2tagdir);

foreach my $pv (sort keys %{$tx2promoter{$tx}}) {
	#promoter version
	
	print STDERR "Printing commands to annotate promoter based on ",$tx2promoter{$tx}{$pv},"\n\n" if $verbose;
	print LOG "Printing commands to annotate promoter based on ",$tx2promoter{$tx}{$pv},"\n\n";
	
	#command to get counting for promoters
	print S1 "$annotatepeaks ",$tx2promoter{$tx}{$pv}," ",$tx2gv{$tx}," -gtf ",$tx2gtf{$tx}," -d $tagdirs > $outputfolder/$pv\_$promotermergedcount_norm 2> $outputfolder/$pv\_promoter_annotatepeaks_norm_run.log\n";
	print S1 "$annotatepeaks ",$tx2promoter{$tx}{$pv}," ",$tx2gv{$tx}," -gtf ",$tx2gtf{$tx}," -d $tagdirs -raw > $outputfolder/$pv\_$promotermergedcount_raw 2> $outputfolder/$pv\_promoter_annotatepeaks_raw_run.log\n";
}

#Merge called peaks for sample groups
my @groupselectedbeds;

foreach my $group (sort keys %group2samples) {
	#1) merge bed first
	#2) identify reproducible bed
	#3) get counting and annotation
	
	#get peak bed files for each sample
	
	my @peakfiles;
	my @renamed;
	my $cmd;

	my $mergedbed="$outputfolder/$group\_merged.bed";
	my $selectedbed="$outputfolder/$group\_reprod.bed";
	my $selectedbed_sorted="$outputfolder/$group\_reprod_sorted.bed";
	my $selectedbb="$outputfolder/$group\_reprod_sorted.bb";
	
	push @groupselectedbeds,$selectedbed;

	print STDERR "Producing reproducible merged bed for $group, as $selectedbed\n" if $verbose;
	print LOG "Producing reproducible merged bed for $group, as $selectedbed\n";	

	#rename bed
	foreach my $sample (sort keys %{$group2samples{$group}}) {
		my $peakfile=$sample2folder{$sample}."/$sample\_Peaks.bed";
		my $peakfilerenamed="$tempfolder/".basename($peakfile);
		$peakfilerenamed=~s/.bed$/_renamed.bed/;
		
		push @peakfiles, $peakfile;
		push @renamed, $peakfilerenamed;
		
		print S1 "$processmergebed -i $peakfile -o $peakfilerenamed -s $sample;"
		#output is all.reprod.peak.merged.summary.txt
	}
	
	
	
	if(keys %{$group2samples{$group}} >1) {
		#for group with biological replicate
		
		#merge bed 
		print S1 "cat ",join(" ",@renamed)," | sort -k1,1 -k2,2n | $mergebed -c 4 -o collapse > $mergedbed;";
		
		#select reproducible bed
		print S1 "$selectbed $mergedbed $selectedbed;";
	}
	else {
		print STDERR "WARNING: $group doesn't have biologial replicate. Only 1 sample ",join(",",sort keys %{$group2samples{$group}})," identified.\n\n";
		print LOG "WARNING: $group doesn't have biologial replicate. Only 1 sample ",join(",",sort keys %{$group2samples{$group}})," identified.\n\n";
		
		print S1 "cut -f 1-4 ",join(",",@renamed)," > $mergedbed;";
		
		print S1 "$selectbed $mergedbed $selectedbed 1;";
	}

	#produce bb
	print S1 "cat $selectedbed | grep -v \"#\" | sort -k1,1 -k2,2n > $selectedbed_sorted;$bedtobigbed $selectedbed_sorted ".$tx2ref{$tx}{"chrsize"}." $selectedbb;";
	
	#convert 2 pos? and annotate ...
	

	print S1 "\n";
}



close S1;

#All merge
#Merge called peaks for sample groups

print STDERR "Producing merged bed for all groups, as $allmergedbed\n" if $verbose;
print LOG "Producing merged bed for all groups, as $allmergedbed\n";	


open(S2,">$scriptfile2") || die "Error writing $scriptfile2. $!";

#merge all reprod peaks
print S2 "cat ",join(" ",@groupselectedbeds)," | sort -k1,1 -k2,2n | $mergebed -c 4 -o collapse > $tempfolder/$allmergedbed\_original;";

#rename
print S2 "$processmergebed -i $tempfolder/$allmergedbed\_original -r ",join(",",@groupselectedbeds)," -o $outputfolder/$allmergedbed -c;";

#bb
print S2 "cat $outputfolder/$allmergedbed | grep -v \"#\" | sort -k1,1 -k2,2n > $outputfolder/$allmergedbed_sorted;$bedtobigbed $outputfolder/$allmergedbed_sorted ".$tx2ref{$tx}{"chrsize"}." $outputfolder/$allmergedbb;";
	
#convert to pos
print S2 "$bed2pos $outputfolder/$allmergedbed -o $outputfolder/$allmergedpos;";

#annotate peaks
print S2 "$annotatepeaks $outputfolder/",$allmergedpos," ",$tx2gv{$tx}," -gtf ",$tx2gtf{$tx}," -d $tagdirs > $outputfolder/$allmergedcount_norm 2> $outputfolder/allmergedcount_annotatepeaks_norm_run.log;";
print S2 "$annotatepeaks $outputfolder/",$allmergedpos," ",$tx2gv{$tx}," -gtf ",$tx2gtf{$tx}," -d $tagdirs -raw > $outputfolder/$allmergedcount_raw 2> $outputfolder/allmergedcount_annotatepeaks_raw_run.log;";

#expr-qc

print S2 "$reformatpeakcount -i $outputfolder/$allmergedcount_norm -o $outputfolder/$allmergedcount_norm_reformated -a $allmergedcount_norm_reformated_anno -c $newconfigfile;";
print S2 "$rscript $expr_qc --input $outputfolder/$allmergedcount_norm_reformated --config $newconfigfile --out $outputfolder/expr-qc-all --group $group --geneanno $outputfolder/$allmergedcount_norm;";


print S2 "$reformatpeakcount -i $outputfolder/1000u0d_longest\_$promotermergedcount_norm -o $outputfolder/1000u0d_longest\_$promotermergedcount_norm_reformated -a $outputfolder/1000u0d_longest\_$promotermergedcount_norm_reformated_anno -c $newconfigfile;";
print S2 "$rscript $expr_qc --input $outputfolder/1000u0d_longest\_$promotermergedcount_norm_reformated --config $newconfigfile --out $outputfolder/expr-qc-promoter --group $group --geneanno $outputfolder/1000u0d_longest\_$promotermergedcount_norm;";


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
my $jobname="chipseq-merge-$timestamp";

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


