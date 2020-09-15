#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);
use List::Util qw(sum);


########
#Interface
########


my $version="0.4";

#v0.1a, revised merge.txt title
#v0.2, add tfbs
#v0.2a, improved tfbs
#v0.3, add runmode
#v0.31, solves screen envinroment problem
#v0.4 Firefly, versioning, etc.

my $usage="

chipseq-summary
version: $version
Usage: sbptools chipseq-summary [parameters]

Description: Summarize chipseq-de results and recalculate significance if needed
1. signficant calls from DE results using new cutoffs
2. merge all peaks DE by Signal and by Peak calling
3. for IPA


Parameters:

    --in|-i           Input folder(s). chipseq-de parent folder.
    --output|-o       Output folder

    --fccutoff        Log2 FC cutoff, optional
    --qcutoff         Corrected P cutoff, optional

    --tfbs            Calculate TFBS for revised DE files [F]
	
    --tx|-t           Transcriptome
                        Current support Human.B38.Ensembl84, Mouse.B38.Ensembl84

    --runmode|-r      Where to run the scripts, cluster, local, or none [none]

    Parallel computating parameters
    --task            Number of tasks to be paralleled. By default 4 tasks for local mode, 8 tasks for cluster mode.
    --ncpus           No. of cpus for each task [4]
    --mem|-m          Memory usage for each process, e.g. 100mb, 100gb [40gb]
	
	
";

#--config|-c       Configuration file match the samples in the chipseq-merge folder
#					   first column as sample name.
#--group|-g        Group name in config file for avg calculation

#--chipseq-merge    chipseq-merge folder to retrieve TPM/FPKM data
	
	
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
my $outputfolder;
my $configfile;
my $group;
my $chipseqmerge;
my $fccutoff;
my $qcutoff;
my $geneinput;
my $txinput;
my $verbose=1;
my $tx;
my $tfbs="F";
my $pcol=4; #hidden param
my $runmode="none";
my $jobs=3;

my $task;
my $ncpus=4;
my $mem="40gb";
my $dev=0; #developmental version



GetOptions(
	"in|i=s" => \$inputfolders,
	"output|o=s" => \$outputfolder,
	"config|c=s" => \$configfile,

	"group|g=s" => \$group,
	"chipseq-merge=s" => \$chipseqmerge,
	
	"fccutoff=s" => \$fccutoff,
	"qcutoff=s" => \$qcutoff,
	"pcol=i"=> \$pcol,
	"tfbs=s"=> \$tfbs,
	
	"tx|t=s" => \$tx,	
	"runmode|r=s" => \$runmode,	

	"task=s" => \$task,
	"ncpus=s" => \$ncpus,
	"mem=s" => \$mem,

	"dev" => \$dev,	

	"verbose|v" => \$verbose,
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
	$sbptoolsfolder=get_parent_folder(abs_path(dirname($0)));
}


my $parallel_job="$sbptoolsfolder/parallel-job/parallel-job_caller.pl";
my $mergefiles="$sbptoolsfolder/mergefiles/mergefiles_caller.pl";
my $desummary="$sbptoolsfolder/chipseq-de/summarize_dm_peaks.pl";
my $convertdetobed="$sbptoolsfolder/chipseq-de/convert_chipseq_de_to_bed.pl";
my $motiffinder="$sbptoolsfolder/motif-finder/motif-finder_caller.pl";


########
#Code begins
########

#Create folders

if(!-e $outputfolder) {
	mkdir($outputfolder);
}

$outputfolder = abs_path($outputfolder);

if(!-e "$outputfolder/temp") {
	mkdir("$outputfolder/temp");
}

if(!-e "$outputfolder/forIPA") {
	mkdir("$outputfolder/forIPA");
}

my $logfile="$outputfolder/chipseq-summary_run.log";

my $scriptlocalrun="$outputfolder/chipseq-summary_local_submission.sh";
my $scriptclusterrun="$outputfolder/chipseq-summary_cluster_submission.sh";



#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();
my $timestamp=build_timestamp($now,"long");


print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";
print LOG "\n";


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

#Files for merging
#Examples
#JY_315_Isl.genes.results
#JY_315_Isl.isoforms.results


#Files for merging
#all.reprod.peak.merged.raw.count.DESeq2.Wald.FC0.585.BH.P0.05.summary.txt #DM summary by Signal
#all.reprod.peak.merged.raw.count.DESeq2.Wald.FC0.585.BH.P0.05.anno_rev.txt #DM by Signal
#all.reprod.peak.merged.dm.bycalling.summary.txt #DM Summary by calling

my $scriptfolder="$outputfolder/scripts";

if(!-e $scriptfolder) {
	mkdir($scriptfolder);
}

my $scriptfile1="$scriptfolder/chipseq-summary_run1.sh";



########
#Process
########

print STDERR "\nsbptools chipseq-summary $version running ...\n\n" if $verbose;
print LOG "\nsbptools chipseq-summary $version running ...\n\n";

#open input folders to find comparisons

my %folder2allbycalling;
my %folder2allbysignal;

my %folder2allbycallingsum;
my %folder2allbysignalsum;

my %bysignalsumfiles;
my %bycallingsumfiles;


my %folder2dir;

foreach my $inputfolder (split(",",$inputfolders)) {
	#print STDERR $inputfolder,"#1\n";
	my @folders=glob("$inputfolder/*");
	
	foreach my $folder (@folders) {
		if(-d $folder) {
			#print STDERR $folder,"#2\n";
			if(-e "$folder/chipseq-de_run.log") {
				#chipseq-de folder
				my $foldername=basename($folder);
				
				#foldername needs to be unique
				if(defined $folder2dir{$foldername}) {
					print STDERR "ERROR:$foldername has been used twice:\n",abs_path($folder),"\n",$folder2dir{$foldername},"\n\n";
					print LOG "ERROR:$foldername has been used twice:\n",abs_path($folder),"\n",$folder2dir{$foldername},"\n\n";
					exit;
				}
				
				$folder2dir{$foldername}=abs_path($folder);
				
				my @files=map {basename($_)} glob("$folder/*");
				
				foreach my $file (@files) {
					#by signal
				
					if($file=~/all.reprod.peak.merged.raw.count.+anno_rev.txt/) {
						$folder2allbysignal{$foldername}=$folder2dir{$foldername}."/".$file;
					}
					
					if($file=~/all.reprod.peak.merged.raw.count.+summary.txt/) {
						$folder2allbysignalsum{$foldername}=$folder2dir{$foldername}."/".$file;
						#$bysignalsumfiles{$folder2dir{$foldername}."/".$file}++;
					}					
					
					if($file=~/all.reprod.peak.merged.dm.bycalling.summary.txt/) {
						$folder2allbycallingsum{$foldername}=$folder2dir{$foldername}."/".$file;
						#$bycallingsumfiles{$folder2dir{$foldername}."/".$file}++;
					}
				}
			}
		}
	}
}


unless(keys %folder2allbysignal) {
	print STDERR "ERROR:No chipseq-de results identifed from --input $inputfolders.\n";
	print LOG "ERROR:No chipseq-de results identifed from --input $inputfolders.\n";
	exit;
}

#print out what were found

print STDERR "Printing files identified from input folder(s).\n\n" if $verbose;
print LOG "Printing files identified from input folder(s).\n\n";

foreach my $folder (sort keys %folder2allbysignal) {
	print STDERR $folder,"\t",$folder2allbysignal{$folder},"\n" if $verbose;
	print LOG $folder,"\t",$folder2allbysignal{$folder},"\n";
}

foreach my $folder (sort keys %folder2allbysignalsum) {
	print STDERR $folder,"\t",$folder2allbysignalsum{$folder},"\n" if $verbose;
	print LOG $folder,"\t",$folder2allbysignalsum{$folder},"\n";
}

foreach my $folder (sort keys %folder2allbycallingsum) {
	print STDERR $folder,"\t",$folder2allbycallingsum{$folder},"\n" if $verbose;
	print LOG $folder,"\t",$folder2allbycallingsum{$folder},"\n";
}


##
#Defined gene/peak list
##

#to be implemented




#####
#Work on gene DE files
#####

#check fccutoff and qcutoff
if(defined $fccutoff && length($fccutoff)>0) {
	unless(defined $qcutoff && length($qcutoff)>0) {
		print STDERR "ERROR:Both --fccutoff and --qcutoff need to be redefined.\n";
		print LOG "ERROR:Both --fccutoff and --qcutoff need to be redefined.\n";
		exit;
	}
}

#read file
my %folder2genesig;
my %folder2geneinfo;
my %folder2genetitle;

open(S1,">$scriptfile1") || die "ERROR:Can't write into $scriptfile1.$!\n";

foreach my $folder (sort keys %folder2allbysignal) {

	print STDERR "Processing ",$folder2allbysignal{$folder},"\n" if $verbose;
	print LOG "Processing ",$folder2allbysignal{$folder},"\n";
	
	#my $foripafile="$folder\_GeneDE.foripa.txt";
	
	#print out file ready for IPA
	#system("cut -f 1,2,5 ".$folder2dir{$folder}."/".$folder2allbysignal{$folder}."> $outputfolder/forIPA/$foripafile");
	
	if( defined $fccutoff && length($fccutoff)>0 && defined $qcutoff && length($qcutoff)>0) {
		open(IN,$folder2allbysignal{$folder}) || die $!;
		open(OUT,">$outputfolder/$folder\_BySignal_Reformated.txt") || die $!;

		print STDERR "Convert $folder2allbysignal{$folder} to $outputfolder/$folder\_BySignal_Reformated.txt\n";
		print LOG "Convert $folder2allbysignal{$folder} to $outputfolder/$folder\_BySignal_Reformated.txt\n";
		
		while(<IN>) {
			tr/\r\n//d;
			my @array=split/\t/;
			
			if($_=~/^\W/) {

				if(defined $fccutoff && length($fccutoff)>0 && defined $qcutoff && length($qcutoff)>0) {
					$array[5]="Significance: Log2FC $fccutoff BHP $qcutoff";
				}
				elsif((defined $fccutoff && length($fccutoff)>0) || (defined $qcutoff && length($qcutoff)>0)) {
					print STDERR "ERROR: --fccutoff and --qcutoff have to be both defined.\n\n";
					print LOG "ERROR: --fccutoff and --qcutoff have to be both defined.\n\n";
					exit;
				}
				
				#$folder2genetitle{$folder}=join("\t",@array);
				print OUT join("\t",@array),"\n";
			}
			else {
				
				#record genes
				#unless(defined $geneinput && length($geneinput)>0) {
				#	$genes{$array[0]}++;
				#}
				
				#new fc and q
				if($array[$pcol] ne "NA" && $array[$pcol] ne " ") {
					if(defined $fccutoff && length($fccutoff)>0 && defined $qcutoff && length($qcutoff)>0) {
						if($array[1] >= $fccutoff && $array[$pcol] < $qcutoff) {
							$array[5]=1;
						}
						elsif($array[1] <= -$fccutoff && $array[$pcol] < $qcutoff) {
							$array[5]=-1;
						}
						else{
							$array[5]=0;
						}
					}
				}
				
				#$folder2genesig{$folder}{$array[0]}=$array[5];
				#$folder2geneinfo{$folder}{$array[0]}=join("\t",@array);
				
				print OUT join("\t",@array),"\n";
			}
		}
		close IN;
		close OUT;

		#regenerate By Signal Sum
		system("$desummary -i $outputfolder/$folder\_BySignal_Reformated.txt --tx $tx --out $outputfolder/$folder\_BySignal_Reformated_summary.txt");

		#redefine sum file
		$folder2allbysignalsum{$folder}="$outputfolder/$folder\_BySignal_Reformated_summary.txt";

	
		#$bysignalsumfiles{"$outputfolder/$folder\_BySignal_Reformated_summary.txt"}++;

	}
	else {
		system("cp $folder2allbysignal{$folder} $outputfolder/$folder\_BySignal_Reformated.txt");
		system("cp $folder2allbysignalsum{$folder} $outputfolder/$folder\_BySignal_Reformated_summary.txt");
		
		print STDERR "cp $folder2allbysignal{$folder} $outputfolder/$folder\_BySignal_Reformated.txt\n";
		print LOG "cp $folder2allbysignal{$folder} $outputfolder/$folder\_BySignal_Reformated.txt\n";
		
		print STDERR "cp $folder2allbysignalsum{$folder} $outputfolder/$folder\_BySignal_Reformated_summary.txt\n";
		print LOG "cp $folder2allbysignalsum{$folder} $outputfolder/$folder\_BySignal_Reformated_summary.txt\n";
	}
						
	
	#for IPA here
	
	#Directional total
	system("printf \"Gene\tLog2FC\n\" > $outputfolder/forIPA/$folder\_BySignal_DirectionalTotal_forIPA.txt;cut -f 1,8 $outputfolder/$folder\_BySignal_Reformated_summary.txt | sed '1d' >> $outputfolder/forIPA/$folder\_BySignal_DirectionalTotal_forIPA.txt");
	print LOG "printf \"Gene\tLog2FC\n\" > $outputfolder/forIPA/$folder\_BySignal_DirectionalTotal_forIPA.txt;cut -f 1,8 $outputfolder/$folder\_BySignal_Reformated_summary.txt | sed '1d' >> $outputfolder/forIPA/$folder\_BySignal_DirectionalTotal_forIPA.txt\n";
	
	#Promoter
	system("printf \"Gene\tLog2FC\n\" > $outputfolder/forIPA/$folder\_BySignal_Promoter_forIPA.txt;cut -f 1,6 $outputfolder/$folder\_BySignal_Reformated_summary.txt | sed '1d' >> $outputfolder/forIPA/$folder\_BySignal_Promoter_forIPA.txt");
	print LOG "printf \"Gene\tLog2FC\n\" > $outputfolder/forIPA/$folder\_BySignal_Promoter_forIPA.txt;cut -f 1,6 $outputfolder/$folder\_BySignal_Reformated_summary.txt | sed '1d' >> $outputfolder/forIPA/$folder\_BySignal_Promoter_forIPA.txt\n";
	
	
	#TFBS
	if($tfbs eq "T") {
		#TFBS
		mkdir("$outputfolder/$folder\_deboth_tfbs");
		mkdir("$outputfolder/$folder\_deup_tfbs");
		mkdir("$outputfolder/$folder\_dedown_tfbs");
		
		#change the implementation, so that it can be better paralleled
		system("$convertdetobed -i $outputfolder/$folder\_BySignal_Reformated.txt -o $outputfolder/$folder\_BySignal_Reformated_de.bed;");
	
		#needs to be updated by adding background?
		print S1 "$motiffinder -i $outputfolder/$folder\_BySignal_Reformated_de.bed -o $outputfolder/$folder\_deboth_tfbs/ -t $tx;\n";
		print S1 "$motiffinder -i $outputfolder/$folder\_BySignal_Reformated_de_up.bed -o $outputfolder/$folder\_deup_tfbs/ -t $tx;\n";
		print S1 "$motiffinder -i $outputfolder/$folder\_BySignal_Reformated_de_down.bed -o $outputfolder/$folder\_dedown_tfbs/ -t $tx;\n";
		
	}					
}

close S1;

#copy gene names 
system("cut -f 1 ".$tx2ref{$tx}{"geneanno"}." > $outputfolder/$tx.genes.txt");

######
#merge bysignal summary
######
#Intergenic	TTS	exon	intron	promoter-TSS	Undirectional Total	Directional Total
my @cates=("Intergenic","TTS","exon","intron","promoter-TSS","Undirectional_Total","Directional_Total");

open(OUT,">$outputfolder/temp/merge_title.txt") || die $!;
print OUT "Gene\t",join("\t",map {join("\t",$_,@cates,@cates[0..4],@cates[0..4])} sort keys %folder2allbysignal),"\n";
close OUT;


#Merge BySignal summary	
print STDERR "$mergefiles -m $outputfolder/$tx.genes.txt -i ".join(",",map {$folder2allbysignalsum{$_}} sort keys %folder2allbysignal)." -o $outputfolder/temp/chipseq-summary_BySignal_Reformated_summary_merged_temp.txt"."\n";
#
system("$mergefiles -m $outputfolder/$tx.genes.txt -i ".join(",",map {$folder2allbysignalsum{$_}} sort keys %folder2allbysignal)." -o $outputfolder/temp/chipseq-summary_BySignal_Reformated_summary_merged_temp.txt"); #create a temporary file first
#make title
system("cat $outputfolder/temp/merge_title.txt > $outputfolder/chipseq-summary_BySignal_Reformated_summary_merged.txt");

system("awk '{if(NR>1)print}' $outputfolder/temp/chipseq-summary_BySignal_Reformated_summary_merged_temp.txt >> $outputfolder/chipseq-summary_BySignal_Reformated_summary_merged.txt");



#folder title
open(OUT,">$outputfolder/temp/samples.txt") || die $!;
print OUT "Gene\t",join("\t",sort keys %folder2allbysignal),"\n";
close OUT;


#bysignal summary, by category

for(my $num=0;$num<@cates;$num++) {
	#get colnums for each category
	my @colnums=map {$_*18-15+$num} 1..scalar(keys %folder2allbysignal);
	system("cut -f 1,".join(",",@colnums)." $outputfolder/chipseq-summary_BySignal_Reformated_summary_merged.txt > $outputfolder/temp/chipseq-summary_BySignal_Reformated_summary_merged_$cates[$num]_temp.txt");
	
	system("cat $outputfolder/temp/samples.txt > $outputfolder/chipseq-summary_BySignal_Reformated_summary_merged_$cates[$num].txt");
	
	system("awk '{if(NR>1)print}' $outputfolder/temp/chipseq-summary_BySignal_Reformated_summary_merged_$cates[$num]_temp.txt >> $outputfolder/chipseq-summary_BySignal_Reformated_summary_merged_$cates[$num].txt");
}


######
#merge bycalling summary
######

system("$mergefiles -m $outputfolder/$tx.genes.txt -i ".join(",",map {$folder2allbycallingsum{$_}} sort keys %folder2allbysignal)." -o $outputfolder/temp/chipseq-summary_ByCalling_Reformated_summary_merged_temp.txt");

system("cat $outputfolder/temp/merge_title.txt > $outputfolder/chipseq-summary_ByCalling_Reformated_summary_merged.txt");

system("awk '{if(NR>1)print}' $outputfolder/temp/chipseq-summary_ByCalling_Reformated_summary_merged_temp.txt >> $outputfolder/chipseq-summary_ByCalling_Reformated_summary_merged.txt");


#bysignal summary, by category
#Intergenic	TTS	exon	intron	promoter-TSS	Undirectional Total	Directional Total

for(my $num=0;$num<@cates;$num++) {
	#get colnums for each category
	my @colnums=map {$_*18-15+$num} 1..scalar(keys %folder2allbysignal);
	system("cut -f 1,".join(",",@colnums)." $outputfolder/chipseq-summary_ByCalling_Reformated_summary_merged.txt > $outputfolder/temp/chipseq-summary_ByCalling_Reformated_summary_merged_$cates[$num]_temp.txt");
	
	system("cat $outputfolder/temp/samples.txt > $outputfolder/chipseq-summary_ByCalling_Reformated_summary_merged_$cates[$num].txt");
	
	system("awk '{if(NR>1)print}' $outputfolder/temp/chipseq-summary_ByCalling_Reformated_summary_merged_$cates[$num]_temp.txt >> $outputfolder/chipseq-summary_ByCalling_Reformated_summary_merged_$cates[$num].txt");	
}




####
#Print OUT gene annotation
####

print STDERR "Copy gene annotation.\n" if $verbose;
print LOG "Copy gene annotation.\n";

system("cp ".$tx2ref{$tx}{"geneanno"}." $outputfolder/geneanno.txt");
#system("$mergefiles -m $outputfolder/genes_sel.txt -i ".$tx2ref{$tx}{"geneanno"}." -o $outputfolder/geneanno_sel.txt");



#######
#Run mode
#######


if($tfbs eq "T") {

	open(LOUT,">$scriptlocalrun") || die "ERROR:can't write to $scriptlocalrun. $!";
	open(SOUT,">$scriptclusterrun") || die "ERROR:can't write to $scriptclusterrun. $!";


	my @scripts_all=($scriptfile1);


	#print out command for local and cluster parallel runs
	my $jobnumber=0;
	my $jobname="chipseq-summary-$timestamp";

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
}




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

sub cal_group_avg {
	my ($nums,$names)=@_;
	
	my @avgs;
	my @names_uniq=sort(uniq(@$names));
	my %name2pos;
	
	for(my $n=0;$n<@$names;$n++) {
		push @{$name2pos{$names->[$n]}},$n;
	}
	
	foreach my $name (@names_uniq) {
		my @selnums=@{$nums}[@{$name2pos{$name}}];
		push @avgs,mean(@selnums);
	}
	
	return @avgs;
}

sub mean {
	return sum(@_)/@_;
}

sub uniq {
	my @array=@_;
	my %hash;
	
	foreach my $item (@array) {
		$hash{$item}++;
	}
	
	return sort keys %hash;
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




