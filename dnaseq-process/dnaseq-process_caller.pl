#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);

#CutAdapt+FASTQC+BWA+GATK4


########
#Interface
########


my $version="0.13";

#v0.11, fix file removing bugs
#v0.12, add TMB summary
#v0.13, rm temp files earlier


my $usage="

dnaseq-process
version: $version
Usage: sbptools dnaseq-process [parameters]

Description: 


Parameters:

    --config|-c       Configuration file
                           first column as sample name.
                           Controled headers include:
                               fastq1,fastq2 columns for fastqs
                               index for sample index used for folder names
    --samplepair|-p   Tumor and Normal pairing (Optional)
    --output|-o       Output folder

    --tx|-t           Transcriptome [Human.B38.Ensembl84]
                        Currently only support Human.B38.Ensembl84

    --dp|-d           Sequencing coverage, DP, filter >= [5]
    --af|-a           Allele Frequency, AF, filter >= [0.05]
	
    --task            Number of tasks to be paralleled. By default 4 tasks for local mode, 8 tasks for cluster mode.
    --ncpus           No. of cpus for each task [8]
    --mem|-m          Memory usage for each process, e.g. 100mb, 100gb [40gb]
	
    --runmode|-r      Where to run the scripts, local, cluster or none [none]
                            local is to run locally using \"parallel\", recommended for Falco
                            cluster is to submit jobs to PBS queue in the HPC, recommended for Firefly
                            none is to generate scripts only, after that,
                                   you can use \"sh rnaseq-merge_local_submission.sh\" to run locally, or
                                   you can use \"sh rnaseq-merge_cluster_submission.sh\" to submit job to PBS

								   
";

#--verbose|-v      Verbose

unless (@ARGV) {
	print STDERR $usage;
	exit;
}

my $params=join(" ",@ARGV);
#then call different scripts



########
#Parameters
########


my $configfile;
my $samplepairfile;
my $outputfolder;
my $verbose=1;
my $tx="Human.B38.Ensembl84";
my $dpfilter=5;
my $affilter=0.05;
my $task;
my $ncpus=8;
my $runbamcoverage="F";
my $mem="40gb";
my $runmode="none";

my $dev=0; #developmental version


GetOptions(
	"config|c=s" => \$configfile,
	"samplepair|p=s" => \$samplepairfile,
	"output|o=s" => \$outputfolder,
	"tx|t=s" => \$tx,
	"dp|d=s" => \$dpfilter,	
	"af|a=s" => \$affilter,		
	"task=s" => \$task,
	"ncpus=s" => \$ncpus,
	"mem=s" => \$mem,
	"bamcoverage=s" => \$runbamcoverage,	
	"runmode|r=s" => \$runmode,		
	"verbose|v" => \$verbose,
	"dev" => \$dev,		
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
my $rnaseq_var_filter="$sbptoolsfolder/rnaseq-var/rnaseq-var_filter.pl";
my $reformatsnpeffvcf="$sbptoolsfolder/rnaseq-var/reformat_snpeff_vcf.pl";

#my $cutadapt=find_program("/apps/python-3.5.2/bin/cutadapt");
my $fastqc=find_program("/apps/FastQC/fastqc");

my $trimmomatic=find_program("/apps/Trimmomatic-0.38/trimmomatic-0.38.jar");
my $bwa=find_program("/apps/bwa-0.7.15/bwa");
my $samtools=find_program("/apps/samtools-1.3.1/bin/samtools");
my $gatk=find_program("/apps/gatk-4.1.8.1/gatk");
my $picard="/apps/gatk-4.1.8.1/picard.jar";

my $java=find_program("/home/apps/jdk1.8.0_231/bin/java");

my $bamcoverage=find_program("/apps/python-3.5.2/bin/bamCoverage");

my $snpeff="/apps/snpEff-5/snpEff.jar";
my $snpsift="/apps/snpEff-5/SnpSift.jar";



#######
#Output folder
#######

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


my $logfile="$outputfolder/dnaseq-process_run.log";
my $newconfigfile="$outputfolder/dnaseq-process_config.txt";

my $scriptfile1="$scriptfolder/dnaseq-process_run1.sh";
my $scriptfile2="$scriptfolder/dnaseq-process_run2.sh";


my $scriptlocalrun="$outputfolder/dnaseq-process_local_submission.sh";
my $scriptclusterrun="$outputfolder/dnaseq-process_cluster_submission.sh";



#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();
my $timestamp=build_timestamp($now,"long");

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";
#print LOG "$cutadapt version:", getsysoutput("$cutadapt --version"),"\n";
print LOG "$fastqc version:", getsysoutput("$fastqc --version"),"\n";
print LOG "$bwa version:", getsysoutput("$bwa 2>&1 | grep Version "),"\n";
print LOG "$gatk version:", getsysoutput("$gatk --version 2>&1 | grep \"The Genome Analysis Toolkit (GATK)\""),"\n";
print LOG "$picard version:", getsysoutput("$java -jar $picard MergeBamAlignment --version 2>&1 "),"\n";


print LOG "\n";

#Annotation files

my %tx2ref=(
	"Human.B38.Ensembl84"=> { 
		"star"=>"/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl84_STAR",
		"rsem"=>"/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl84_STAR/Human_RSEM",
		"chrsize"=>"/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl84_STAR/chrNameLength.txt",
		"fasta"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.dna.primary_assembly_ucsc.fa",
		"gtf"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc.gtf",
		"homeranno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_homeranno.txt",
		"geneanno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_gene_annocombo_rev.txt",
		"txanno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_tx_annocombo.txt",
		"gnomad_common001"=>"/data/jyin/Databases/gnomad/fromGATK4/af-only-gnomad.hg38.common001_refchrs.vcf.gz",
		"gnomad_common1e4"=>"/data/jyin/Databases/gnomad/fromGATK4/af-only-gnomad.hg38.common1e4_refchrs.vcf.gz",
		"knownsnp"=>"/data/jyin/Databases/SNP/00-All.vcf.gz",
		"genomeversion"=>"hg38",
	}
);


########
#Process
########

print STDERR "\nsbptools dnaseq-process $version running ...\n\n" if $verbose;
print LOG "\nsbptools dnaseq-process $version running ...\n\n";


if(defined $tx2ref{$tx}) {
	print STDERR "Starting analysis using $tx.\n\n" if $verbose;
	print LOG "Starting analysis using $tx.\n\n";
}
else {
	print STDERR "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n" if $verbose;
	print LOG "ERROR:$tx not defined. Currently only supports ",join(",",sort keys %tx2ref),"\n\n";
	exit;
}


#open config files to find fastqs

my %sample2fastq;
my %sample2indexname;
my %configattrs;

my $fileline=0;
my $newconfigfiletitle;
open(IN,$configfile) || die "Error reading $configfile. $!";
open(OUT,">$newconfigfile") || die "Error reading $newconfigfile. $!";

while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	if($fileline==0) {
		for(my $num=0;$num<@array;$num++) {
			#case sensitive match to key words, Sample, FASTQ1, FASTQ2, Index (for foldernames)
			$configattrs{uc $array[$num]}=$num;
		}
		
		#should always use first column as sample name and index to maintain continuous workflow
		
		#first column as sample
		
		#check title
		unless(defined $configattrs{"FASTQ1"}) {
			#SAMPLE and FASTQ1 have to be defined. The others are optional
			print STDERR "ERROR: FASTQ1 need to be defined in $configfile. Input config includes:",join(",",map {uc $_} @array),"\n";
			print LOG "ERROR: FASTQ1 need to be defined in $configfile. Input config includes:",join(",",map {uc $_} @array),"\n";
			exit;
		}
		else {
			print STDERR "Input config $configfile includes:",join(",",map {uc $_} @array),"\n\n" if $verbose;
			print LOG "Input config $configfile includes:",join(",",map {uc $_} @array),"\n\n";
		}
		
		#version 0.2, removed self defined INDEX
		
		if(defined $configattrs{"INDEX"}) {
			$newconfigfiletitle=join("\t","INDEX",map {uc $_} (splice @array,$configattrs{"INDEX"},1));
		}
		else {
			$newconfigfiletitle=join("\t","INDEX",map {uc $_} @array);
		}
		
		print OUT $newconfigfiletitle,"\n";
		
		print STDERR "New config $newconfigfile includes:",join(",",split("\t",$newconfigfiletitle)),"\n\n" if $verbose;
		print LOG "New config $newconfigfile includes:",join(",",split("\t",$newconfigfiletitle)),"\n\n";		
	
	}
	else {
		
		my $indexname;
		my $newconfigline;
		
		#need to change non word chars to words
		#if defined INDEX, use index, otherwise, use first column
		if(defined $configattrs{"INDEX"}) {
			$indexname=$array[$configattrs{"INDEX"}];
		}
		else {
			$indexname=$array[0]; #change first column as sample column in v0.2
			$indexname=~s/[^\w-]/_/g;
		}
		
		
		$sample2indexname{$array[0]}=$indexname;

		
		#fastq1
		if(-e $array[$configattrs{"FASTQ1"}]) {
			push @{$sample2fastq{$indexname}},$array[$configattrs{"FASTQ1"}];
		}
		else {
			print STDERR "\n\nERROR:",$array[$configattrs{"FASTQ1"}]," doesn't exist. Please check $configfile setting.\n\n";
			print LOG "\n\nERROR:",$array[$configattrs{"FASTQ1"}]," doesn't exist. Please check $configfile setting.\n\n";
		}
		
		#fastq2 (Read2)
		if(defined $configattrs{"FASTQ2"}) {
			if(-e $array[$configattrs{"FASTQ2"}]) {
				push @{$sample2fastq{$indexname}},$array[$configattrs{"FASTQ2"}];
			}
			else {
				print STDERR "\n\nERROR:",$array[$configattrs{"FASTQ2"}]," doesn't exist. Please check $configfile setting.\n\n";
				print LOG "\n\nERROR:",$array[$configattrs{"FASTQ2"}]," doesn't exist. Please check $configfile setting.\n\n";	
			}
		}
		
		if(defined $configattrs{"INDEX"}) {
			$newconfigline=join("\t",$indexname,(splice @array,$configattrs{"INDEX"},1));
		}
		else {
			$newconfigline=join("\t",$indexname,@array);
		}
		
		print OUT $newconfigline,"\n";
		
	}
	
	$fileline++;
}

close IN;
close OUT;


#----------------
my %samplepairs;

if(defined $samplepairfile && length($samplepairfile) >0) {
	open(IN,$samplepairfile) || die "ERROR:Can't read $samplepairfile.\n";
	while(<IN>) {
		tr/\r\n//d;
		my @array=split/\t/;
		$samplepairs{$array[2]}=[@array[0,1]]; #Tumor, Normal, Name #three columns file
	}
	close IN;
}

#----------------
#create folder
print STDERR scalar(keys %sample2indexname)," samples identified from $configfile, including:\nINDEX\tSAMPLE\n",join("\n",map {$sample2indexname{$_}."\t".$_} sort keys %sample2indexname),"\n\n" if $verbose;
print LOG scalar(keys %sample2indexname)," samples identified from $configfile, including:\nINDEX\tSAMPLE\n",join("\n",map {$sample2indexname{$_}."\t".$_} sort keys %sample2indexname),"\n\n" if $verbose;

#print STDERR scalar(keys %sample2fastq)," samples identified from $configfile, including:",join("\n",sort keys %sample2fastq),"\n\n" if $verbose;
#print LOG scalar(keys %sample2fastq)," samples identified from $configfile, including:",join("\n",sort keys %sample2fastq),"\n\n";


print STDERR "Make folders for different samples.\n\n" if $verbose;
print LOG "Make folders for different samples.\n\n";

foreach my $sample (sort keys %sample2fastq) {
	my $samplefolder="$outputfolder/$sample"; #edited to enable space in folder name
	if(!-e	$samplefolder) {
		mkdir $samplefolder;
	}
}

#----------------
# Sequence alignment and calibration for BAM files for each sample
#----------------

my %sample2workflow;
my %tempfiles2rm1;

foreach my $sample (sort keys %sample2fastq) {
	#my $samplefolder="$outputfolder/$sample";


	mkdir("$outputfolder/$sample/alignment");
	mkdir("$outputfolder/$sample/gatk4");
	mkdir("$outputfolder/$sample/snpanno");

	print STDERR "Printing variant calling script.\n\n" if $verbose;
	print LOG "Printing variant calling script.\n\n";

	my $bwacommand;

	if(defined $configattrs{"FASTQ2"}) {
		#PE

		my $fastq1=$sample2fastq{$sample}[0];
		my $fastq2=$sample2fastq{$sample}[1]; #
		
		my $fastq1trim="$outputfolder/$sample/alignment/".basename($fastq1);
		$fastq1trim=~s/\.fastq\.gz/_trimmed.fastq.gz/;
		push @{$sample2fastq{$sample}},$fastq1trim;

		my $fastq2trim="$outputfolder/$sample/alignment/".basename($fastq2);
		$fastq2trim=~s/\.fastq\.gz/_trimmed.fastq.gz/;
		push @{$sample2fastq{$sample}},$fastq2trim;
		
		my $fastq1unpaired="$outputfolder/$sample/alignment/".basename($fastq1);
		$fastq1unpaired=~s/\.fastq\.gz/_unpaired.fastq.gz/;
		push @{$sample2fastq{$sample}},$fastq1unpaired;

		my $fastq2unpaired="$outputfolder/$sample/alignment/".basename($fastq2);
		$fastq2unpaired=~s/\.fastq\.gz/_unpaired.fastq.gz/;
		push @{$sample2fastq{$sample}},$fastq2unpaired;
		
		#Trimmomatic, light weight trimming
		#------------------
		$sample2workflow{$sample}.="$java -jar $trimmomatic PE ".$sample2fastq{$sample}[0]." ".$sample2fastq{$sample}[1]." ".$sample2fastq{$sample}[2]." ".$sample2fastq{$sample}[4]." ".$sample2fastq{$sample}[3]." ".$sample2fastq{$sample}[5]." ILLUMINACLIP:/apps/Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10 > $outputfolder/$sample/alignment/trimmomatic.log 2>&1;";
		
		$tempfiles2rm1{$sample}{"$outputfolder/$sample/alignment/*.fastq.gz"}++;
		
		
		#Fastqc
		#------------------
		$sample2workflow{$sample}.="$fastqc --nogroup -o $outputfolder/$sample/alignment/ -f fastq ".$sample2fastq{$sample}[2]." ".$sample2fastq{$sample}[3].";";

		#Convert Fastq to uBAM for mergealignment step
		#------------------
		$sample2workflow{$sample}.="$java -Xmx8G -jar $picard FastqToSam  FASTQ=".$sample2fastq{$sample}[2]." FASTQ2=".$sample2fastq{$sample}[3]." OUTPUT=$outputfolder/$sample/alignment/$sample\_fastq2bam.bam  READ_GROUP_NAME=ReadGroupName SAMPLE_NAME=$sample LIBRARY_NAME=LibraryName PLATFORM=illumina > $outputfolder/$sample/alignment/FastqToSam.log 2>&1;";

		#temp file: $outputfolder/$sample/$sample\_fastq2bam.bam
		$tempfiles2rm1{$sample}{"$outputfolder/$sample/alignment/$sample\_fastq2bam.bam"}++;


		#BWA alignment and convert to bam
		#------------------
		$sample2workflow{$sample}.="$bwa mem -v 3 -t $ncpus -Y ".$tx2ref{$tx}{"fasta"}." ".$sample2fastq{$sample}[2]." ".$sample2fastq{$sample}[3]." 2> $outputfolder/$sample/alignment/bwa.log | $samtools view -1 - > $outputfolder/$sample/alignment/$sample\_unmerged.bam;";
		
		#temp file $outputfolder/$sample/$sample\_unmerged.bam
		$tempfiles2rm1{$sample}{"$outputfolder/$sample/alignment/$sample\_unmerged.bam"}++;
		
		$bwacommand="$bwa mem -v 3 -t $ncpus -Y ".$tx2ref{$tx}{"fasta"}." ".$sample2fastq{$sample}[2]." ".$sample2fastq{$sample}[3];
	}
	else {
		#SE
	
	
	}
	
	#now move to gatk4 folder
	#MergeBamAlignment
	#------------------
	$sample2workflow{$sample}.="$java -Dsamjdk.compression_level=5 -Xmx8G -jar $picard MergeBamAlignment --VALIDATION_STRINGENCY SILENT  --EXPECTED_ORIENTATIONS FR  --ATTRIBUTES_TO_RETAIN X0  --ALIGNED_BAM $outputfolder/$sample/alignment/$sample\_unmerged.bam --UNMAPPED_BAM $outputfolder/$sample/alignment/$sample\_fastq2bam.bam  --OUTPUT $outputfolder/$sample/gatk4/$sample.aligned.unsorted.bam  --REFERENCE_SEQUENCE ".$tx2ref{$tx}{"fasta"}."  --PAIRED_RUN true  --SORT_ORDER \"unsorted\"  --IS_BISULFITE_SEQUENCE false  --ALIGNED_READS_ONLY false  --CLIP_ADAPTERS false  --MAX_RECORDS_IN_RAM 2000000  --ADD_MATE_CIGAR true  --MAX_INSERTIONS_OR_DELETIONS -1  --PRIMARY_ALIGNMENT_STRATEGY MostDistant  --PROGRAM_RECORD_ID \"bwamem\"  --PROGRAM_GROUP_VERSION \"1.7\"  --PROGRAM_GROUP_COMMAND_LINE \"$bwa mem -v 3 -t $ncpus -Y ".$tx2ref{$tx}{"fasta"}." ".$sample2fastq{$sample}[2]." ".$sample2fastq{$sample}[3]."\"  --PROGRAM_GROUP_NAME \"bwamem\"  --UNMAPPED_READ_STRATEGY COPY_TO_TAG  --ALIGNER_PROPER_PAIR_FLAGS true  --UNMAP_CONTAMINANT_READS true > $outputfolder/$sample/gatk4/mergebamalignment.log 2>&1;";
	
	#temp file $outputfolder/$sample/$sample.aligned.unsorted.bam
	$tempfiles2rm1{$sample}{"$outputfolder/$sample/gatk4/$sample.aligned.unsorted.bam"}++;
	
	
	#MarkDuplicates
	#------------------
	$sample2workflow{$sample}.="$java -Dsamjdk.compression_level=5 -Xms4000m -jar $picard MarkDuplicates INPUT=$outputfolder/$sample/gatk4/$sample.aligned.unsorted.bam OUTPUT=$outputfolder/$sample/gatk4/$sample.aligned.unsorted.dupliates_marked.bam METRICS_FILE=$outputfolder/$sample/gatk4/$sample.dupliates_metrics.txt VALIDATION_STRINGENCY=SILENT OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 ASSUME_SORT_ORDER=queryname CLEAR_DT=false ADD_PG_TAG_TO_READS=false > $outputfolder/$sample/gatk4/markduplicates.log 2>&1;";
	
	#temp file $outputfolder/$sample/gatk4/$sample.aligned.unsorted.dupliates_marked.bam
	$tempfiles2rm1{$sample}{"$outputfolder/$sample/gatk4/$sample.aligned.unsorted.dupliates_marked.bam"}++;
	
	
	#Sort and fix tags
	#------------------
	$sample2workflow{$sample}.="$java -Dsamjdk.compression_level=5 -Xms4g -jar $picard SortSam --INPUT $outputfolder/$sample/gatk4/$sample.aligned.unsorted.dupliates_marked.bam --OUTPUT /dev/stdout --SORT_ORDER coordinate --CREATE_INDEX false --CREATE_MD5_FILE false | $java -Dsamjdk.compression_level=5 -Xms4g -jar $picard SetNmMdAndUqTags --INPUT /dev/stdin --OUTPUT $outputfolder/$sample/gatk4/$sample.aligned.dupliates_marked.sorted.bam --CREATE_INDEX true --CREATE_MD5_FILE true --REFERENCE_SEQUENCE ".$tx2ref{$tx}{"fasta"}." > $outputfolder/$sample/gatk4/sortandfixtags.log 2>&1;";
	
	#temp file? $outputfolder/$sample/gatk4/$sample.aligned.dupliates_marked.sorted.bam
	#may be final file if no calibration is needed
	$tempfiles2rm1{$sample}{"$outputfolder/$sample/gatk4/$sample.aligned.dupliates_marked.sorted.bam"}++;

	#BaseRecalibrator
	#------------------
	$sample2workflow{$sample}.="$gatk --java-options \"-Xms4g\" BaseRecalibrator -R ".$tx2ref{$tx}{"fasta"}." -I $outputfolder/$sample/gatk4/$sample.aligned.dupliates_marked.sorted.bam --use-original-qualities -O $outputfolder/$sample/gatk4/recal_data.csv --known-sites ".$tx2ref{$tx}{"gnomad_common1e4"}." > $outputfolder/$sample/gatk4/BaseRecalibrator.log 2>&1;";
	
	#Apply BQSR
	#------------------
	$sample2workflow{$sample}.="$gatk --java-options \"-Xms4g\" ApplyBQSR -R ".$tx2ref{$tx}{"fasta"}." -I $outputfolder/$sample/gatk4/$sample.aligned.dupliates_marked.sorted.bam -O $outputfolder/$sample/gatk4/$sample.aligned.dupliates_marked.recalibrated.bam -bqsr $outputfolder/$sample/gatk4/recal_data.csv --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 --add-output-sam-program-record --create-output-bam-md5 --use-original-qualities > $outputfolder/$sample/gatk4/ApplyBQSR.log 2>&1;";
	
	#final bam, used for Mutect2. Only leave this file for record
	#$outputfolder/$sample/$sample.aligned.dupliates_marked.recalibrated.bam
	
	#rm temporary files
	$sample2workflow{$sample}.="rm ".join(" ",sort keys %{$tempfiles2rm1{$sample}}).";";
	
}


#----------------
# Call variants and annotation
#----------------

my %sample2workflow2; #mutect2 for tumor only or tumor/normal
my %tempfiles2rm2;

#The second step of the workflow needs the first steps to be all finished 

if(keys %samplepairs>0) {
	#To be implemented


}
else {

	foreach my $sample (sort keys %sample2fastq) {
		#For tumor, Mutect2. For Germline, Haplotyper
		#Current workflow doesn't have normal match...

		#Mutect2
		#------------------
		#for tumor only
		$sample2workflow2{$sample}.="$gatk --java-options \"-Xms4g\" Mutect2 -R ".$tx2ref{$tx}{"fasta"}." -I $outputfolder/$sample/gatk4/$sample.aligned.dupliates_marked.recalibrated.bam -O $outputfolder/$sample/gatk4/$sample.vcf.gz --germline-resource ".$tx2ref{$tx}{"gnomad_common1e4"}." --f1r2-tar-gz $outputfolder/$sample/gatk4/f1r2.tar.gz > $outputfolder/$sample/gatk4/Mutect2.log 2>&1;";
		
		#GetPileupSummaries #variants for contamination/common snp needed #minimum 30g memory ...
		#------------------
		
		$sample2workflow2{$sample}.="$gatk --java-options \"-Xms40g\" GetPileupSummaries -R ".$tx2ref{$tx}{"fasta"}." -I $outputfolder/$sample/gatk4/$sample.aligned.dupliates_marked.recalibrated.bam -V ".$tx2ref{$tx}{"gnomad_common1e4"}." -L ".$tx2ref{$tx}{"gnomad_common1e4"}." -O $outputfolder/$sample/gatk4/tumor-pileups.table > $outputfolder/$sample/gatk4/GetPileupSummaries.log 2>&1;";

		
		#Learn Read Orientation Model
		#------------------
		$sample2workflow2{$sample}.="$gatk --java-options \"-Xms4g\" LearnReadOrientationModel -I $outputfolder/$sample/gatk4/f1r2.tar.gz -O $outputfolder/$sample/gatk4/artifact-priors.tar.gz > $outputfolder/$sample/gatk4/LearnReadOrientationModel.log 2>&1;";
		
		#Calculate contamination
		#------------------
		$sample2workflow2{$sample}.="$gatk --java-options \"-Xmx4000m\" CalculateContamination -I $outputfolder/$sample/gatk4/tumor-pileups.table -O $outputfolder/$sample/gatk4/contamination.table --tumor-segmentation $outputfolder/$sample/gatk4/segments.table > $outputfolder/$sample/gatk4/CalculateContamination.log 2>&1;";
		
		#Filter
		#------------------
		$sample2workflow2{$sample}.="$gatk --java-options \"-Xmx4000m\" FilterMutectCalls -V $outputfolder/$sample/gatk4/$sample.vcf.gz -R ".$tx2ref{$tx}{"fasta"}." -O $outputfolder/$sample/gatk4/$sample-filtered.vcf.gz --contamination-table $outputfolder/$sample/gatk4/contamination.table --tumor-segmentation $outputfolder/$sample/gatk4/segments.table --ob-priors $outputfolder/$sample/gatk4/artifact-priors.tar.gz --filtering-stats $outputfolder/$sample/gatk4/filtering.stats.txt > $outputfolder/$sample/gatk4/FilterMutectCalls.log 2>&1;";

		#Funcotator (not working)
		#------------------
		
		#$sample2workflow2{$sample}.="$gatk --java-options \"-Xmx4000m\" Funcotator --data-sources-path /data/jyin/Pipeline_Test/dnaseq-process/funcotator_dataSources.v1.6.20190124s --ref-version hg38 --output-file-format VCF -R ".$tx2ref{$tx}{"fasta"}." -V $outputfolder/$sample/$sample-filtered.vcf.gz -O $outputfolder/$sample/$sample-funcotated.vcf;";
		
		#Sift
		#------------------
		
		#Outputname is same with RNASeq var
		
		$sample2workflow2{$sample}.="$java -jar $snpsift annotate -noInfo -v ".$tx2ref{$tx}{"knownsnp"}." $outputfolder/$sample/gatk4/$sample-filtered.vcf.gz > $outputfolder/$sample/snpanno/${sample}.filtered.sift.annotated.vcf 2> $outputfolder/$sample/snpanno/snpsift.log;";
		
		$tempfiles2rm2{$sample}{"$outputfolder/$sample/snpanno/${sample}.filtered.sift.annotated.vcf"}++;
		
		#Manual filter, because SNPeff doesn't differentiate filtered/unfiltered
		#filter by common snp, 
		#------------------
		$sample2workflow2{$sample}.="perl $rnaseq_var_filter -i $outputfolder/$sample/snpanno/${sample}.filtered.sift.annotated.vcf -d $dpfilter -a $affilter -c ".substr($tx2ref{$tx}{"gnomad_common1e4"},0,length($tx2ref{$tx}{"gnomad_common1e4"})-3)." -o $outputfolder/$sample/snpanno/${sample}.filtered-cleaned.sift.annotated.vcf 2> $outputfolder/$sample/snpanno/var_filter.log;";
		
		$tempfiles2rm2{$sample}{"$outputfolder/$sample/snpanno/${sample}.filtered-cleaned.sift.annotated.vcf"}++;
		
		#SNPeff
		#------------------
		#filtered-cleaned
		$sample2workflow2{$sample}.="$java -jar $snpeff -v ".$tx2ref{$tx}{"genomeversion"}." $outputfolder/$sample/snpanno/${sample}.filtered-cleaned.sift.annotated.vcf -stats $outputfolder/$sample/snpanno/$sample.filtered-cleaned.snpEff_summary.html > $outputfolder/$sample/snpanno/${sample}.filtered-cleaned.snpeff.sift.annotated.vcf 2> $outputfolder/$sample/snpanno/snpeff1.log;";

		$sample2workflow2{$sample}.="mv $outputfolder/$sample/snpanno/snpEff_summary.genes.txt $outputfolder/$sample/snpanno/$sample.filtered-cleaned.snpEff_summary.genes.txt;";
		
		#reformat snpeff ANN column
		$sample2workflow2{$sample}.="perl $reformatsnpeffvcf -i $outputfolder/$sample/snpanno/${sample}.filtered-cleaned.snpeff.sift.annotated.vcf -o $outputfolder/$sample/snpanno/${sample}.filtered-cleaned.snpeff.sift.annotated.edited.txt;";

		
		#before filtered
		#--
		$sample2workflow2{$sample}.="$java -jar $snpeff -v ".$tx2ref{$tx}{"genomeversion"}." $outputfolder/$sample/snpanno/${sample}.filtered.sift.annotated.vcf -stats $outputfolder/$sample/snpanno/$sample.filtered.snpEff_summary.html > $outputfolder/$sample/snpanno/${sample}.filtered.snpeff.sift.annotated.vcf 2> $outputfolder/$sample/snpanno/snpeff2.log;";

		$sample2workflow2{$sample}.="mv $outputfolder/$sample/snpanno/snpEff_summary.genes.txt $outputfolder/$sample/snpanno/$sample.filtered.snpEff_summary.genes.txt;";
		
		#reformat snpeff ANN column
		$sample2workflow2{$sample}.="perl $reformatsnpeffvcf -i $outputfolder/$sample/snpanno/${sample}.filtered.snpeff.sift.annotated.vcf -o $outputfolder/$sample/snpanno/${sample}.filtered.snpeff.sift.annotated.edited.txt;";
		
		
		#Remove temporary files
		#------------------
		$sample2workflow2{$sample}.="rm ".join(" ",sort keys %{$tempfiles2rm2{$sample}}).";";
	}
}


########
#Print out commands, for local and server run
########


open(S1,">$scriptfile1") || die "Error writing $scriptfile1. $!";

foreach my $sample (sort keys %sample2workflow) {
	print S1 $sample2workflow{$sample},"\n";
}

close S1;

open(S2,">$scriptfile2") || die "Error writing $scriptfile2. $!";

foreach my $sample (sort keys %sample2workflow2) {
	print S2 $sample2workflow2{$sample},"\n";
}

close S2;


#######
#Run mode
#######

open(LOUT,">$scriptlocalrun") || die "ERROR:can't write to $scriptlocalrun. $!";
open(SOUT,">$scriptclusterrun") || die "ERROR:can't write to $scriptclusterrun. $!";


my @scripts_all=($scriptfile1,$scriptfile2);


#print out command for local and cluster parallel runs
my $jobnumber=0;
my $jobname="dnaseq-process-$timestamp";

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

my $clustercommand="perl $parallel_job -i ".join(",", @scripts_all)." -o $scriptfolder -n ".join(",",@script_names)." --tandem -t $task --ncpus $ncpus --env"; #changed here for none version

if(defined $mem && length($mem)>0) {
	$clustercommand.=" -m $mem";	
}

print SOUT "sh $outputfolder/scripts/parallel-job_submit.sh\n"; #submit step
close SOUT;


system("$clustercommand");
print LOG "$clustercommand;\n\n";



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



