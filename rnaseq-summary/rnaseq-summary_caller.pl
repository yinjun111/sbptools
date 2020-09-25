#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);
use List::Util qw(sum);

#CutAdapt+FASTQC+RSEM+STAR




########
#Interface
########


my $version="0.62";

#v0.1b, changed DE match pattern
#v0.1c, add first line recognition in DE results

#v0.4, to be consistent with other rnaseq scripts
#v0.5, add metascape and gsea support
#v0.51, pass --dev to other programs
#v0.52, new gsea-gen params
#v0.6, gsea-gen/summary and rnaseq-motif/summary will be ran by parallel-job 
#v0.61, versioning
#v0.62, DE folder signature changed
#v0.63, Add note to summary file

my $usage="

rnaseq-summary
version: $version
Usage: sbptools rnaseq-summary [parameters]

Description: Summarize rnaseq-de results and recalculate significance if needed
1. signficant calls from DE results
2. reformat DE results using union gene lists
3. report group avg TPM, and FPKM
4. recalculate significance using new fc or new q
5. merged excel file for the project
6. Files ready for IPA/Metascape/GSEA/rnaseq-motif analyses


Parameters:

    --in|-i           Input rnaseq-de folder(s)
    --output|-o       Output folder

    --config|-c       Configuration file match the samples in the rnaseq-merge folder
                           first column as sample name.
    --group|-g        Group name in config file for avg calculation

    --rnaseq-merge    rnaseq-merge folder to retrieve TPM/FPKM data

    --gi              use selected gene list to summarize results	
    --ti              use selected tx list to summarize results

    --fccutoff        Log2 FC cutoff, optional
    --qcutoff         Corrected P cutoff, optional
	
    --tx|-t           Transcriptome
                        Current support Human.B38.Ensembl84, Mouse.B38.Ensembl84

    --run_rnaseq-motif  Whether to run sbptools rnaseq-motif, only generate script by default [none]
                        use \"cluster\" to run in Firefly

    --run_gsea-gen      Whether to run sbptools gsea-gen, run by default [cluster]
                        use \"none\" to turn off

    --gseadbs           Default dbs to run are [h.all.v7.1,c5.bp.v7.1]
                        Other popular dbs include c2.cp.kegg.v7.1,c3.tft.v7.1,c5.cc.v7.1,c5.mf.v7.1
	
";

#    --verbose|-v      Verbose
#    --runmode|-r      Where to run the scripts, local, cluster or none [local]
#    --jobs|-j         Number of jobs to be paralleled. By default 5 jobs. [5]


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
my $rnaseqmerge;
my $fccutoff;
my $qcutoff;
my $geneinput;
my $txinput;
my $gseadbs="h.all.v7.1,c5.bp.v7.1";
my $runrnaseqmotif="none";
my $rungseagen="cluster";
my $verbose=1;
my $jobs=5;
my $tx;
my $runmode="local";

my $dev=0; #developmental version

GetOptions(
	"in|i=s" => \$inputfolders,
	"output|o=s" => \$outputfolder,
	"config|c=s" => \$configfile,

	"group|g=s" => \$group,
	"rnaseq-merge=s" => \$rnaseqmerge,
	
	"fccutoff=s" => \$fccutoff,
	"qcutoff=s" => \$qcutoff,
	
	"tx|t=s" => \$tx,	
	"runmode|r=s" => \$runmode,		
	"jobs|j=s" => \$jobs,
	"verbose|v" => \$verbose,
	
	"run_rnaseq-motif=s"=>\$runrnaseqmotif,
	"run_gsea-gen=s"=>\$rungseagen,
	"gseadbs=s"=>\$gseadbs,
	
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
	$sbptoolsfolder=get_parent_folder(abs_path(dirname($0)));
}


my $mergefiles="$sbptoolsfolder/mergefiles/mergefiles_caller.pl";
my $text2excel="$sbptoolsfolder/text2excel/text2excel.pl";
my $metascape_gen="perl $sbptoolsfolder/metascape-gen/metascape-gen_caller.pl";
my $gsea_gen="perl $sbptoolsfolder/gsea-gen/gsea-gen_caller.pl";
my $gsea_gen_summary="perl $sbptoolsfolder/gsea-gen-summary/gsea-gen-summary.pl".scalar(add_dev($dev));
my $rnaseq_motif="perl $sbptoolsfolder/rnaseq-motif/rnaseq-motif_caller.pl".scalar(add_dev($dev));
my $rnaseq_motif_summary="perl $sbptoolsfolder/rnaseq-motif-summary/rnaseq-motif-summary.pl";
my $parallel_job="perl $sbptoolsfolder/parallel-job/parallel-job_caller.pl"; 

my $rnaseqsummarynote="$sbptoolsfolder/rnaseq-summary/rnaseq-summary_note.txt";


my $zip=find_program("/usr/local/bin/zip");


########
#Code begins
########

#Create folders

if(!-e $outputfolder) {
	mkdir($outputfolder);
}

$outputfolder = abs_path($outputfolder);
my $outputfoldername = basename($outputfolder);


#folders for gs analysis
if(!-e "$outputfolder/forIPA") {
	mkdir("$outputfolder/forIPA");
}

if(!-e "$outputfolder/forMetascape") {
	mkdir("$outputfolder/forMetascape");
}

my $gsea_gen_run="$outputfolder/forGSEA/gsea-gen_cluster_run.sh";
if(!-e "$outputfolder/forGSEA") {
	mkdir("$outputfolder/forGSEA");
}


if(!-e "$outputfolder/for_rnaseq-motif/") {
	mkdir("$outputfolder/for_rnaseq-motif/");
}


my $scriptfolder="$outputfolder/scripts";

if(!-e $scriptfolder) {
	mkdir($scriptfolder);
}

my $scriptfile1="$scriptfolder/rnaseq-summary_run1.sh";
my $scriptfile2="$scriptfolder/rnaseq-summary_run2.sh";

my $logfile="$outputfolder/rnaseq-summary_run.log";

#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();

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
		#"geneannogsea"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_gene_annocombo_rev.txt",
		"txanno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_tx_annocombo.txt"},
	"Mouse.B38.Ensembl84"=>{ 
		"star"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl84_STAR",
		"chrsize"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl84_STAR/chrNameLength.txt",
		"fasta"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.dna.primary_assembly_ucsc.fa",
		"gtf"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc.gtf",
		"homeranno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_homeranno.txt",
		"geneanno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_gene_annocombo_rev.txt",
		#"geneannogsea"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_tohumansymmbol_one2one.txt",		
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


#cutoffs to report
##current cutoff
##FC2, bhp005
##FC1.5, bhp005
##FC4, bhp005
##FC2, bhp001

my %report_cutoffs=(
	"FC2BHP005"=> [1,0.05],
	"FC15BHP005"=> [0.585,0.05],
	"FC4BHP005"=> [2,0.05],
	"FC2BHP001"=> [1,0.01],
	"FC4BHP001"=> [2,0.01],
);


#Files for merging
#Examples
#JY_315_Isl.genes.results
#JY_315_Isl.isoforms.results


########
#Process
########

print STDERR "\nsbptools rnaseq-summary $version running ...\n\n" if $verbose;
print LOG "\nsbptools rnaseq-summary $version running ...\n\n";

#open input folders to find comparisons

my %folder2genede;
my %folder2txde;
my %folder2dir;
my %comparisons;

foreach my $inputfolder (split(",",$inputfolders)) {
	#print STDERR $inputfolder,"#1\n";
	my @folders=glob("$inputfolder/*");
	foreach my $folder (@folders) {
		if(-d $folder) {
			#print STDERR $folder,"#2\n";
			#this signature file is changed in v0.6
			if(-e "$folder/gene_de_test_run.log") {
				#rnaseq-de folder
				my $foldername=basename($folder);
				
				#foldername needs to be unique
				if(defined $folder2dir{$foldername}) {
					print STDERR "ERROR:$foldername has been used twice:\n",abs_path($folder),"\n",$folder2dir{$foldername},"\n\n";
					print LOG "ERROR:$foldername has been used twice:\n",abs_path($folder),"\n",$folder2dir{$foldername},"\n\n";
					exit;
				}

				$comparisons{$foldername}++;
				
				$folder2dir{$foldername}=abs_path($folder);
				my @files=map {basename($_)} glob("$folder/*");
				
				foreach my $file (@files) {
					if($file=~/gene.results.merged.count.[\.\w\d]+FC[\.\d]+.\w+.P[\.\d]+.txt/) {
						$folder2genede{$foldername}=$file;
					}
					if($file=~/tx.results.merged.count.[\.\w\d]+FC[\.\d]+.\w+.P[\.\d]+.txt/) {
						$folder2txde{$foldername}=$file;
					}
				}
			}
		}
	}
}

#error message

if(keys %folder2genede ==0) {
	print STDERR "\nERROR:No rnaseq-de folder was detected in $inputfolders\n\n";
	print LOG "\nERROR:No rnaseq-de folder was detected in $inputfolders\n\n";
	exit;
}


#print out what were found

print STDERR "Printing files identified from input folder(s).\n\n" if $verbose;
print LOG "Printing files identified from input folder(s).\n\n";


foreach my $folder (sort keys %folder2genede) {
	print STDERR $folder,"\t",$folder2genede{$folder},"\n" if $verbose;
	print LOG $folder,"\t",$folder2genede{$folder},"\n";
}

foreach my $folder (sort keys %folder2txde) {
	print STDERR $folder,"\t",$folder2txde{$folder},"\n" if $verbose;
	print LOG $folder,"\t",$folder2txde{$folder},"\n";
}

#####
#read input gene/tx list
#####

my %genes;

if(defined $geneinput && length($geneinput)>0) {
	print STDERR "\nReading predefined gene list from --gi $geneinput.\n" if $verbose;
	print LOG "\nReading predefined gene list from --gi $geneinput.\n";
	
	open(IN,$geneinput) || die "$geneinput not readable.$!\n";
	
	while(<IN>) {
		tr/\r\n//d;
		next if $_=~/^(Gene|Feature)/; 
		$genes{$_}++;
	}
	close IN;
	
	print STDERR scalar(keys %genes)," genes found from $geneinput.\n" if $verbose;
	print LOG scalar(keys %genes)," genes found from $geneinput.\n" 

}


my %txs;

if(defined $txinput && length($txinput)>0) {
	print STDERR "\nReading predefined tx list from --gi $txinput.\n" if $verbose;
	print LOG "\nReading predefined tx list from --gi $txinput.\n";
	
	open(IN,$txinput) || die "$txinput not readable.$!\n";
	
	while(<IN>) {
		tr/\r\n//d;
		next if $_=~/^tx$/i; 
		$txs{$_}++;
	}
	close IN;
	
	print STDERR scalar(keys %txs)," txs found from $txinput.\n" if $verbose;
	print LOG scalar(keys %txs)," txs found from $txinput.\n" 

}

#####
#Config file
#####

#read config info

my %sample2group;
my $firstline=0;
my @covtitle;
my $groupcolnum=0;


#push @outfiles,"$cov";
#push @outfilenames,"SampleIndex";

print STDERR "\nReading config file $configfile\n" if $verbose;
print LOG "\nReading config file $configfile\n";

open(IN,$configfile) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	if($firstline==0) {
		@covtitle=split(/\t/,$_);
		for(my $num=0;$num<@covtitle;$num++) {
			if($covtitle[$num] eq $group) {
				$groupcolnum=$num;
				last;
			}
		}
		$firstline++;
	}
	else {
		$sample2group{$array[0]}=$array[$groupcolnum];
	}
}
close IN;


#####
#Work on gene DE files
#####

#read file
my %folder2genesig;
my %folder2geneinfo;
my %folder2genetitle;

my %cutoff_summary;

foreach my $folder (sort keys %folder2genede) {

	print STDERR "Processing ",$folder2dir{$folder}."/".$folder2genede{$folder},"\n" if $verbose;
	print LOG "Processing ",$folder2dir{$folder}."/".$folder2genede{$folder},"\n";
	
	my $foripafile="$folder\_GeneDE.foripa.txt";
	
	#print out file ready for IPA
	system("cut -f 1,2,5 ".$folder2dir{$folder}."/".$folder2genede{$folder}."> $outputfolder/forIPA/$foripafile");
	
	open(IN,$folder2dir{$folder}."/".$folder2genede{$folder}) || die $!;
	
	my $linenum=0;
	
	while(<IN>) {
		tr/\r\n//d;
		my @array=split/\t/;
		
		if($linenum==0) {

			if(defined $fccutoff && length($fccutoff)>0 && defined $qcutoff && length($qcutoff)>0) {
				$array[5]="Significance: Log2FC $fccutoff BHP $qcutoff";
			}
			elsif((defined $fccutoff && length($fccutoff)>0) || (defined $qcutoff && length($qcutoff)>0)) {
				print STDERR "ERROR: --fccutoff and --qcutoff have to be both defined.\n\n";
				print LOG "ERROR: --fccutoff and --qcutoff have to be both defined.\n\n";
				exit;
			}
			
			$folder2genetitle{$folder}=join("\t",@array);
		}
		else {
			
			#record genes
			unless(defined $geneinput && length($geneinput)>0) {
				$genes{$array[0]}++;
			}
			
			#new fc and q
			if($array[4] ne "NA" && $array[4] ne " ") {
				if(defined $fccutoff && length($fccutoff)>0 && defined $qcutoff && length($qcutoff)>0) {
					if($array[1] >= $fccutoff && $array[4] < $qcutoff) {
						$array[5]=1;
					}
					elsif($array[1] <= -$fccutoff && $array[4] < $qcutoff) {
						$array[5]=-1;
					}
					else{
						$array[5]=0;
					}
				}
				
				#summary for different cutoffs
				foreach my $cutoff (sort keys %report_cutoffs) {
					if($array[1] >= $report_cutoffs{$cutoff}[0] && $array[4] < $report_cutoffs{$cutoff}[1]) {
						$cutoff_summary{$cutoff}{$folder}{1}++;
					}
					elsif($array[1] <= -$report_cutoffs{$cutoff}[0] && $array[4] < $report_cutoffs{$cutoff}[1]) {
						$cutoff_summary{$cutoff}{$folder}{-1}++;
					}
				}
				
			}
			
			$cutoff_summary{"current"}{$folder}{$array[5]}++;
			$folder2genesig{$folder}{$array[0]}=$array[5];
			$folder2geneinfo{$folder}{$array[0]}=join("\t",@array);
		}
		
		$linenum++;
	}
	close IN;
	
	
}


#Summary of No. of DE Genes

foreach my $cutoff (sort keys %report_cutoffs,"current") {
	my $cutoffsumfile="$outputfolder/rnaseq-summary_GeneDESigs_summary_$cutoff.txt";
	open(OUT,">$cutoffsumfile") || die $!;
	print OUT "Comparisons_$cutoff\tNo. of Up-Regulated Genes\tNo. of Down-Regulated Genes\n";
	foreach my $folder (sort keys %folder2genede) {
		print OUT $folder,"\t";
		if(defined $cutoff_summary{$cutoff}{$folder}{1}) {
			print OUT $cutoff_summary{$cutoff}{$folder}{1},"\t";
		}
		else {
			print OUT "0\t";
		}
		if(defined $cutoff_summary{$cutoff}{$folder}{-1}) {
			print OUT $cutoff_summary{$cutoff}{$folder}{-1},"\n";
		}
		else {
			print OUT "0\n";
		}
	}
	close OUT;
}


#print OUT used genes

open(OUT,">$outputfolder/genes_sel.txt") || die $!;
print OUT join("\n","Gene",sort keys %genes),"\n";
close OUT;


#print OUT sig file
open(OUT,">$outputfolder/rnaseq-summary_GeneDESigs.txt") || die $!;
print OUT "Gene\t",join("\t",sort keys %folder2genede) ,"\n";

foreach my $gene (sort keys %genes) {
	print OUT $gene,"\t";
	my @marks;
	
	foreach my $folder (sort keys %folder2genede) {
		if(defined $folder2genesig{$folder}{$gene}) {
			push @marks,$folder2genesig{$folder}{$gene};
		}
		else {
			push @marks," ";
		}
	}
	
	print OUT join("\t",@marks),"\n";
}
close OUT;


#reformat gene de
my @defiles;
foreach my $folder (sort keys %folder2genede) {
	
	push @defiles,"$outputfolder/$folder\_GeneDEreformated.txt";
	open(OUT,">$outputfolder/$folder\_GeneDEreformated.txt") || die $!;
	print OUT "Gene",$folder2genetitle{$folder},"\n";
	
	my @title=split("\t",$folder2genetitle{$folder});
	my $colnum=@title;
	
	foreach my $gene (sort keys %genes) {
		if(defined $folder2geneinfo{$folder}{$gene}) {
			print OUT $folder2geneinfo{$folder}{$gene},"\n";
		}
		else {
			print OUT $gene,"\t",join("\t",(" ") x ($colnum-1)),"\n";
		}
	}
	close OUT;
}


#merge gene de results
system("$mergefiles -m $outputfolder/genes_sel.txt -i ".join(",",@defiles)." -o $outputfolder/rnaseq-summary_GeneDEMerged.txt -n T");

print STDERR "Calculating group average for gene.\n" if $verbose;
print LOG "Calculating group average for gene.\n";


#read FPKM/TPM,cal group
foreach my $file ("gene.results.merged.fpkm.txt","gene.results.merged.tpm.txt") {

	my $outfile0=$file;
	$outfile0=~s/\.txt/.sel.txt/;
	
	#copy files
	system("cp $rnaseqmerge/$file $outputfolder/");
	
	system("$mergefiles -m $outputfolder/genes_sel.txt -i $rnaseqmerge/$file -o $outputfolder/$outfile0 -n T"); #edited here to remove extra rownames
	
	my $outfile1=$file;
	my $outfile2=$file;
	my @groupnames;
	
	$outfile1=~s/\.txt/.groupavg.txt/;
	$outfile2=~s/\.txt/.groupavg.sel.txt/;
	
	open(IN,"$rnaseqmerge/$file") || die $!;
	open(OUT1,">$outputfolder/$outfile1") || die $!;
	open(OUT2,">$outputfolder/$outfile2") || die $!;

	while(<IN>) {
		tr/\r\n//d;
		my @array=split/\t/;
		
		if($_=~/^Gene\t/) {
			#first line
			foreach my $sample (@array[1..$#array]) {
				push @groupnames,$sample2group{$sample};
			}
			
			print OUT1 "Gene\t",join("\t",sort(uniq(@groupnames))),"\n";
			print OUT2 "Gene\t",join("\t",sort(uniq(@groupnames))),"\n";
		}
		else {
			my @groupavg=cal_group_avg([@array[1..$#array]],\@groupnames);
			
			print OUT1 $array[0],"\t",join("\t",@groupavg),"\n";
			
			if(defined $genes{$array[0]}) {
				print OUT2 $array[0],"\t",join("\t",@groupavg),"\n";
			}	
		}
	}
	close IN;
	close OUT1;
	close OUT2;
}


####
#Print OUT gene annotation
####

print STDERR "Copy gene annotation.\n" if $verbose;
print LOG "Copy gene annotation.\n";

system("cp ".$tx2ref{$tx}{"geneanno"}." $outputfolder/geneanno.txt;");
system("cp ".$tx2ref{$tx}{"txanno"}." $outputfolder/txanno.txt;");

system("$mergefiles -m $outputfolder/genes_sel.txt -i ".$tx2ref{$tx}{"geneanno"}." -o $outputfolder/geneanno_sel.txt -n T");


#####
#GS Analysis
#####

#IPA
#-----------
print STDERR "Files ready for IPA analysis are in: $outputfolder/forIPA\n" if $verbose;
print LOG "Files ready for IPA analysis are in: $outputfolder/forIPA\n";

#Metascape
#-----------
print STDERR "Generate files for Metascape analysis.\n" if $verbose;
print LOG "Generate files for Metascape analysis.\n";
print LOG "$metascape_gen -i $outputfolder/rnaseq-summary_GeneDESigs.txt -o $outputfolder/forMetascape/$outputfoldername\_forMetascape.txt\n";

system("$metascape_gen -i $outputfolder/rnaseq-summary_GeneDESigs.txt -o $outputfolder/forMetascape/$outputfoldername\_forMetascape.txt");

print STDERR "Files ready for Metascape analysis are in: $outputfolder/forMetascape/\n" if $verbose;
print LOG "Files ready for Metascape analysis are in: $outputfolder/forMetascape/\n";

#GSEA
#-----------
print STDERR "Generate files for GSEA analysis.\n" if $verbose;
print LOG "Generate files for GSEA analysis.\n";

my $gsearunnum=0;
my @gseascripts;

#Only two tests are performed due to size of the GSEA results
#,c2.cp.kegg.v7.1,c3.tft.v7.1,c5.cc.v7.1,c5.mf.v7.1


#open(OUT,">$gsea_gen_run") || die "ERROR:Can't write into $gsea_gen_run.$!\n\n";

foreach my $comparison (sort keys %comparisons) {
	if($comparison=~/(.+)_vs_(.+)/) {
		print STDERR "Comparison:$comparison recognized. Perform GSEA analysis.\n";
		print LOG "Comparison:$comparison recognized. Perform GSEA analysis.\n";
		
		
		foreach my $gseadb (split(",",$gseadbs)) {
			my $gseacmd="$gsea_gen --tx $tx -e $outputfolder/gene.results.merged.tpm.sel.txt -s ".abs_path($configfile)." -n $group -o $outputfolder/forGSEA/$comparison\-$gseadb -d $gseadb -c $1\_versus_$2";
			
			#print OUT "$gseacmd -r cluster\n";		
		
			
			#just generate script and run later with gsea-summary
			#if($rungseagen eq "cluster") {
			#	system("$gseacmd -r cluster");
			#	print LOG "$gseacmd -r cluster\n";
			#}
			#else {
				system("$gseacmd -r none");
				print LOG "$gseacmd -r none\n";
			#}
			
			push @gseascripts,"$outputfolder/forGSEA/$comparison\-$gseadb/gsea-gen_run.sh";
			
			$gsearunnum++;
		}
	}
}
#close OUT;

if($gsearunnum==0) {
	#failed to recognize any comparison, then just generate cls and gct
	system("$gsea_gen --tx $tx -e $outputfolder/gene.results.merged.tpm.sel.txt -s ".abs_path($configfile)." -n $group -o $outputfolder/forGSEA/GSEARaw/ -r none");
	print LOG "$gsea_gen --tx $tx -e $outputfolder/gene.results.merged.tpm.sel.txt -s ".abs_path($configfile)." -n $group -o $outputfolder/forGSEA/GSEARaw/ -r none\n";
}

#merge all GSEA scripts together
system("cat ".join(" ",@gseascripts)." > $outputfolder/forGSEA/run_gsea_1.sh");

#script for gsea-gen-summary and zip the folder
open(OUT,">$outputfolder/forGSEA/run_gsea_2.sh") || die $!;
print OUT "$gsea_gen_summary -i $outputfolder/forGSEA/ -o $outputfolder/forGSEA/GSEASummary;cd $outputfolder/;$zip -r GSEA.zip forGSEA/\n";
close OUT;

my $gseaclustercommand="$parallel_job -i $outputfolder/forGSEA/run_gsea_1.sh,$outputfolder/forGSEA/run_gsea_2.sh -o $outputfolder/forGSEA/ -n run_gsea_1,run_gsea_2 --tandem -r";

print STDERR "$gseaclustercommand\n";
print LOG "$gseaclustercommand\n";
	
if($rungseagen eq "cluster" && $gsearunnum>0) {
	system($gseaclustercommand);
}

print STDERR "Files ready for GSEA analysis are in: $outputfolder/forGSEA/\n" if $verbose;
print LOG "Files ready for GSEA analysis are in: $outputfolder/forGSEA/\n";

#rnaseq-motif
#-----------
print STDERR "Generate files for rnaseq-motif analysis.\n" if $verbose;
print LOG "Generate files for rnaseq-motif analysis.\n";


foreach my $folder (sort keys %folder2genede) {
	#use new de file
	my $rnaseqmotifcmd="$rnaseq_motif --ge $outputfolder/$folder\_GeneDEreformated.txt --promoter longesttx -o $outputfolder/for_rnaseq-motif/$folder --tx $tx -v";

	#if($runrnaseqmotif eq "cluster") {
	#	system($rnaseqmotifcmd." -r cluster");
	#	print LOG $rnaseqmotifcmd." -r cluster\n";
	#}
	#else {
		system($rnaseqmotifcmd);
		print LOG $rnaseqmotifcmd,"\n";
	#}
}

#merge all running scripts
system("cat $outputfolder/for_rnaseq-motif/*/scripts/*.sh > $outputfolder/for_rnaseq-motif/run_rnaseq-motif_1.sh");

#script for rnaseq-motif-summary
open(OUT,">$outputfolder/for_rnaseq-motif/run_rnaseq-motif_2.sh") || die $!;
print OUT "$rnaseq_motif_summary -i $outputfolder/for_rnaseq-motif/ -o $outputfolder/for_rnaseq-motif/RNASeqMotifSummary\n";
close OUT;

#use parallel-job to run the jobs in tandem
my $rnaseqmotifcommand="$parallel_job -i $outputfolder/for_rnaseq-motif/run_rnaseq-motif_1.sh,$outputfolder/for_rnaseq-motif/run_rnaseq-motif_2.sh -o $outputfolder/for_rnaseq-motif/ -n run_rnaseq-motif_1,run_rnaseq-motif_2 --tandem -r";

print STDERR "$rnaseqmotifcommand\n";
print LOG "$rnaseqmotifcommand\n";


if($runrnaseqmotif eq "cluster" && $gsearunnum>0) {
	system($gseaclustercommand);
}

print STDERR "Files ready for rnaseq-motif analysis are in: $outputfolder/for_rnaseq-motif/\n" if $verbose;
print LOG "Files ready for rnaseq-motif analysis are in: $outputfolder/for_rnaseq-motif/\n";


####
#Add gene annotation
####

#gene de merge
system("$mergefiles -m $outputfolder/rnaseq-summary_GeneDEMerged.txt -i ".$tx2ref{$tx}{"geneanno"}." -o $outputfolder/rnaseq-summary_GeneDEMerged_anno.txt -n T");

#gene sig merge
system("$mergefiles -m $outputfolder/rnaseq-summary_GeneDESigs.txt -i ".$tx2ref{$tx}{"geneanno"}." -o $outputfolder/rnaseq-summary_GeneDESigs_anno.txt -n T");

#gene fpkm
system("$mergefiles -m $outputfolder/gene.results.merged.fpkm.sel.txt -i ".$tx2ref{$tx}{"geneanno"}." -o $outputfolder/gene.results.merged.fpkm.sel_anno.txt -n T");

#gene tpm
system("$mergefiles -m $outputfolder/gene.results.merged.tpm.sel.txt -i ".$tx2ref{$tx}{"geneanno"}." -o $outputfolder/gene.results.merged.tpm.sel_anno.txt -n T");


#gene fpkm group avg
system("$mergefiles -m $outputfolder/gene.results.merged.fpkm.groupavg.sel.txt -i ".$tx2ref{$tx}{"geneanno"}." -o $outputfolder/gene.results.merged.fpkm.groupavg.sel_anno.txt -n T");


#gene tpm group avg
system("$mergefiles -m $outputfolder/gene.results.merged.tpm.groupavg.sel.txt -i ".$tx2ref{$tx}{"geneanno"}." -o $outputfolder/gene.results.merged.tpm.groupavg.sel_anno.txt -n T");



####
#Produce Excel file
####

#merge the previous files together

my $timestamp=substr($now,0,10);

system("$text2excel -i $rnaseqsummarynote,$outputfolder/rnaseq-summary_GeneDEMerged_anno.txt,$outputfolder/rnaseq-summary_GeneDESigs_anno.txt,$outputfolder/gene.results.merged.fpkm.groupavg.sel_anno.txt,$outputfolder/gene.results.merged.tpm.groupavg.sel_anno.txt -n Note,GeneDE,GeneDESigs,FPKM_Group,TPM_Group -o $outputfolder/rnaseq-summary_all_$timestamp.xlsx --theme theme2");

print STDERR "\nResults summary is in $outputfolder/rnaseq-summary_all_$timestamp.xlsx.\n";

####
#work on tx later
####
				


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

sub add_dev {
	my $dev=shift @_;
	
	if($dev) {
		return " --dev";
	}
	else {
		return "";
	}
}


sub get_parent_folder {
	my $dir=shift @_;
	
	if($dir=~/^(.+\/)[^\/]+\/?/) {
		return $1;
	}
}
