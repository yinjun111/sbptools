#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);
use List::Util qw(sum);

#CutAdapt+FASTQC+RSEM+STAR


########
#Prerequisites
########

my $mergefiles="/apps/sbptools/mergefiles/mergefiles_caller.pl";
my $text2excel="perl /apps/sbptools/text2excel/text2excel.pl";



########
#Interface
########


my $version="0.1c";

#v0.1b, changed DE match pattern
#v0.1c, add first line recognition in DE results

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


Parameters:

    --in|-i           Input folder(s)
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

    --runmode|-r      Where to run the scripts, local, server or none [local]
    --jobs|-j         Number of jobs to be paralleled. By default 5 jobs. [5]
	
	
";

#    --verbose|-v      Verbose


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
my $verbose=1;
my $jobs=5;
my $tx;
my $runmode="local";

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
);


########
#Code begins
########

#Create folders

if(!-e $outputfolder) {
	mkdir($outputfolder);
}

$outputfolder = abs_path($outputfolder);


if(!-e "$outputfolder/forIPA") {
	mkdir("$outputfolder/forIPA");
}

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

foreach my $inputfolder (split(",",$inputfolders)) {
	print STDERR $inputfolder,"#1\n";
	my @folders=glob("$inputfolder/*");
	foreach my $folder (@folders) {
		if(-d $folder) {
			print STDERR $folder,"#2\n";
			if(-e "$folder/rnaseq-de_run.log") {
				#rnaseq-de folder
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
			}
			
			$folder2genesig{$folder}{$array[0]}=$array[5];
			$folder2geneinfo{$folder}{$array[0]}=join("\t",@array);
		}
		
		$linenum++;
	}
	close IN;
	
	
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
system("$mergefiles -m $outputfolder/genes_sel.txt -i ".join(",",@defiles)." -o $outputfolder/rnaseq-summary_GeneDEMerged.txt");

print STDERR "Calculating group average for gene.\n" if $verbose;
print LOG "Calculating group average for gene.\n";




#read FPKM/TPM,cal group
foreach my $file ("gene.results.merged.fpkm.txt","gene.results.merged.tpm.txt") {

	my $outfile0=$file;
	$outfile0=~s/\.txt/.sel.txt/;
	
	#copy files
	system("cp $rnaseqmerge/$file $outputfolder/");
	#system("$mergefiles -m $outputfolder/$file -i "..");
	
	system("$mergefiles -m $outputfolder/genes_sel.txt -i $rnaseqmerge/$file -o $outputfolder/$outfile0");
	
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

system("$mergefiles -m $outputfolder/genes_sel.txt -i ".$tx2ref{$tx}{"geneanno"}." -o $outputfolder/geneanno_sel.txt");




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

system("$text2excel -i $outputfolder/rnaseq-summary_GeneDEMerged_anno.txt,$outputfolder/rnaseq-summary_GeneDESigs_anno.txt,$outputfolder/gene.results.merged.fpkm.groupavg.sel_anno.txt,$outputfolder/gene.results.merged.tpm.groupavg.sel_anno.txt -n GeneDE,GeneDESigs,FPKM_Group,TPM_Group -o $outputfolder/rnaseq-summary_all_$timestamp.xlsx --theme theme2");



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

