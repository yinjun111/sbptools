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


my $version="0.41";

#v0.2, add system run
#v0.21, add locate_cmd
#V0.3, add --dev
#v0.4, add background for motif finder
#v0.41, versioning

my $usage="

rnaseq-motif
version: $version
Usage: sbptools rnaseq-motif [parameters]

Description: RNA-Seq motif analysis for DE gene promoters


Parameters:

    --gl              Gene list
    --ge              Gene DE result from rnaseq-de

    For --promoter detx option (not implemeted yet)
    --te              Tx DE result from rnaseq-de
	
    --up              Upstream region length (not implemeted yet) [1000]
    --down            Downstream region length (not implemeted yet) [100]

    --promoter        alltx,longesttx,detx [alltx,longesttx]
    --de              both,up,down [both,up,down]

    --out|-o          Output folder
	
    --tx|-t           Transcriptome
                        Current support Human.B38.Ensembl84, Mouse.B38.Ensembl84

    --runmode|-r      Where to run the scripts, local, system, server or none [none]
                                  none, only prints scripts
                                  local, run the scripts using parallel in local workstation(Falco)	
                                  system, run the scripts without paralleling (Falco)
                                  server, run the scripts in cluster (Firefly)
								  
    --jobs|-j         Number of jobs to be paralleled. By default 5 jobs. [5]
	
	
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

my $genelist;
my $genede;
my $up=1000;
my $down=100;
my $promoter="alltx,longesttx";
my $de="both,up,down";
my $tx;
my $outputfolder;

my $runmode="none";
my $jobs=5;
my $dev=0;
my $verbose;

GetOptions(
	"gl|l=s" => \$genelist,
	"ge|e=s" => \$genede,

	"up=s" => \$up,	
	"down=s" => \$down,	
	
	"promoter=s" => \$promoter,	
	"de=s" => \$de,	
	
	"output|o=s" => \$outputfolder,

	"tx|t=s" => \$tx,	
	"runmode|r=s" => \$runmode,	
	"jobs|j=s" => \$jobs,	
	"verbose|v" => \$verbose,
	
	"dev" => \$dev,
);


my %promoters=map {$_,1} split(",",$promoter);
my %des=map {$_,1} split(",",$de);


########
#Prerequisites
########

#my $sbptools=locate_cmd("sbptools","/usr/bin/sbptools");
#my $motiffinder="$sbptools motif-finder"; #improve compatibility

my $sbptoolsfolder="/apps/sbptools/";

#Dev version
if($dev) {
	$sbptoolsfolder="/home/jyin/Projects/Pipeline/sbptools/";
}
else {
	#the tools called will be within the same folder of the script
	$sbptoolsfolder=get_parent_folder(abs_path(dirname($0)));
}


my $motiffinder="perl $sbptoolsfolder/motif-finder/motif-finder_caller.pl"; #improve compatibility
my $parallel_job="perl $sbptoolsfolder/parallel-job/parallel-job_caller.pl"; 
my $rnaseq_motif_annotate="perl $sbptoolsfolder/rnaseq-motif/rnaseq-motif_annotate.pl";


########
#Code begins
########

#Create folders

if(!-e $outputfolder) {
	mkdir($outputfolder);
}

$outputfolder = abs_path($outputfolder);


my $logfile="$outputfolder/rnaseq-motif_run.log";



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
		"txanno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_tx_annocombo.txt",
		"promoter_alltx"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_promoter_1kup100down_alltxs_nr.bed",
		"promoter_longesttx"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_promoter_1kup100down_longesttxs_nr.bed",
		"promoter_alltxanno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_promoter_1kup100down_alltxs_nr_intersect_homer.Knownmotifs_anno.txt"
		},
	
	"Mouse.B38.Ensembl84"=>{ 
		"star"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl84_STAR",
		"chrsize"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl84_STAR/chrNameLength.txt",
		"fasta"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.dna.primary_assembly_ucsc.fa",
		"gtf"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc.gtf",
		"homeranno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_homeranno.txt",
		"geneanno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_gene_annocombo_rev.txt",
		"txanno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_tx_anno.txt",
		"promoter_alltx"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_promoter_1kup100down_alltxs_nr.bed",		
		"promoter_longesttx"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_promoter_1kup100down_longesttxs_nr.bed",
		"promoter_alltxanno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_promoter_1kup100down_alltxs_nr_intersect_homer.Knownmotifs_anno.txt"
		},
		
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

my $scriptfile1="$scriptfolder/rnaseq-motif_run1.sh";



########
#Process
########

print STDERR "\nsbptools rnaseq-motif $version running ...\n\n" if $verbose;
print LOG "\nsbptools rnaseq-motif $version running ...\n\n";


########
#annotation of gene to tx
########

my %gene2tx;
my %tx2gene;

open(IN,$tx2ref{$tx}{"geneanno"}) || die $!;
while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^Gene/;
	my @array=split/\t/;
	
	foreach my $tx (split(",",$array[4])) {
		$gene2tx{$array[0]}{$array[1]."|".$tx}++; #symbol+tx
		$tx2gene{$array[1]."|".$tx}=$array[0];
	}
}
close IN;

########
#Gene to promoter
########


my %tx2promoter;
my %gene2promoter_alltx;
my %gene2promoter_longesttx;

#all txs promoter
open(IN,$tx2ref{$tx}{"promoter_alltx"}) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	
	foreach my $tx (split(",",$array[3])) {
		$tx2promoter{$tx}=$_;
		if(defined $tx2gene{$tx}) {
			$gene2promoter_alltx{$tx2gene{$tx}}{$_}++;
		}
		else {
			print STDERR "ERROR:$tx not defined.\n";
		}
	}
}
close IN;


#longest tx promoter
open(IN,$tx2ref{$tx}{"promoter_longesttx"}) || die $!;
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	
	foreach my $tx (split(",",$array[3])) {
		$tx2promoter{$tx}=$_;
		if(defined $tx2gene{$tx}) {
			$gene2promoter_longesttx{$tx2gene{$tx}}=$_;
		}
		else {
			print STDERR "ERROR:$tx not defined.\n";
		}
	}
}
close IN;

#de tx promoter ....



########
#Process Gene list --gl, Gene DE --ge
########


my %bedfiles;

if(defined $genelist && length($genelist)>0) {
	my %genes_list;

	open(IN,$genelist) || die $!;
	while(<IN>) {
		tr/\r\n//d;
		$genes_list{$_}++;
	}
	close IN;
	
	print STDERR scalar(keys %genes_list)," genes identified from $genelist.\n";
	
	
	#print out bed files 
	
	my $outfile;
	
	if(defined $promoters{"alltx"}) {
		#$outfile=basename($genelist);
		#$outfile=~s/\.\w+$/_alltx_promoter.bed/;

		$outfile=basename($outputfolder)."_alltx_promoter.bed";

		open(OUT,">$outputfolder/$outfile") || die $!;
		foreach my $gene (sort keys %genes_list) {
			if(defined $gene2promoter_alltx{$gene}) {
				print OUT join("\n",sort keys %{$gene2promoter_alltx{$gene}}),"\n";
			}
			else {
				print STDERR "ERROR:Gene $gene not defined.\n";
			}
		}
		close OUT;
		
		$bedfiles{"$outputfolder/$outfile"}=$tx2ref{$tx}{"promoter_alltx"};
	}
	
	if(defined $promoters{"longesttx"}) {
		#$outfile=basename($genelist);
		#$outfile=~s/\.\w+$/_longesttx_promoter.bed/;
		
		$outfile=basename($outputfolder)."_longesttx_promoter.bed";

		open(OUT,">$outputfolder/$outfile") || die $!;
		foreach my $gene (sort keys %genes_list) {
			if(defined $gene2promoter_longesttx{$gene}) {
				print OUT $gene2promoter_longesttx{$gene},"\n";
			}
			else {
				print STDERR "ERROR:Gene $gene not defined.\n";
			}
		}
		close OUT;
		
		$bedfiles{"$outputfolder/$outfile"}=$tx2ref{$tx}{"promoter_longesttx"};
	}	
}
elsif(defined $genede && length($genede)>0) {
	my %genes_de;

	open(IN,$genede) || die $!;
	while(<IN>) {
		tr/\r\n//d;
		next unless $_=~/ENS/; #no title row
		my @array=split/\t/;
		
		if($array[5] eq "1") {
			$genes_de{"up"}{$array[0]}++;
			$genes_de{"both"}{$array[0]}++;
		}
		elsif($array[5] eq "-1") {
			$genes_de{"down"}{$array[0]}++;
			$genes_de{"both"}{$array[0]}++;
		}
	}
	close IN;
	
	print STDERR scalar(keys %{$genes_de{"up"}})," up-regulated genes identified from $genede.\n";
	print STDERR scalar(keys %{$genes_de{"down"}})," down-regulated genes identified from $genede.\n";
	print STDERR scalar(keys %{$genes_de{"both"}})," up/down-regulated genes identified from $genede.\n";
	
	print LOG scalar(keys %{$genes_de{"up"}})," up-regulated genes identified from $genede.\n";
	print LOG scalar(keys %{$genes_de{"down"}})," down-regulated genes identified from $genede.\n";
	print LOG scalar(keys %{$genes_de{"both"}})," up/down-regulated genes identified from $genede.\n";	
	
	#implement DE tx here .....
	
	
	
	
	
	#print out bed files 
	
	my $outfile;
	
	
	foreach my $dechoice (sort keys %des) {
		if(defined $genes_de{$dechoice}) {
	
			if(defined $promoters{"alltx"}) {
				#$outfile=basename($genede);
				#$outfile=~s/\.\w+$/_$dechoice\_alltx_promoter.bed/;
				$outfile=basename($outputfolder)."_$dechoice\_alltx_promoter.bed";
				
				open(OUT,">$outputfolder/$outfile") || die $!;
				foreach my $gene (sort keys %{$genes_de{$dechoice}}) {
					if(defined $gene2promoter_alltx{$gene}) {
						print OUT join("\n",sort keys %{$gene2promoter_alltx{$gene}}),"\n";
					}
					else {
						print STDERR "ERROR:Gene $gene not defined.\n";
					}
				}
				close OUT;
				
				$bedfiles{"$outputfolder/$outfile"}=$tx2ref{$tx}{"promoter_alltx"};
			}
			
			if(defined $promoters{"longesttx"}) {
				#$outfile=basename($genede);
				#$outfile=~s/\.\w+$/_$dechoice\_longesttx_promoter.bed/;
				
				$outfile=basename($outputfolder)."_$dechoice\_longesttx_promoter.bed";

				open(OUT,">$outputfolder/$outfile") || die $!;
				foreach my $gene (sort keys %{$genes_de{$dechoice}}) {
					if(defined $gene2promoter_longesttx{$gene}) {
						print OUT $gene2promoter_longesttx{$gene},"\n";
					}
					else {
						print STDERR "ERROR:Gene $gene not defined.\n";
					}
				}
				close OUT;
				
				$bedfiles{"$outputfolder/$outfile"}=$tx2ref{$tx}{"promoter_longesttx"};
			}
		}
	}
}


########
#Print out motif-finder script for different bed files
########



open(OUT,">$scriptfile1") || die $!;

foreach my $bedfile (sort keys %bedfiles) {
	my $outfoldername=basename($bedfile);
	
	$outfoldername=~s/\.\w+$//;
	
	if(!-e "$outputfolder/$outfoldername") {
		mkdir("$outputfolder/$outfoldername");
		mkdir("$outputfolder/$outfoldername/findmotifsbed");
	}


	#print OUT "$motiffinder -i $bedfile -o $outputfolder/$outfoldername --tx $tx\n";	
	#v0.4 update, add all promoters to bg
	#v0.4, get annotation from pre-computed file
	print OUT "$rnaseq_motif_annotate -i $bedfile -a ",$tx2ref{$tx}{"promoter_alltxanno"}," -o $outputfolder/$outfoldername/findmotifsbed/$outfoldername\_motif_homer_anno.txt;";
	print OUT "$motiffinder -i $bedfile -b ",$bedfiles{$bedfile}," -o $outputfolder/$outfoldername --tx $tx --annobed F\n";
}

close OUT;





#######
#Run mode
#######

my $jobnumber=0;
my $jobname="rnaseq-motif-$timestamp";

if($jobs eq "auto") {
	$jobnumber=0;
}
else {
	$jobnumber=$jobs;
}

my $localcommand="screen -S $jobname -dm bash -c \"source ~/.bashrc;cat $scriptfile1 | parallel -j $jobnumber;\"";


if($runmode eq "none") {
	print STDERR "\nTo run locally, in shell type: $localcommand\n\n";
	print LOG "\nTo run locally, in shell type: $localcommand\n\n";
	
	print STDERR "\nTo run in cluster, in shell type:$parallel_job -i $scriptfile1 -o $outputfolder/scripts/ -r\n\n";
	print LOG "\nTo run in cluster, in shell type:$parallel_job -i $scriptfile1 -o $outputfolder/scripts/ -r\n\n";
}
elsif($runmode eq "local") {
	#local mode
	
	#need to replace with "sbptools queuejob" later

	system($localcommand);
	print LOG "$localcommand;\n\n";

	print STDERR "Starting local paralleled processing using $jobnumber tasks. To monitor process, use \"screen -r $jobname\".\n\n";
	print LOG "Starting local paralleled processing using $jobnumber tasks. To monitor process, use \"screen -r $jobname\".\n\n";
	
}
if($runmode eq "system") {
	print STDERR "\nRunning locally without parallel. sh $scriptfile1.\n\n";
	print LOG "\nRunning locally without parallel. sh $scriptfile1.\n\n";
	
	system("sh $scriptfile1");
}
elsif($runmode eq "cluster") {
	#cluster mode	
	print STDERR "\nStart running in cluster using:$parallel_job -i $scriptfile1 -o $outputfolder/scripts/ -r\n\n";
	system("$parallel_job -i $scriptfile1 -o $outputfolder/scripts/ -r");
}


close LOG;


########
#Functions
########

sub locate_cmd {
	my ($cmdname,$defaultlocation)=@_;
	my $cmdlocation;
	
	if(-e $defaultlocation) {
		$cmdlocation=$defaultlocation;
	}
	else {
		use File::Which qw(which where); 
		$cmdlocation = which $cmdname;
	}
	
	return $cmdlocation;
}


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



sub get_parent_folder {
	my $dir=shift @_;
	
	if($dir=~/^(.+\/)[^\/]+\/?/) {
		return $1;
	}
}