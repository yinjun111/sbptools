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

#my $convertdetopos="perl /home/jyin/Projects/Pipeline/sbptools/chipseq-de/convert_chipseq_de_to_pos.pl";


#my $homer="/home/jyin/Programs/Homer/bin";
my $homer="/apps/homer/bin/";
my $findmotifsgenome=locate_cmd("findMotifsGenome.pl","$homer/findMotifsGenome.pl";
my $bed2pos=locate_cmd("bed2pos.pl","$homer/bed2pos.pl");

my $bedtools=locate_cmd("bedtools","/apps/bedtools2-2.26.0/bin/bedtools");
my $intersectbed="$bedtools intersect";

my $intersect_multi_bed="/apps/sbptools/motif-finder/intersect_multi_bed.sh";
my $motif_intersect_to_txt="/apps/sbptools/motif-finder/motif_intersect_to_txt.pl";


########
#Interface
########


my $version="0.14";

#v0.11, add mouse motifs
#v0.12, update script directory
#v0.13, add results folder
#v0.14, add locate_cmd

my $usage="

motif-finder
version: $version
Usage: sbptools motif-finder [parameters]

Description: Motif finder for bed files using Homer


Parameters:

    --in|-i           Input bed file
    --output|-o       Output folder
	
    --tx|-t           Transcriptome
                        Current support Human.B38.Ensembl84, Mouse.B38.Ensembl84
						
    --verbose|-v      Verbose
	
	
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

my $inputfile;
my $outputfolder;
my $tx;
my $tfbs="F";
my $verbose=1;
my $pcol=4; #hidden param
#my $runmode="none";


GetOptions(
	"in|i=s" => \$inputfile,
	"output|o=s" => \$outputfolder,

	"tx|t=s" => \$tx,	
	#"runmode|r=s" => \$runmode,		
	"verbose|v=s" => \$verbose,
);


########
#Code begins
########

#Create folders

if(!-e $outputfolder) {
	mkdir($outputfolder);
}

$outputfolder = abs_path($outputfolder);


#sub folders
mkdir("$outputfolder/findmotifsgenome"); #enrichment
mkdir("$outputfolder/findmotifsbed"); #bed themselves

#input file
$inputfile=abs_path($inputfile);
my $inputfile_basename=basename($inputfile);


my $logfile="$outputfolder/motif-finder_run.log";



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
		"txanno"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_tx_annocombo.txt",
		"genepromoter"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_tx_annocombo.txt",
		"homermotif"=>"/data/jyin/Databases/Homer/homer.KnownMotifs.hg38.170917_short.bed",
		"homermotifs"=>"/data/jyin/Databases/Homer/homer.KnownMotifs.hg38.170917/*.bed"},
	"Mouse.B38.Ensembl84"=>{ 
		"star"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl84_STAR",
		"chrsize"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl84_STAR/chrNameLength.txt",
		"fasta"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.dna.primary_assembly_ucsc.fa",
		"gtf"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc.gtf",
		"homeranno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_homeranno.txt",
		"geneanno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_gene_annocombo_rev.txt",
		"txanno"=>"/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.84_ucsc_tx_anno.txt",
		"genepromoter"=>"/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.84_ucsc_tx_annocombo.txt",
		"homermotif"=> "/data/jyin/Databases/Homer/homer.KnownMotifs.mm10.170917_short.bed",
		"homermotifs"=>"/data/jyin/Databases/Homer/homer.KnownMotifs.mm10.170917/*.bed"}
);

my $genomeversion;

if($tx=~/Human.B38/) {
	$genomeversion="hg38";
}
elsif($tx=~/Mouse.B38/) {
	$genomeversion="mm10";
}


#######
#Motif for each peak
#######

my $inputfilename;
if($inputfile_basename=~/(.+)\.bed/) {
	$inputfilename=$1;
}


print STDERR "Identify motifs for $inputfile.\n\n" if $verbose;
print LOG "Identify motifs for $inputfile.\n\n";

print LOG "cut -f 1-4 $inputfile > $outputfolder/$inputfilename\_selected.bed\n\n";
system("cut -f 1-4 $inputfile > $outputfolder/$inputfilename\_selected.bed");

#intersect bed
print LOG "$intersect_multi_bed $outputfolder/$inputfilename\_selected.bed \"".$tx2ref{$tx}{"homermotifs"}."\" $outputfolder/findmotifsbed/$inputfilename\_intersect_homer.txt\n\n";
system("$intersect_multi_bed $outputfolder/$inputfilename\_selected.bed \"".$tx2ref{$tx}{"homermotifs"}."\" $outputfolder/findmotifsbed/$inputfilename\_intersect_homer.txt");

#summarize result
print LOG "$motif_intersect_to_txt $outputfolder/findmotifsbed/$inputfilename\_intersect_homer.txt $outputfolder/findmotifsbed/$inputfilename\_intersect_homer_anno.txt\n\n";
system("$motif_intersect_to_txt $outputfolder/findmotifsbed/$inputfilename\_intersect_homer.txt $outputfolder/findmotifsbed/$inputfilename\_intersect_homer_anno.txt");



########
#bed to pos conversion
########

print STDERR "Convert $inputfile into $outputfolder/$inputfilename.pos\n\n" if $verbose;
print LOG "Convert $inputfile into $outputfolder/$inputfilename.pos\n\n";

print LOG "$bed2pos $inputfile -o $outputfolder/$inputfilename.pos\n\n";
system("$bed2pos $inputfile -o $outputfolder/$inputfilename.pos");


#######
#Motif enrichment analysis
#######

print STDERR "Perform motif enrichment test for $outputfolder/$inputfilename.pos\n\n" if $verbose;
print LOG "Perform motif enrichment test for $outputfolder/$inputfilename.pos\n\n";


print LOG "$findmotifsgenome $outputfolder/$inputfilename.pos $genomeversion $outputfolder/findmotifsgenome -size given -preparsedDir $outputfolder/findmotifsgenome/preparsed\n\n";
system("$findmotifsgenome $outputfolder/$inputfilename.pos $genomeversion $outputfolder/findmotifsgenome -size given -preparsedDir $outputfolder/findmotifsgenome/preparsed");



#remove temporary folder
print LOG "rm -R $outputfolder/findmotifsgenome/preparsed\n\n";
system("rm -R $outputfolder/findmotifsgenome/preparsed");


#######
#Move to results folder
#######

system("mkdir $outputfolder/results");
#motif to region
system("cp $outputfolder/findmotifsbed/$inputfilename\_intersect_homer_anno.txt $outputfolder/results");
#motif enrichment
system("cp $outputfolder/findmotifsgenome/homerResults.html $outputfolder/results/$inputfilename\_homerResults.html");
system("cp $outputfolder/findmotifsgenome/knownResults.html $outputfolder/results/$inputfilename\_knownResults.html");


########
#Functions
########

sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}

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

