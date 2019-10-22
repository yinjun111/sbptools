#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);


########
#Prerequisites
########


my $mergefiles="perl /home/jyin/Projects/Pipeline/sbptools/mergefiles/mergefiles_caller.pl";

########
#Interface
########


my $version="0.3";

#v0.3, adding extra annotation for DE file

my $usage="

summarize_dm_peaks
version: $version
Usage: perl  [parameters]

Description: summarize_dm_peaks


Mandatory Parameters:
    --in|-i           Input folder from chipseq-merge
    --out|--o         Output file, summary of up,down,ns
	
    --tx|-t           Transcriptome version, Human.B38.Ensembl84 or Mouse.B38.Ensembl84
	
    --fccutoff        Log2 FC cutoff [1]
    --qcutoff         Corrected P cutoff [0.05]

    --runmode|-r      Where to run the scripts, local, server or none [none]
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

my ($infile,$annofile,$outfile);

my $fccutoff;
my $qcutoff;
	

my $verbose=1;
my $tx;
my $runmode="none";

GetOptions(
	"in|i=s" => \$infile,
	"anno|a=s" => \$annofile,
	
	"out|o=s" => \$outfile,
	
	"fccutoff=s" => \$fccutoff,
	"qcutoff=s" => \$qcutoff,

	"tx|t=s" => \$tx,	
	
	"runmode|r=s" => \$runmode,		
	"verbose|v" => \$verbose,
);



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


########
#Process
########

#read tx annotation
my %tx2gene;
my %tx2info;
my $txtitle;
my $tlinenum=0;
open(IN,$tx2ref{$tx}{"txanno"}) || die $!;
while(<IN>) {
	tr/\r\n//d;
	if($tlinenum==0) {
		$txtitle=$_;
	}
	else {
		my @array=split/\t/;
		$tx2gene{$array[0]}=$array[1];
		$tx2info{$array[0]}=$_;
	}
	$tlinenum++;
}
close IN;


#read DM file
my %gene2cate;
my %cates;

my $infile_rev=$infile;
$infile_rev=~s/.txt/_rev.txt/;

my $linenum=0;
open(IN,$infile) || die $!;
open(OUT,">$infile_rev") || die $!;
while(<IN>) {
	tr/\r\n//d;
	if ($linenum==0) {
	
		print OUT $_,"\tCate\tTxAssigned\tGeneAssigned\t$txtitle\n";
	}
	else {
	
		my @array=split/\t/;
		
		my $sig;
		
		#redefined significance
		if(defined $fccutoff && length($fccutoff)>0 && defined $qcutoff && length($qcutoff)>0) {
			if(2**$array[1]>= $fccutoff && $array[3] ne "NA" && $array[3] < $qcutoff) {
				$sig=1;
			} elsif (2**$array[1]<= 1/$fccutoff && $array[3] ne "NA" && $array[3] < $qcutoff) {
				$sig=-1;
			}
			else {
				$sig=0;
			}
		}
		else {
			$sig=$array[5];
		}
		
		my $cate;
		my $txassigned;
		
		#TTS (ENST00000484188)
		
		if($array[13]=~/(.+) \((\w+)/) {
			$cate=$1;
			$txassigned=$2;
		}
		elsif($array[13] eq "Intergenic") {
			$cate="Intergenic";
			$txassigned=$array[16];
		}
		else{
			print STDERR $array[13],"\n";
			exit;
		}
		
		$cates{$cate}++;
		$gene2cate{$tx2gene{$txassigned}}{$cate}{"$sig"}{$array[0]}++;
		
		#print STDERR $txassigned,"#\t#$cate\n";
		if(@array==25) {
			print OUT join("\t",@array[0..4],$sig,@array[6..$#array],$cate,$txassigned,$tx2gene{$txassigned},$tx2info{$txassigned}),"\n";
		}
		else {
			print OUT join("\t",@array[0..4],$sig,@array[6..$#array],(" ") x (25-@array),$cate,$txassigned,$tx2gene{$txassigned},$tx2info{$txassigned}),"\n";
		}
	}
	
	$linenum++;
}
close IN;
close OUT;

#print STDERR Dumper(%gene2cate);


my @cates=("Intergenic","TTS","exon","intron","promoter-TSS");


open(OUT,">$outfile") || die $!;

#1> -1 > 0

#print OUT1 "Gene\t",join("\t",@cates),"\n";
print OUT "Gene\t",join("\t",@cates),"\tUndirectional Total\tDirectional Total\t",join("\t",@cates),"\t",join("\t",@cates),"\n";

foreach my $gene (sort keys %gene2cate) {
	print OUT $gene,"\t";
	
	my @marks1;
	my @marks2;
	my @peaks;
	
	my $total_undir=0;
	my $total_dir=0;
	
	foreach my $cate (@cates) {
		my @cmarks1=();
		my $cmarks2=0;
		my @cpeaks;
		if(defined $gene2cate{$gene}{$cate}{"1"}) {
			push @cmarks1,scalar(keys %{$gene2cate{$gene}{$cate}{"1"}});
			$cmarks2=1;
			$total_undir+=1;
			$total_dir+=1;
			push @cpeaks,"1|".join(",",sort keys %{$gene2cate{$gene}{$cate}{"1"}});
		}
		else {
			push @cmarks1,0;
		}


		if (defined $gene2cate{$gene}{$cate}{"-1"}) {
			push @cmarks1,scalar(keys %{$gene2cate{$gene}{$cate}{"-1"}});
			if($cmarks2==0) {
				$cmarks2=-1;
				$total_undir+=1;
				$total_dir+=-1;				
			}
			push @cpeaks,"-1|".join(",",sort keys %{$gene2cate{$gene}{$cate}{"-1"}});
		}
		else {
			push @cmarks1,0;
		}
		
		if (defined $gene2cate{$gene}{$cate}{"0"}) {
			push @cmarks1,scalar(keys %{$gene2cate{$gene}{$cate}{"0"}});
			push @cpeaks,"0|".join(",",sort keys %{$gene2cate{$gene}{$cate}{"0"}});
		}
		else {
			push @cmarks1,0;
		}
		
		push @marks1,join(",",@cmarks1);
		push @marks2,$cmarks2;
		if(@cpeaks>0) {
			push @peaks,join(";",@cpeaks);
		}
		else {
			push @peaks," ";
		}
	}
	
	print OUT join("\t",@marks2),"\t$total_undir\t$total_dir\t";
	print OUT join("\t",@marks1),"\t";
	print OUT join("\t",@peaks),"\n";
}
close OUT;


########
#Functions
########

sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}
