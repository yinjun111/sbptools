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


my $version="0.2";

#v0.2, add treatment and reference option


my $usage="

summarize_dm_peaks_bycalling
version: $version
Usage: perl  [parameters]

Description: summarize_dm_peaks


Mandatory Parameters:
    --in|-i           Input folder from chipseq-merge	
    --anno|-a         Peak annotation
    --out1|--o1       Output file1, annotated peak calling summary
    --out2|--o2       Output file2, summary of up,down,ns
	
    --tx              Transcriptome version, Human.B38.Ensembl84 or Mouse.B38.Ensembl84

    --treatment|-t    Treatment group name
    --reference|-r    Reference group name

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

my ($infile,$annofile,$outfile1,$outfile2);


my $verbose=1;
my $tx;
my $runmode="none";
my $treatment;
my $reference;

GetOptions(
	"in|i=s" => \$infile,
	"anno|a=s" => \$annofile,
	
	"out1|o1=s" => \$outfile1,
	"out2|o2=s" => \$outfile2,
	
	"tx|t=s" => \$tx,	

	"treatment=s" => \$treatment,
	"reference=s" => \$reference,
	
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

#read peak calling summary

my %peak2info;
my %peak2sig;

my $ilinenum=0;
my $ititle;
my %sample2pos;
my @selsamplepos;

open(IN,$infile) || die $!;
#infile: "resultbycalling"=> "all.reprod.peak.merged.summary.txt"

while(<IN>) {

	tr/\r\n//d;
	my @array=split/\t/;
	
	if($ilinenum==0) {

		
		for(my $num=1;$num<=(@array-1)/2;$num++) {
			if($array[$num]=~/(.+)_reprod.bed/) {
				$sample2pos{$1}=$num;
			}
		}
		
		print STDERR join(",",sort keys %sample2pos)," samples are defined in $infile.\n\n";
		
		#treatment, ref
		
		if(defined $sample2pos{$treatment}) {
			push @selsamplepos,$sample2pos{$treatment};
		}
		else {
			print STDERR "ERROR:--treatment $treatment not defined in $infile.\n\n";
			exit;
		}
		
		if(defined $sample2pos{$reference}) {
			push @selsamplepos,$sample2pos{$reference};
		}
		else {
			print STDERR "ERROR:--reference $reference not defined in $infile.\n\n";
			exit;
		}
		
		print STDERR "$treatment, $reference are identfied in columns ",join(",",map {$_+1} @selsamplepos),"\n";
		
		
		$ititle=join("\t",@array[0,@selsamplepos,map {$_+ (@array-1)/2} @selsamplepos]);		
		
		$ititle.="\tSignificance by Peak Calling:$treatment vs $reference";
	}
	else {
		$peak2info{$array[0]}=join("\t",@array[0,@selsamplepos,map {$_+ (@array-1)/2} @selsamplepos]);
				
		my $sig;	
		#called in first, 1, in second, -1, both, 0
		
		#treatment vs reference 
		
		if($array[$selsamplepos[0]]>0 && $array[$selsamplepos[1]] ==0) {
			$sig=1;
		}
		elsif($array[$selsamplepos[0]]==0 && $array[$selsamplepos[1]] >0) {
			$sig=-1;
		}
		else {
			$sig=0;
		}		
		
		$peak2sig{$array[0]}=$sig;
		
	}
	$ilinenum++;
}
close IN;
	

#read annotation file
my %gene2cate;
my %cates;

my $linenum=0;
open(IN,$annofile) || die $!;
open(OUT1,">$outfile1") || die $!;
#"resultbycallinganno"=> "all.reprod.peak.merged.dm.bycalling.anno.txt",

while(<IN>) {
	tr/\r\n//d;
	if ($linenum==0) {
		print OUT1 $ititle,"\t",$_,"\tCate\tTxAssigned\tGeneAssigned\t$txtitle\n";
	}
	else {
		
		my @array=split/\t/;
		
		if(defined $peak2info{$array[0]}) {
			
			#deal with annotation
			#0-18 are annotations 
			
			my $cate;
			my $txassigned;
			
			#TTS (ENST00000484188)
			
			if($array[7]=~/(.+) \((\w+)/) {
				$cate=$1;
				$txassigned=$2;
			}
			elsif($array[7] eq "Intergenic") {
				$cate="Intergenic";
				$txassigned=$array[10];
			}
			else{
				print STDERR $array[7],"\n";
				exit;
			}
			
			$cates{$cate}++;
			
			my $sig=$peak2sig{$array[0]};
			$gene2cate{$tx2gene{$txassigned}}{$cate}{"$sig"}{$array[0]}++;
			
			
			#print STDERR $txassigned,"#\t#$cate\n";
		
			print OUT1 join("\t",$peak2info{$array[0]},$peak2sig{$array[0]},@array,$cate,$txassigned,$tx2gene{$txassigned},$tx2info{$txassigned}),"\n";
		}
		else {
			print STDERR "ERROR:$array[0] not defined in $infile.\n$_\n";
			#print LOG "ERROR:$array[0] not defined in $infile.\n$_\n";
			exit;
		}

	}
	
	$linenum++;
}
close IN;
close OUT1;

#print STDERR Dumper(%gene2cate);


my @cates=("Intergenic","TTS","exon","intron","promoter-TSS");


open(OUT,">$outfile2") || die $!;
#"summarybycalling"=> "all.reprod.peak.merged.dm.bycalling.summary.txt",

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



