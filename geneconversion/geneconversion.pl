#!/usr/bin/perl -w

use strict;
use Getopt::Long;


########
#Prerequisites
########

#access to gene conversion anntoation file


########
#Interface
########


my $version="0.2";

#v0.2, supports --anno

#

my $usage="

geneconversion
version: $version
Usage: sbptools geneconversion -i infile.txt -t Human -o outfile.txt

Description: Convert symbols/ids in infile to symbols/ids in annotation file

Parameters:

    --in|-i           input file, first column is gene symbol to be converted
    --out|-o          output file

    --anno|-a         annotation file. first column is gene symbol/id to be used. second column is gene symbol/id to be matched with input file.
                      By default, it is HUGO gene symbol and gene alias match (/data/jyin/Pipeline_Dev/geneconversion/geneconversion_anno_human.txt).

    --genecolumn|-g   Column number of gene symbol to be converted in input file [1]
    --inputdelim      Input file gene symbol/id delimiter [;]
    --annodelim       Annotation file gene symbol/id delimiter [;]
    --outputdelim     Output file gene symbol/id delimiter [;]
	
    --splitgene|-p    If synonym to gene is multiple match, whether to split results into multiple rows [1]
    --casesensitive|-c    Casesensitive match or not, 1 as case sensitive [1]
	
    --verbose|-v      Verbose


";

#    --species|-s      Species, currently supports Human, later Mouse to be implemented [Human]


unless (@ARGV) {
	print STDERR $usage;
	exit;
}

my $params=join(" ",@ARGV);
#then call different scripts


########
#Parameters
########


my $infile;
my $outfile;
my $annofile;
my $species="Human";
my $inputdelim=";";
my $annodelim=";";
my $outputdelim=";";

my $genecol=1;
my $casesensitive=0;
my $splitgene=1;
my $verbose;


GetOptions(
	"in|i=s" => \$infile,
	"out|o=s" => \$outfile,
	"anno|a=s" => \$annofile,	
	"t=s" => \$species,
	
	"genecolumn|g=s" => \$genecol,
	"inputdelim=s" => \$inputdelim,
	"annodelim=s" => \$annodelim,
	"outputdelim=s" => \$outputdelim,	
	
	"splitgene|p=s" => \$splitgene,
	"casesensitive|c=s"=>\$casesensitive,

	"verbose|v" => \$verbose,
);


my $logfile=$outfile;
$logfile=~s/\.txt/_geneconversion.log/;

#write log file
open(LOG, ">$logfile") || die "Error write $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";

print LOG "\n";




########
#Process
########

print STDERR "\nsbptools genecoversion $version\n\n";
print LOG "\nsbptools genecoversion $version running ...\n\n";

#need to change file location to /database
my %species2anno=(
	"Human"=>"/data/jyin/Pipeline_Dev/geneconversion/geneconversion_anno_human.txt",
	"Mouse"=>"",
);


unless(defined $species2anno{$species}) {
	print STDERR "ERROR:$species is not supported yet. Currently only supports ",join(",",sort keys %species2anno),"\n";
	print LOG "ERROR:$species is not supported yet. Currently only supports ",join(",",sort keys %species2anno),"\n";
	exit;
}

unless(defined $annofile && length($annofile)>0) {
	$annofile=$species2anno{$species};
}

if(!-e $annofile) {
	print STDERR "ERROR:$annofile doesn't exist.\n\n";
	exit;
}



#read  annotation
my %syn2gene;
open(IN,$annofile) || die "ERROR:Can't read ",$annofile,".$!";
while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^Gene/;
	
	my @array=split/\t/;
	
	#sometimes there may be no match column content
	if(defined $array[1] && length($array[1])>0) {	
		foreach my $syn (split($annodelim,$array[1])) {
			if($casesensitive) {
				$syn2gene{$syn}{$array[0]}++;
				$syn2gene{$array[0]}{$array[0]}++;
			}
			else {
				#not casesensitve 
				$syn2gene{uc $syn}{$array[0]}++;
				$syn2gene{uc $array[0]}{$array[0]}++;
			}
		}
	}
}
close IN;

#read  file

open(IN,$infile) || die "ERROR:Can't read $infile.$!";
open(OUT,">$outfile") || die "ERROR:Can't write $outfile.$!";

my $linenum=0; #file should contain title
while(<IN>) {
	tr/\r\n//d;
	my @array=split/\t/;
	
	if($linenum==0) {
		print OUT "Gene\t$_\n";
	}
	else {
		unless (defined $array[$genecol-1] && length($array[$genecol-1])>0) {
			print STDERR "ERROR:Column $genecol not defined in $infile.\n\n";
			exit;
		}
		
		#multiple gene symbols/ids
		my %genes;
		
		foreach  my $syn (split($inputdelim,$array[$genecol-1])) {
		
			if(!$casesensitive) {
				$syn=uc $syn;
			}
		
			if(defined $syn2gene{$syn}) {
				foreach my $gene (keys %{$syn2gene{$syn}}) {
					$genes{$gene}++;
				}
			}
		}
		
		if(keys %genes>1) {
			if($splitgene==1) {
				foreach my $gene (sort keys %genes) {
					print OUT $gene,"\t",$_,"\n";
				}
			}
			else {
				print OUT join($outputdelim,sort keys %genes),"\t",$_,"\n";
			}
		}
		elsif(keys %genes==1) {
			#unique match
			print OUT join($outputdelim,sort keys %genes),"\t",$_,"\n";
		}
		else {
			#no match
			print OUT $array[$genecol-1],"\t",$_,"\n";
		}
	}
	
	$linenum++;
}
close IN;
close OUT;

	
	


########
#Functions
########

sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}
