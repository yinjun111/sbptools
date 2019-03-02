#!/usr/bin/perl -w
use strict;
use Getopt::Long;


########
#Prerequisites
########

#none


########
#Interface
########


my $version="0.1";

my $usage="

ensembl2ucsc
version: $version
Usage: sbptools ensembl2ucsc -i genome.fa -t fastq -o genome_ucsc.fa

Description: Convert Ensembl genome fasta/gtf/bed files into UCSC format

Parameters:

    --in|-i           input file
    --out|-o          out file
    --type|-t         FASTA, GTF or BED

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

my $infile;
my $outfile;
my $type;
my $verbose;

GetOptions(
	"in|i=s" => \$infile,
	"out|o=s" => \$outfile,
	"type|t=s" => \$type,

	"verbose|v" => \$verbose,
);


my $logfile=$outfile;
$logfile=~s/\.\w+$/_ensembl2ucsc.log/;


#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";
print LOG "\n";



########
#Process
########

print STDERR "\nsbptools ensembl2ucsc $version\n\n";
print LOG "\nsbptools ensembl2ucsc $version\n\n";

print STDERR "Converting $infile into $outfile\n\n";
print LOG "Converting $infile into $outfile\n\n";

#

if($type =~/^fasta/i) {
	#fasta file
	open(IN,$infile) || die "ERROR:Can't open $infile.$!.\n";
	open(OUT,">$outfile") || die "ERROR:Can't write into $infile.$!.\n";
	while(<IN>) {
		#need to double check the pattern for different species
		if($_=~/^>(\w+)/) {
			print OUT ">chr$1\n";
		}
		else {
			print OUT $_;
		}
	}
	close IN;
	close OUT;
}
elsif($type =~/^gtf/i) {
	#gtf file
	open(IN,$infile) || die "ERROR:Can't open $infile.$!.\n";
	open(OUT,">$outfile") || die "ERROR:Can't write into $infile.$!.\n";
	while(<IN>) {
		#need to double check the pattern for different species
		if($_=~/^#/) {
			print OUT $_;
		}
		else {
			print OUT "chr$_";
		}
	}
	close IN;
	close OUT;
}
elsif (	$type =~/^bed/i) {
	#gtf file
	open(IN,$infile) || die "ERROR:Can't open $infile.$!.\n";
	open(OUT,">$outfile") || die "ERROR:Can't write into $infile.$!.\n";
	while(<IN>) {
		#need to double check the pattern for different species
		if($_=~/^#/) {
			print OUT $_;
		}
		else {
			print OUT "chr$_";
		}
	}
	close IN;
	close OUT;
}
else {
	print STDERR "ERROR:--type $type not supported. Only support FASTA,GTF,BED for now.\n\n";

}

print STDERR "\nDone.\n\n";
print LOG "\nDone.\n\n";

close LOG;

#####
#Functions
#####

sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}





