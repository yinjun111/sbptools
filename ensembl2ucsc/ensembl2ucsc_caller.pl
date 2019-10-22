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


my $version="0.2";

#v0.2, not fully assembled chrs are removed, due to difficulty of converting to UCSC


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
	
	my %chr2size;
	
	open(IN,$infile) || die "ERROR:Can't open $infile.$!.\n";
	open(OUT,">$outfile") || die "ERROR:Can't write into $outfile.$!.\n";
	my $chrprint=0;
	
	while(<IN>) {
		#need to double check the pattern for different species
		
		#only print \d+, X, Y, M
		if($_=~/^>(\w+)/) {
			my $chrname=$1;
			if( $chrname eq "MT") {
				print OUT ">chrM\n";	#change Mitochodria chr name
				$chrprint=1;
				print STDERR "chrM printed\n";
				print LOG "chrM printed.\n";
			}
			elsif( $chrname=~/^(\d+|X|Y)$/) {
				print OUT ">chr$chrname\n";
				$chrprint=1;
				print STDERR "chr$chrname printed.\n";
				print LOG "chr$chrname printed.\n";				
			}
			else {
				$chrprint=0;
				print STDERR "$chrname passed.\n";
				print LOG "$chrname passed.\n";
			}
		}
		else {
			if($chrprint==1) {
				print OUT $_;
			}	
		}
	}
	close IN;
	close OUT;
	
}
elsif($type =~/^gtf/i || $type =~/^bed/i) {
	#gtf file
	open(IN,$infile) || die "ERROR:Can't open $infile.$!.\n";
	open(OUT,">$outfile") || die "ERROR:Can't write into $infile.$!.\n";
	
	my $gtfprint=0;
	my %chrstat;
	
	while(<IN>) {
		#need to double check the pattern for different species
		my @array=split/\t/;
		
		if($_=~/^#/) {
			print OUT $_;
		}
		else {
			if($array[0]=~/^(\d+|X|Y)$/) {
				print OUT "chr$_";
				$chrstat{"printed"}{"chr$array[0]"}=1;
			}
			elsif($array[0]=~/^MT/) {
				print OUT join("\t","chrM",@array[1..$#array]);
				$chrstat{"printed"}{"chrM"}=1;
			}
			else {
				$chrstat{"passed"}{$array[0]}=1;
			}
		}
	}
	close IN;
	close OUT;
	
	print LOG join("\n", map {$_."\tprinted\n"} sort keys %{$chrstat{"printed"}}),"\n";
	print LOG join("\n", map {$_."\tpassed\n"} sort keys %{$chrstat{"passed"}}),"\n";

	print STDERR join("\n", map {$_."\tprinted\n"} sort keys %{$chrstat{"printed"}}),"\n";
	print STDERR join("\n", map {$_."\tpassed\n"} sort keys %{$chrstat{"passed"}}),"\n";	
	
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





