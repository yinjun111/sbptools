#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Basename qw(basename);


########
#Prerequisites
########



########
#Interface
########


my $version="0.3";

#v0.3 changed script name, and summarize mergebed output

my $usage="

process_mergebed
version: $version
Usage: perl rename_bed.pl -i file.bed -o file_rename.bed -s samplename

Description: Can be used to rename bed file, or summarize merged bed

Parameters:	
    --in|-i           Input file
    --source|-r       Source bed files for mergebed (optional)
    --out|-o          Output file
    --bysample|-s     Rename by sample name
    --bycoord|-c      Rename by coordinate
	
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

my ($infile,$outfile,$sourcefiles);
my $bysample;
my $bycoord;


GetOptions(
	"in|i=s" => \$infile,
	"source|r=s" => \$sourcefiles,	
	"out|o=s" => \$outfile,
	"bysample|s=s" => \$bysample,
	"bycoord|c" => \$bycoord,	
);


unless((defined $bysample && length($bysample)>0) || $bycoord) {
	print STDERR "ERROR:either --bysample or --bycoord needs to be defined.\n\n";
	exit;
}

########
#Process
#######

#read source file, find which peak exist in the source
my %peak2file;
my @sourcefiles_array;

if(defined $sourcefiles && length($sourcefiles)>0) {
	foreach my $file (split(",",$sourcefiles)) {
		my $filename=basename($file);
		push @sourcefiles_array,$filename;

		open(IN,$file) || die "ERROR:Can't open $file.$!\n";
		while(<IN>) {
			tr/\r\n//d;
			my @array=split/\t/;
			
			$peak2file{$array[3]}{$filename}++;
		}
		close IN;
	}
}	


#open merged peak, rename peak and write summary

my $summaryfile=$outfile;
$summaryfile=~s/\.\w+$/.summary.txt/;

open(IN,$infile) || die $!;
open(OUT,">$outfile") || die $!;

if(defined $sourcefiles && length($sourcefiles)>0) {
	open(SOUT,">$summaryfile") || die $!;
	print SOUT "MergedPeak\t",join("\t",@sourcefiles_array,@sourcefiles_array),"\n";
}

while(<IN>) {
	tr/\r\n//d;
	next if $_=~/^#/;
	my @array=split/\t/;
	
	#at least 4 columns need to be defined
	if(@array<4) {
		print STDERR "ERORR:at least 4 columns need to be defined in bedfile.\n$_\n";
		exit;
	}
		
	my $newname;
	my @marks;
	my @peaks;
	my %file2peak;
	
	if(defined $bysample && length($bysample)>0) {
		$newname=$array[3]."|".$bysample;
	}
	elsif($bycoord) {
		$newname=$array[0].":".($array[1]+1)."-".$array[2];
	}
	
	if(@array==4) {
		print OUT join("\t",@array[0..2],$newname),"\n";
	}
	else {
		print OUT join("\t",@array[0..2],$newname,@array[4..$#array]),"\n";
	}
	
	if(defined $sourcefiles && length($sourcefiles)>0) {
		#summarize peaks
		foreach my $peak (split(",",$array[3])) {
			foreach my $file (keys %{$peak2file{$peak}}) {
				$file2peak{$file}{$peak}++;
			}
		}
		
		for(my $num=0;$num<@sourcefiles_array;$num++) {
			if(defined $file2peak{$sourcefiles_array[$num]}) {
				$marks[$num]=scalar(keys %{$file2peak{$sourcefiles_array[$num]}});
				$peaks[$num]=join(",",sort keys %{$file2peak{$sourcefiles_array[$num]}});
			}
			else {
				$marks[$num]=0;
				$peaks[$num]=" ";
			}
		}
		
		print SOUT join("\t",$newname,@marks,@peaks),"\n";
	}

}
close IN;
close OUT;

if(defined $sourcefiles && length($sourcefiles)>0) {
	close SOUT;
}
