#!/usr/bin/perl -w
use strict;
use Getopt::Long;

#CutAdapt+FASTQC+RSEM+STAR


########
#Prerequisites
########

#none


########
#Interface
########


my $version="0.9";

my $usage="

mergefiles
version: $version
Usage: sbptools mergefiles -m model.txt -i file1.txt,file2.txt,file3.txt -o result.txt

Description: use a model file to merge all the files together. default case insensitive, and files have title

Parameters:

    --model|-m        model file
    --in|-i           input file(s) separated by \",\"
    --out|-o          output file

    --case|-c         model case sensitive or not [F]
    --title|-t        files have title or not[T]
    --uniq|-u         only use unique items from model [F]
	
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


my $modelfile;
my $infiles;
my $outfile;
my $case="F";
my $withtitle="T";
my $verbose;
my $unique="F";
my $runmode="none";

GetOptions(
	"model|m=s" => \$modelfile,
	"in|i=s" => \$infiles,
	"out|o=s" => \$outfile,

	"case|c=s" => \$case,
	"title|t=s" => \$withtitle,
	"uniq|u=s" => \$unique,

	"verbose|v" => \$verbose,
);


my $logfile=$outfile;
$logfile=~s/\.txt/_mergefiles.log/;



########
#Process
########

print STDERR "\nsbptools mergefiles $version\n\n";

#open(LOG,">$logfile") || die "Error writing $logfile. $!";
#read model file


my %models;
my @models;

open(IN,$modelfile) || die $!;

my $mline=0;
my $mmaxcol=0; #max col
my $mtitle;

while(<IN>) {
	tr/\r\n//d;
	my @array=split(/\t/,$_,-1);
	
	my $key;
	
	if($case eq "F") {
		$key=uc $array[0];
	}
	else {
		$key=$array[0];
	}	
	
	if($mline==0 && $withtitle eq "T") {
		$mtitle=$_;
	}
	else {
		$models{$key}++;
	}	
	
	if(@array>$mmaxcol) {
		$mmaxcol=@array;#max col defined by the longest row
	}

	
	$mline++;
}
close IN;

print STDERR "Model file $modelfile:$mline row(s) and $mmaxcol column(s) defined.\n" if $verbose;
#print LOG "Model file $modelfile:$mline row(s) and $mmaxcol column(s) defined.\n";


#reading infiles
my %imaxcols;
my %ititles;
my %inputs;

foreach my $infile (split(",",$infiles)) {

	my $iline=0;
	#my $ititle;
	my $imaxcol=0;
	my %infile;
	
	open(IN,$infile) || die $!;

	while(<IN>) {
		tr/\r\n//d;
		my @array=split(/\t/,$_,-1);
		
		my $key;
		if($case eq "F") {
			$key=uc $array[0];
		}
		else {
			$key=$array[0];
		}
		
		
		if($iline==0 && $withtitle eq "T") {
			#$ititle=$_;
			$ititles{$infile}=$_;
		}
		else {
			if(defined $models{$key}) {
				$inputs{$key}{$infile}=$_;
			}
		}
		
		if(@array>$imaxcol) {
			$imaxcol=@array;
			$imaxcols{$infile}=$imaxcol; #max col defined by the longest row
		}
		
		$iline++;
	}
	close IN;
	
	print STDERR "Input file $infile:$iline row(s) and $imaxcol column(s) defined.\n" if $verbose;
	#print LOG "Input file $infile:$iline row(s) and $imaxcol column(s) defined.\n";
}



#writing output files
open(OUT,">$outfile") || die $!;

#print title
my $title;
if($withtitle eq "T") {
	if(split("\t",$mtitle,-1)<$mmaxcol) {
		$title=join("\t",$mtitle,(" ") x ($mmaxcol - scalar(split("\t",$mtitle,-1))));
	}
	else {
		$title=$mtitle;
	}
	
	foreach my $infile (split(",",$infiles)) {
		my $ititle=$ititles{$infile};
		if(split("\t",$ititle,-1)<$imaxcols{$infile}) {
			$title.="\t".join("\t",$ititle,(" ") x ($imaxcols{$infile} - scalar(split("\t",$ititle,-1))));
		}
		else {
			$title.="\t".$ititle;
		}
	}
	print OUT "$title\n";
}

#print content
open(IN,$modelfile) || die $!;
$mline=0;
my %models2;
while(<IN>) {
	tr/\r\n//d;
	my @array=split(/\t/,$_,-1);
	
	my $key;
	my $content;
	
	if($case eq "F") {
		$key=uc $array[0];
	}
	else {
		$key=$array[0];
	}	
	
	unless($mline==0 && $withtitle eq "T") {
		unless($unique eq "T" && defined $models2{$key}) {
			$models2{$key}++;
			
			if(@array<$mmaxcol) {
				$content=join("\t",@array,(" ") x ($mmaxcol -@array));
			}
			else {
				$content=$_;
			}
					
			foreach my $infile (split(",",$infiles)) {
				if(defined $inputs{$key}{$infile}) {
					my $inputline=$inputs{$key}{$infile};
					if(split("\t",$inputline,-1)<$imaxcols{$infile}) {
						print STDERR $inputline,"#",scalar(split("\t",$inputline,-1)),"\n";
						$content.="\t".join("\t",$inputline,(" ") x ($imaxcols{$infile} -scalar(split("\t",$inputline,-1))));
					}
					else {
						$content.="\t".$inputline;
					}
				}
				else {
					$content.="\t".join("\t",(" ") x $imaxcols{$infile});
				}
			}
			
			print OUT $content,"\n";
		}
	}
	$mline++;
}
close IN;
close OUT;
#close LOG;

print STDERR "\nDone.\n";
