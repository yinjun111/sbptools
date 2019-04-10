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


my $version="0.93";

#future to be implemented, 1. don't retain input files rownames. 2. print to standard output for piping
#v0.92, better notice for file IO errors
#v0.93 add norownames

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

    --norownames|-n   Dont retain input file rownames [F]

    --column|-l       only retrieve specific column(s) from inputs
	
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
my $column;
my $norownames="F";

GetOptions(
	"model|m=s" => \$modelfile,
	"in|i=s" => \$infiles,
	"out|o=s" => \$outfile,

	"case|c=s" => \$case,
	"title|t=s" => \$withtitle,
	"uniq|u=s" => \$unique,
	"norownames|nr=s"=> \$norownames,
	
	"column|l=s" => \$column,

	"verbose|v" => \$verbose,
);


my $logfile=$outfile;
$logfile=~s/\.txt/_mergefiles.log/;

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

print STDERR "\nsbptools mergefiles $version\n\n";
print LOG "\nsbptools mergefiles $version running ...\n\n";

#read model file

my %models;
my @models;

open(IN,$modelfile) || die "ERROR:Can't read $modelfile.$!\n";

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
print LOG "Model file $modelfile:$mline row(s) and $mmaxcol column(s) defined.\n";


#reading infiles
my %imaxcols;
my %ititles;
my %inputs;

foreach my $infile (split(",",$infiles)) {

	my $iline=0;
	#my $ititle;
	my $imaxcol=0;
	my %infile;
	
	open(IN,$infile) || die "ERROR:Can't read $infile.$!\n";

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
			if($norownames eq "F") {
				$ititles{$infile}=$_;
			}
			else {
				$ititles{$infile}=join("\t",@array[1..$#array]);
			}
		}
		else {
			if(defined $models{$key}) {
				if($norownames eq "F") {
					$inputs{$key}{$infile}=$_;
				}
				else {
					$inputs{$key}{$infile}=join("\t",@array[1..$#array]);
				}
			}
		}
		
		if(@array>$imaxcol) {
			if($norownames eq "F") {
				$imaxcol=@array;
			}
			else {
				$imaxcol=@array-1;
			}
			$imaxcols{$infile}=$imaxcol; #max col defined by the longest row
		}
		
		$iline++;
	}
	close IN;
	
	print STDERR "Input file $infile:$iline row(s) and $imaxcol column(s) defined.\n" if $verbose;
	print LOG "Input file $infile:$iline row(s) and $imaxcol column(s) defined.\n";
}



#writing output files
open(OUT,">$outfile") || die "ERROR:Can't write $outfile.$!\n";

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
		my @ititle_array;
		
		if(split("\t",$ititle,-1)<$imaxcols{$infile}) {
			#implement specific cols here
			@ititle_array=(split("\t",$ititle,-1),(" ") x ($imaxcols{$infile} - scalar(split("\t",$ititle,-1))));
		}
		else {
			@ititle_array=split("\t",$ititle,-1);
		}
		
		#only use defined columns		
		if(defined $column && length($column)>0) {
			$title.="\t".join("\t",@ititle_array[map {$_-1} split(",",$column)]);
		}
		else {
			$title.="\t".join("\t",@ititle_array);
		}
		
	}
	print OUT "$title\n";
}

#print content
open(IN,$modelfile) || die "ERROR:Can't read $modelfile.$!\n";
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
			
			my @content_array;
			
			foreach my $infile (split(",",$infiles)) {
				if(defined $inputs{$key}{$infile}) {
					my $inputline=$inputs{$key}{$infile};
					if(split("\t",$inputline,-1)<$imaxcols{$infile}) {
						@content_array=(split("\t",$inputline,-1),(" ") x ($imaxcols{$infile} -scalar(split("\t",$inputline,-1))));
					}
					else {
						@content_array=split("\t",$inputline,-1);
					}
				}
				else {
					@content_array=((" ") x $imaxcols{$infile});
				}
				
				if(defined $column && length($column)>0) {
					$content.="\t".join("\t",@content_array[map {$_-1} split(",",$column)]);
				}
				else {
					$content.="\t".join("\t",@content_array);
				}
			}
			
			print OUT $content,"\n";
		}
	}
	$mline++;
}
close IN;
close OUT;


print STDERR "\nDone.\n";

print LOG "\nDone.\n";

close LOG;

########
#Functions
########

sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}