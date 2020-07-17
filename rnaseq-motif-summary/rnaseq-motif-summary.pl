#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);

########
#Updates
########



########
#Prerequisites
########

#none


########
#Interface
########


my $version="0.1";


my $usage="

rnaseq-motif-summary
version: $version
Usage: 

Description: 

Parameters:
    --in|-i           Input folder, containing rnaseq-motif results. Can be multiple folders separated by \",\"
    --out|-o          Output folder
	
    --pcutoff|-p      pvalue cutoff [0.05]
    --qcutoff|-q      Adjusted pvalue(q) cutoff. Use \"-q 2\" to ensure q of 1 is not selected [2]
	
";

#In Falco, BSR Linux server, use screen+parallel to parallelly running jobs in background.


unless (@ARGV) {
	print STDERR $usage;
	exit;
}

my $params=join(" ",@ARGV);
#then call different scripts


########
#Parameters
########

my $infolders;
my $outputfolder;

my $pcutoff=0.05;
my $qcutoff=2;

my $runmode=0;
my $verbose=1;


GetOptions(
	"in|i=s" => \$infolders,
	"out|o=s" => \$outputfolder,
	"pcutoff|p=s" => \$pcutoff,
	"qcutoff|q=s" => \$qcutoff,

	"verbose" => \$verbose,
	"help|h" => sub {print STDERR $usage;exit;}
	
);

my $homermotifannotation="/data/jyin/Databases/Homer/all.motifs.anno.togene.txt";



######
#Process input file
######


#Create folders

if(!-e $outputfolder) {
	mkdir($outputfolder);
}

$outputfolder = abs_path($outputfolder);


my $logfile="$outputfolder/rnaseq-motif-summary_run.log";

#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();
my $timestamp=build_timestamp($now,"long");


print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";
print LOG "\n";


my $outfile_p="$outputfolder/rnaseq-motif-summary_p.txt";
my $outfile_q="$outputfolder/rnaseq-motif-summary_q.txt";
my $outfile_sig="$outputfolder/rnaseq-motif-summary_sig.txt";
my $outfile_sig_anno="$outputfolder/rnaseq-motif-summary_sig_anno.txt";



# Read input

my @infiles;

foreach my $infolder (split(",",$infolders)) {
	#get all knownResults.txt
	foreach my $infile (glob("$infolder/*/*/findmotifsgenome/knownResults.txt")) {
		my $infile_path=abs_path($infile);
		print LOG "$infile_path\n";
	
		push @infiles,$infile_path;
	}
}


#process input


my %sample2motif_p;
my %sample2motif_q;
my %sample2motif_sig;

my %samples;
my %motifs;
my $title;

foreach my $file (@infiles) {
	
	my $filename;
	
	if($file=~/([^\/]+)\/findmotifsgenome/) {
		$filename=$1;
	}
	
	$samples{$filename}++;
	
	open(IN,$file) || die "ERROR:Can't read $file.$!\n";
	while(<IN>) {
		tr/\r\n//d;
		if ($_=~/^Motif Name/) {
			$title=$_;
		}
		else {
		
			my @array=split/\t/;
			
			#variants_impact_HIGH
			$sample2motif_p{$filename}{$array[0]}=$array[2];
			$sample2motif_q{$filename}{$array[0]}=$array[4];
			
			if($array[2]<$pcutoff && $array[4]<$qcutoff) {
				$sample2motif_sig{$filename}{$array[0]}=1;
			}
			else {
				$sample2motif_sig{$filename}{$array[0]}=0;
			}
			
			$motifs{$array[0]}++;
		}
	}
	close IN;
}


#read motif annotation

my %motif2anno;
my $annotitle;

open(IN,$homermotifannotation) || die $!;
while(<IN>) {
	tr/\r\n//d;
	if($_=~/^Gene\t/) {
		$annotitle=$_;
	}
	else {
		my @array=split/\t/;
		$motif2anno{$array[7]}{$_}++;
	}
}
close IN;


#write summary

open(POUT,">$outfile_p") || die $!;
open(QOUT,">$outfile_q") || die $!;
open(SOUT,">$outfile_sig") || die $!;
open(SAOUT,">$outfile_sig_anno") || die $!;

print POUT "Motif\t",join("\t",sort keys %samples),"\n";
print QOUT "Motif\t",join("\t",sort keys %samples),"\n";
print SOUT "Motif\t",join("\t",sort keys %samples),"\n";
print SAOUT "$annotitle\t",join("\t",sort keys %samples),"\n";


foreach my $motif (sort keys %motifs) {
	my @ps;
	my @qs;
	my @sigs;
	
	foreach my $sample (sort keys %samples) {
		if(defined $sample2motif_p{$sample}) {
			if(defined $sample2motif_p{$sample}{$motif}) {
				push @ps,$sample2motif_p{$sample}{$motif};
				push @qs,$sample2motif_q{$sample}{$motif};
				push @sigs,$sample2motif_sig{$sample}{$motif};
			}
			else {
				push @ps," ";
				push @qs," ";
				push @sigs," ";
			}
		}
		else {
				push @ps," ";
				push @qs," ";
				push @sigs," ";			
		}
	}
	
	print POUT $motif,"\t",join("\t",@ps),"\n";
	print QOUT $motif,"\t",join("\t",@qs),"\n";
	print SOUT $motif,"\t",join("\t",@sigs),"\n";
	
	if(defined $motif2anno{$motif}) {
		foreach my $anno (sort keys %{$motif2anno{$motif}}) {
			print SAOUT $anno,"\t",join("\t",@sigs),"\n";
		}
	}
	else {
		print STDERR "ERROR:No annotation for $motif.\n";
	}	
}
close POUT;
close QOUT;
close SOUT;
close SAOUT;




#########
# Functions
#########


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


