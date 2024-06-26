#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename);

########
#Updates
########

#v0.2, add option to submit each exp separately
#v0.3, combind up & down


########
#Prerequisites
########

#none


########
#Interface
########


my $version="0.3";


my $usage="

metascape-gen
version: $version
Usage: 

Description: 

Parameters:
    --in|-i           Input folder, DE significant matrix from rnaseq-summary
    --out|-o          Output file


	
	
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

my $infile;
my $outfile;

my $runmode=0;
my $verbose=1;


GetOptions(
	"in|i=s" => \$infile,
	"out|o=s" => \$outfile,

	"verbose" => \$verbose,
	"help|h" => sub {print STDERR $usage;exit;}
	
);


my $logfile="metascape-gen_run.log";


######
#Process input file
######

#combine up&down together
my $outfile_both=$outfile;
$outfile_both=~s/(\.\w+)$/_both$1/;

my $outfile_updown=$outfile;
$outfile_updown=~s/(\.\w+)$/_updown$1/;

#my $outfile_down=$outfile;
#$outfile_down=~s/(\.\w+)$/_down$1/;


my %comp2gene;
my @comps;
my %allgenes;

open(IN,$infile) || die "ERROR:Can't read $infile.$!\n";
my $linenum=0;
while(<IN>) {
	tr/\r\n//d;
	
	my @array=split/\t/;
	
	if($linenum==0) {
		@comps=@array[1..$#array];
	}
	else {
		for(my $num=1;$num<@array;$num++) {
			if($array[$num] ne " " && $array[$num] ne "NA" && $array[$num] == 1) {
				$comp2gene{$comps[$num-1]}{"both"}{$array[0]}++;
				$comp2gene{$comps[$num-1]}{"up"}{$array[0]}++;
			}
			elsif($array[$num] ne " " && $array[$num] ne "NA" && $array[$num] == -1) {
				$comp2gene{$comps[$num-1]}{"both"}{$array[0]}++;
				$comp2gene{$comps[$num-1]}{"down"}{$array[0]}++;
			}
		}
		$allgenes{$array[0]}++;
	}
	
	$linenum++;
}
close IN;

#All DE
open(OUT,">$outfile_both") || die "ERROR:Can't write $outfile.$!\n";
print OUT "#Comparison\tDE Genes\n";
foreach my $comp (@comps) {
	if(defined $comp2gene{$comp}{"both"}) {
		print OUT $comp,"\t",join(",",sort keys %{$comp2gene{$comp}{"both"}}),"\n";
		
		#output for each exp
		my $outfile_exp_both=$outfile;
		$outfile_exp_both=~s/(\.\w+)$/_$comp\_both$1/;
		
		open(OUT_EXP,">$outfile_exp_both") || die $!;
		
		print OUT_EXP $comp,"\t",join(",",sort keys %{$comp2gene{$comp}{"both"}}),"\n";
		print OUT_EXP "_BACKGROUND\t",join(",",sort keys %allgenes),"\n";
		
		close OUT_EXP;
		
	}
}
print OUT "_BACKGROUND\t",join(",",sort keys %allgenes),"\n";
close OUT;
		
#Up & Down
open(OUT,">$outfile_updown") || die "ERROR:Can't write $outfile_updown.$!\n";
print OUT "#Comparison\tDE Genes\n";
foreach my $comp (@comps) {
	my $outfile_exp_updown=$outfile;
	$outfile_exp_updown=~s/(\.\w+)$/_$comp\_updown$1/;
	
	open(OUT_EXP,">$outfile_exp_updown") || die $!;

	if(defined $comp2gene{$comp}{"up"}) {
		print OUT $comp,"-Up\t",join(",",sort keys %{$comp2gene{$comp}{"up"}}),"\n";
		
		#output for each exp
		print OUT_EXP $comp,"-Up\t",join(",",sort keys %{$comp2gene{$comp}{"up"}}),"\n";
	}
	
	if(defined $comp2gene{$comp}{"down"}) {
		print OUT $comp,"-Down\t",join(",",sort keys %{$comp2gene{$comp}{"down"}}),"\n";

		#output for each exp
		print OUT_EXP $comp,"-Down\t",join(",",sort keys %{$comp2gene{$comp}{"down"}}),"\n";
	}
	
	print OUT_EXP "_BACKGROUND\t",join(",",sort keys %allgenes),"\n";
	close OUT_EXP;
}
print OUT "_BACKGROUND\t",join(",",sort keys %allgenes),"\n";
close OUT;
