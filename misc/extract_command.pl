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

extract_command
version: $version
Usage: perl extract_command.pl -i command1.sh,command2.sh -c \'command\' -o newcommand.sh

Description: extract a command to run from a sbptools generated script file

Parameters:

    --in|-i           input shell script(s)
    --command|-c      command to extract
    --out|-o          output shell script

    --type|-t         command only or all after this [only]
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

my $infiles;
my $command;
my $outfile;
my $type="only";
my $verbose;

GetOptions(
	"in|i=s" => \$infiles,
	"out|o=s" => \$outfile,
	"command|c=s" => \$command,
	"type|t=s" => \$type,
	"verbose|v" => \$verbose,
);


########
#Process
########


open(OUT,">$outfile") || die "Error writing $outfile. $!";

foreach my $infile (split(",",$infiles)) {

	open(IN,$infile) || die "Error openning $infile. $!";
	while(<IN>) {
		tr/\r\n//d;
		if($type eq "only") {
			#only the selected command
			if($_=~/($command[^;]+)/) {
				print OUT $1,"\n";
			}
		}
		elsif($type=~/^(\d+)/) {
			#several commands after the selection
			my $commands;
			my $line;
			if($_=~/($command.+)/) {
				$line=$1;
				my $step=0;
				
				while($line=~/([^;]+)/g && $step<$type) {
					$commands.=$1.";";
					$step++;
				}
				print OUT $commands,"\n";
			}	
		}
		else {
			#everything after this command
			if($_=~/($command.+)/) {
				print OUT $1,"\n";
			}
		}
	}
	close IN;
}
close OUT;



