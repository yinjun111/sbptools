#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename qw(basename dirname);


########
#Updates
########



########
#Interface
########


my $version="0.12";

#v0.11, add heatmap generation to p/q matrix
#v0.12, R 4.0

my $usage="

gsea-gen-summary
version: $version
Usage: 

Description: 

Parameters:
    --in|-i           Input folder, containing GSEA results. Can be multiple folders separated by \",\"
    --out|-o          Output folder
    --db|-d           Selected databases to summarize, separated by \",\" (optional)
                          e.g. h.all.v7.1,c5.bp.v7.1
                      If not defined, by default, gsea-gen-summary will summarize all the gsea-gen results
    --pcutoff|-p      pvalue cutoff [0.05]
    --qcutoff|-q      FDR qval cutoff. Use \"-q 2\" to ensure q of 1 is not selected [2]
	
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

my $db;
my $dev=0;

my $runmode=0;
my $verbose=1;


GetOptions(
	"in|i=s" => \$infolders,
	"out|o=s" => \$outputfolder,
	"db|d=s" => \$db,
	"pcutoff|p=s" => \$pcutoff,
	"qcutoff|q=s" => \$qcutoff,
	
	"dev" => \$dev,
	"verbose" => \$verbose,
	"help|h" => sub {print STDERR $usage;exit;}
	
);


######
#Prerequisites
######


#my $r=find_program("/apps/R-3.4.1/bin/R");
#my $rscript=find_program("/apps/R-3.4.1/bin/Rscript");

my $r=find_program("/apps/R-4.0.2/bin/R");
my $rscript=find_program("/apps/R-4.0.2/bin/Rscript");

my $sbptoolsfolder="/apps/sbptools/";

#adding --dev switch for better development process
if($dev) {
	$sbptoolsfolder="/home/jyin/Projects/Pipeline/sbptools/";
}
else {
	#the tools called will be within the same folder of the script
	$sbptoolsfolder=get_parent_folder(dirname(abs_path($0)));
}



#sbptools
my $gs_heatmap="$sbptoolsfolder/gsea-gen-summary/gs_heatmap.R";


######
#Process input file
######


#Create folders

if(!-e $outputfolder) {
	mkdir($outputfolder);
}

$outputfolder = abs_path($outputfolder);


my $logfile="$outputfolder/gsea-gen-summary_run.log";

#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();
my $timestamp=build_timestamp($now,"long");


print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";
print LOG "\n";


my $outfile_p="$outputfolder/gsea-gen-summary_p.txt";
my $outfile_q="$outputfolder/gsea-gen-summary_q.txt";
my $outfile_sig="$outputfolder/gsea-gen-summary_sig.txt";
my $outfile_sig_anno="$outputfolder/gsea-gen-summary_sig_anno.txt";






######
#Program begins
#####

print STDERR "gsea-gen-summary $version running...\n\n";
print LOG "gsea-gen-summary $version running...\n\n";

#process input

my %seldbs;

if(defined $db && length($db)>0) {
	foreach my $dbsel (split(",",$db)) {
		my $db2file="/data/jyin/Databases/GSEA/msigdb_v7.1_files_to_download_locally/msigdb_v7.1_GMTs/"."$dbsel.symbols.gmt";

		if(-e $db2file) {
			$seldbs{$dbsel}++;
		}
		else {
			print STDERR "ERROR:--db $db not supported. $db2file was not found.\n\n";
			print LOG "ERROR:--db $db not supported. $db2file was not found.\n\n";
			exit;
		}
	}
}

if(defined $db && length($db)>0) {
	print STDERR "\n",join(",",sort keys %seldbs)," are needed for gsea-gen-summary.\n\n";
	print LOG "\n",join(",",sort keys %seldbs)," are needed for gsea-gen-summary.\n\n";
}
else {
	print STDERR "\nNo --db defined. All dbs are needed for gsea-gen-summary.\n\n";
	print LOG "\nNo --db defined. All dbs are needed for gsea-gen-summary.\n\n";
}


#process infiles

my %sample2gs_p;
my %sample2gs_q;
my %sample2gs_sig;

my %sample2gsbydb_p;
my %sample2gsbydb_q;
my %sample2gsbydb_sig;


my %samples;
my %gss;
my %gssbydb;
my $title;

my @infiles;
my %alldbs;

foreach my $infolder (split(",",$infolders)) {
	#get all knownResults.txt
	foreach my $infile (glob("$infolder/*/*/gsea_report_for_*.xls")) {
		my $infile_path=abs_path($infile);
		
		#print STDERR "$infile_path\n";
		print LOG "$infile_path\n";
	
		my $filename;
		my $dbused;
		
		#Summary2/forGSEA/CT_MG132_vs_CT_DMSO-c2.cp.kegg.v7.1/my_analysis.Gsea.1595291435163/
		if($infile=~/([^\/]+)-([ch][^-]+)\/my_analysis.Gsea.\d+\/gsea_report_for_([^\/]+)_\d+.xls/) {
			$filename=$1."-".$3;
			$dbused=$2;
			
			if(defined $db && length($db)>0) {
				unless(defined $seldbs{$dbused}) {
					print STDERR "In $infile, $dbused is not required by --db $db. Skipping...\n";
					next;
				}
			}
		}
		else {
			print STDERR "ERROR:$infile is not recognized as GSEA result.\n\n";
			next;
		}
		
		$alldbs{$dbused}++;
		
		push @infiles,$infile_path;

		$samples{$filename}++;
		
		open(IN,$infile) || die "ERROR:Can't read $infile.$!\n";
		while(<IN>) {
			tr/\r\n//d;
			if ($_=~/^NAME/) {
				$title=$_;
			}
			else {
			
				my @array=split/\t/;
				
				#p/q/sigs
				$sample2gs_p{$filename}{$array[0]}=$array[6];
				$sample2gs_q{$filename}{$array[0]}=$array[7];
				
				$sample2gsbydb_p{$dbused}{$filename}{$array[0]}=$array[6];
				$sample2gsbydb_q{$dbused}{$filename}{$array[0]}=$array[7];
				
				
				if($array[6]<$pcutoff && $array[7]<$qcutoff) {
					$sample2gs_sig{$filename}{$array[0]}=1;					
					$sample2gsbydb_sig{$dbused}{$filename}{$array[0]}=1;					
				}
				else {
					$sample2gs_sig{$filename}{$array[0]}=0;
					$sample2gsbydb_sig{$dbused}{$filename}{$array[0]}=1;					
				}				

				$gss{$array[0]}++;
				$gssbydb{$dbused}{$array[0]}++;
				
			}
		}
		close IN;
	}
}

print STDERR scalar(@infiles)," files detected for ",scalar(keys %samples)," comparisons and ",scalar(keys %alldbs)," databases.\n\n";


#write summary for all the dbs together

open(POUT,">$outfile_p") || die $!;
open(QOUT,">$outfile_q") || die $!;
open(SOUT,">$outfile_sig") || die $!;

print POUT "Name\t",join("\t",sort keys %samples),"\n";
print QOUT "Name\t",join("\t",sort keys %samples),"\n";
print SOUT "Name\t",join("\t",sort keys %samples),"\n";

foreach my $gs (sort keys %gss) {
	my @ps;
	my @qs;
	my @sigs;
	
	foreach my $sample (sort keys %samples) {
		if(defined $sample2gs_p{$sample}) {
			if(defined $sample2gs_p{$sample}{$gs}) {
				push @ps,$sample2gs_p{$sample}{$gs};
				push @qs,$sample2gs_q{$sample}{$gs};
				push @sigs,$sample2gs_sig{$sample}{$gs};
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
	
	print POUT $gs,"\t",join("\t",@ps),"\n";
	print QOUT $gs,"\t",join("\t",@qs),"\n";
	print SOUT $gs,"\t",join("\t",@sigs),"\n";
	
}
close POUT;
close QOUT;
close SOUT;




#write summary for each db
my @outfiles;

foreach my $dbsel (sort keys %alldbs) {
	if(defined $sample2gsbydb_sig{$dbsel}) {
		my $outfilebydb_p=$outfile_p;		
		$outfilebydb_p=~s/_p.txt/_$dbsel\_p.txt/;
		push @outfiles,$outfilebydb_p;

		my $outfilebydb_q=$outfile_q;
		$outfilebydb_q=~s/_q.txt/_$dbsel\_q.txt/;
		push @outfiles,$outfilebydb_q;
		
		my $outfilebydb_sig=$outfile_sig;
		$outfilebydb_sig=~s/_sig.txt/_$dbsel\_sig.txt/;
		
	
		open(POUT,">$outfilebydb_p") || die $!;
		open(QOUT,">$outfilebydb_q") || die $!;
		open(SOUT,">$outfilebydb_sig") || die $!;

		print POUT "Name\t",join("\t",sort keys %samples),"\n";
		print QOUT "Name\t",join("\t",sort keys %samples),"\n";
		print SOUT "Name\t",join("\t",sort keys %samples),"\n";

		foreach my $gs (sort keys %{$gssbydb{$dbsel}}) {
			my @ps;
			my @qs;
			my @sigs;
			
			foreach my $sample (sort keys %samples) {
				if(defined $sample2gsbydb_p{$dbsel}{$sample}) {
					if(defined $sample2gsbydb_p{$dbsel}{$sample}{$gs}) {
						push @ps,$sample2gsbydb_p{$dbsel}{$sample}{$gs};
						push @qs,$sample2gsbydb_q{$dbsel}{$sample}{$gs};
						push @sigs,$sample2gsbydb_sig{$dbsel}{$sample}{$gs};
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
			
			print POUT $gs,"\t",join("\t",@ps),"\n";
			print QOUT $gs,"\t",join("\t",@qs),"\n";
			print SOUT $gs,"\t",join("\t",@sigs),"\n";
			
		}
		close POUT;
		close QOUT;
		close SOUT;
	}
}


#Generate heatmap figures for p/q

foreach my $outfile (@outfiles) {
	my $outfilefig=$outfile;
	$outfilefig=~s/.txt/_heatmap.png/;	

	print STDERR "Generating heatmap for $outfile.\n";
	print LOG "Generating heatmap for $outfile.\n";

	system("$rscript $gs_heatmap -i $outfile -o $outfilefig");
	print LOG "$rscript $gs_heatmap -i $outfile -o $outfilefig\n";
	
}



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



sub find_program {
	my $fullprogram=shift @_;
	
	my $program;
	if($fullprogram=~/([^\/]+)$/) {
		$program=$1;
	}
	
	if(-e $fullprogram) {
		return $fullprogram;
	}
	else {
		my $sysout=`$program`;
		if($sysout) {
			my $location=`which $program`;
			return $location;
		}
		else {
			print STDERR "ERROR:$fullprogram or $program not found in your system.\n\n";
			exit;
		}
	}
}



sub get_parent_folder {
	my $dir=shift @_;
	
	if($dir=~/^(.+\/)[^\/]+\/?/) {
		return $1;
	}
}


