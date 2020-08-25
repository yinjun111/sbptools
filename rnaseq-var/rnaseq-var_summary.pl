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


my $version="0.2";

#v0.2, summary genotype



my $usage="

rnaseq-var_summary
version: $version
Usage: Summarize variant genotype (GT), allele depth (AD) and impact for rnaseq-var.

Description: 

Parameters:
    --in|-i           Input folder, containing rnaseq-var results
    --out|-o          Output folder
    --cate|-c         Category, variants_impact_HIGH (default), variants_impact_LOW, variants_impact_MODERATE or variants_impact_MODIFIER

	
	
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

my $inputfolder;
my $outputfolder;
my $cate="variants_impact_HIGH";

my $runmode=0;
my $verbose=1;


GetOptions(
	"in|i=s" => \$inputfolder,
	"out|o=s" => \$outputfolder,
	"cate|c=s" => \$cate,

	"verbose" => \$verbose,
	"help|h" => sub {print STDERR $usage;exit;}
	
);



if(!-e $outputfolder) {
	mkdir($outputfolder);
}

$outputfolder = abs_path($outputfolder);

my $logfile="$outputfolder/rnaseq-var_summary_run.log";

#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();
my $timestamp=build_timestamp($now,"long");

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";



#output files

my $snpeffsummary="$outputfolder/rnaseq-var_snpeff_summary_merged.txt";

#genotype
my $gt_nofilter_nr_summary="$outputfolder/rnaseq-var_gt_nofilter_nr_summary.txt";
my $gt_nofilter_anno_summary="$outputfolder/rnaseq-var_gt_nofilter_anno_summary.txt";

#allele depth
my $ad_nofilter_nr_summary="$outputfolder/rnaseq-var_ad_nofilter_nr_summary.txt";
my $ad_nofilter_anno_summary="$outputfolder/rnaseq-var_ad_nofilter_anno_summary.txt";

#genotype
my $gt_filtered_nr_summary="$outputfolder/rnaseq-var_gt_filtered_nr_summary.txt";
my $gt_filtered_anno_summary="$outputfolder/rnaseq-var_gt_filtered_anno_summary.txt";

#allele depth
my $ad_filtered_nr_summary="$outputfolder/rnaseq-var_ad_filtered_nr_summary.txt";
my $ad_filtered_anno_summary="$outputfolder/rnaseq-var_ad_filtered_anno_summary.txt";


######
#Process input file
######


print STDERR "\nsbptools rnaseq-var-summary $version running ...\n\n" if $verbose;
print LOG "\nsbptools rnaseq-var-summary $version running ...\n\n";


print STDERR "Summarizing snpEff summary files:\n";
print LOG "Summarizing snpEff summary files:\n";

#
my %cate2col=qw(variants_impact_HIGH 4 variants_impact_LOW 5 variants_impact_MODERATE 6 variants_impact_MODIFIER 7);


#get all summary.genes.txt
my @snpeff_infiles=glob("$inputfolder/*/snpanno/*snpEff_summary.genes.txt");

my %sample2snp;
my %samples;
my %genes;

foreach my $file (@snpeff_infiles) {
	
	print STDERR "$file\n";
	print LOG "$file\n";
	
	my $filename=basename($file);
	my $samplename;
	
	if($file=~/([^\/]+)\/snpanno/) {
		$samplename=$1;
		$samples{$samplename}++;
	}
	
	open(IN,$file) || die "ERROR:Can't read $file.$!\n";
	while(<IN>) {
		tr/\r\n//d;
		next if $_=~/^#/;
		
		my @array=split/\t/;
		
		#variants_impact_HIGH
		$sample2snp{$samplename}{$array[0]}=$array[$cate2col{$cate}];
		$genes{$array[0]}++;
	}
	close IN;
}

print STDERR "Writing $snpeffsummary\n";
print LOG "Writing $snpeffsummary\n";

open(OUT,">$snpeffsummary") || die $!;
print OUT "Gene\t",join("\t",sort keys %samples),"\n";

foreach my $gene (sort keys %genes) {
	my @marks;
	
	foreach my $sample (sort keys %samples) {
		if(defined $sample2snp{$sample}) {
			if(defined $sample2snp{$sample}{$gene}) {
				push @marks,$sample2snp{$sample}{$gene};
			}
			else {
				push @marks,0;
			}
		}
		else {
			push @marks,0;
		}
	}
	
	print OUT $gene,"\t",join("\t",@marks),"\n";
}
close OUT;

	

####
#GT summary, no filter
####

print STDERR "\nSummarizing genotyping files:\n";
print LOG "\nSummarizing genotyping files:\n";


#get all summary.genes.txt
my @nofilter_infiles=glob("$inputfolder/*/snpanno/*sortedByCoord.outoutput.filtered.snpeff.sift.annotated.edited.txt");

#my %sample2snp;
#my %samples;
my %varids_nofilter;
my %varidannos_nofilter;

my %sample2gt_nofilter_nr;
my %sample2gt_nofilter_anno;

my %sample2ad_nofilter_nr;
my %sample2ad_nofilter_anno;


foreach my $file (@nofilter_infiles) {
	
	print STDERR "$file\n";
	print LOG "$file\n";
	my $samplename;
	
	#should be the same as the previous one
	if($file=~/([^\/]+)\/snpanno/) {
		$samplename=$1;
	}
	
	open(IN,$file) || die "ERROR:Can't read $file.$!\n";
	while(<IN>) {
		tr/\r\n//d;
		next if $_=~/^#/;
		
		my @array=split/\t/;
		
		#chr,start,name,ref,alt
		my $varid=join(",",@array[0..4]);
		
		#gt, 0, 1 (0/1), 2 (1/1)		
		
		my @scores=split(":",$array[9]);
		my $gt=0;
		
		if($scores[0] eq "1/1") {
			$gt=2;
		}
		elsif($scores[0] eq "0/1") {
			$gt=1;
		}
		
		$sample2gt_nofilter_nr{$samplename}{$varid}=$gt;
		$sample2ad_nofilter_nr{$samplename}{$varid}=$scores[1];
		
		$varids_nofilter{$varid}++;
			
		if(defined $array[14]) {		
			my $varid_anno=join(",",@array[0..4,13,11,12]);
			
			$sample2ad_nofilter_anno{$samplename}{$varid_anno}=$scores[1];
			$sample2gt_nofilter_anno{$samplename}{$varid_anno}=$gt;
			$varidannos_nofilter{$varid_anno}++;
		}
	}
	close IN;
}

#GT, nofilter, NR

print STDERR "\nWriting $gt_nofilter_nr_summary\n";
print LOG "\nWriting $gt_nofilter_nr_summary\n";

open(OUT,">$gt_nofilter_nr_summary") || die $!;
print OUT "Variant\t",join("\t",sort keys %samples),"\n";
foreach my $variant (sort keys %varids_nofilter) {
	print OUT $variant,"\t";
	my @marks;
	
	foreach my $sample (sort keys %samples) {
		if(defined $sample2gt_nofilter_nr{$sample}{$variant}) {
			push @marks,$sample2gt_nofilter_nr{$sample}{$variant};
		}
		else {
			push @marks,0;
		}
	}
	print OUT join("\t",@marks),"\n";
}
close OUT;
	

#GT, nofilter, anno
print STDERR "Writing $gt_nofilter_anno_summary\n";
print LOG "Writing $gt_nofilter_anno_summary\n";

open(OUT,">$gt_nofilter_anno_summary") || die $!;
print OUT "Variant\t",join("\t",sort keys %samples),"\n";
foreach my $variant (sort keys %varidannos_nofilter) {
	print OUT $variant,"\t";
	my @marks;
	
	foreach my $sample (sort keys %samples) {
		if(defined $sample2gt_nofilter_anno{$sample}{$variant}) {
			push @marks,$sample2gt_nofilter_anno{$sample}{$variant};
		}
		else {
			push @marks,0;
		}
	}
	print OUT join("\t",@marks),"\n";
}
close OUT;


#AD, nofilter, NR
print STDERR "Writing $ad_nofilter_nr_summary\n";
print LOG "Writing $ad_nofilter_nr_summary\n";

open(OUT,">$ad_nofilter_nr_summary") || die $!;
print OUT "Variant\t",join("\t",sort keys %samples),"\n";
foreach my $variant (sort keys %varids_nofilter) {
	print OUT $variant,"\t";
	my @marks;
	
	foreach my $sample (sort keys %samples) {
		if(defined $sample2ad_nofilter_nr{$sample}{$variant}) {
			push @marks,$sample2ad_nofilter_nr{$sample}{$variant};
		}
		else {
			push @marks,0;
		}
	}
	print OUT join("\t",@marks),"\n";
}
close OUT;

#AD, nofilter, anno
print STDERR "Writing $ad_nofilter_anno_summary\n";
print LOG "Writing $ad_nofilter_anno_summary\n";

open(OUT,">$ad_nofilter_anno_summary") || die $!;
print OUT "Variant\t",join("\t",sort keys %samples),"\n";
foreach my $variant (sort keys %varidannos_nofilter) {
	print OUT $variant,"\t";
	my @marks;
	
	foreach my $sample (sort keys %samples) {
		if(defined $sample2ad_nofilter_anno{$sample}{$variant}) {
			push @marks,$sample2ad_nofilter_anno{$sample}{$variant};
		}
		else {
			push @marks,0;
		}
	}
	print OUT join("\t",@marks),"\n";
}
close OUT;





####
#GT summary, filtered
####

print STDERR "\nSummarizing genotyping files:\n";
print LOG "\nSummarizing genotyping files:\n";


#get all summary.genes.txt
my @filtered_infiles=glob("$inputfolder/*/snpanno/*sortedByCoord.outoutput.afterfiltered.snpeff.sift.annotated.edited.txt");

#my %sample2snp;
#my %samples;
my %varids_filtered;
my %varidannos_filtered;

my %sample2gt_filtered_nr;
my %sample2gt_filtered_anno;

my %sample2ad_filtered_nr;
my %sample2ad_filtered_anno;


foreach my $file (@filtered_infiles) {
	
	print STDERR "$file\n";
	print LOG "$file\n";
	my $samplename;
	
	#should be the same as the previous one
	if($file=~/([^\/]+)\/snpanno/) {
		$samplename=$1;
	}
	
	open(IN,$file) || die "ERROR:Can't read $file.$!\n";
	while(<IN>) {
		tr/\r\n//d;
		next if $_=~/^#/;
		
		my @array=split/\t/;
		
		#chr,start,name,ref,alt
		my $varid=join(",",@array[0..4]);
		
		#gt, 0, 1 (0/1), 2 (1/1)		
		
		my @scores=split(":",$array[9]);
		my $gt=0;
		
		if($scores[0] eq "1/1") {
			$gt=2;
		}
		elsif($scores[0] eq "0/1") {
			$gt=1;
		}
		
		$sample2gt_filtered_nr{$samplename}{$varid}=$gt;
		$sample2ad_filtered_nr{$samplename}{$varid}=$scores[1];
		
		$varids_filtered{$varid}++;
			
		if(defined $array[14]) {		
			my $varid_anno=join(",",@array[0..4,13,11,12]);
			
			$sample2ad_filtered_anno{$samplename}{$varid_anno}=$scores[1];
			$sample2gt_filtered_anno{$samplename}{$varid_anno}=$gt;
			$varidannos_filtered{$varid_anno}++;
		}
	}
	close IN;
}

#GT, filtered, NR

print STDERR "\nWriting $gt_filtered_nr_summary\n";
print LOG "\nWriting $gt_filtered_nr_summary\n";

open(OUT,">$gt_filtered_nr_summary") || die $!;
print OUT "Variant\t",join("\t",sort keys %samples),"\n";
foreach my $variant (sort keys %varids_filtered) {
	print OUT $variant,"\t";
	my @marks;
	
	foreach my $sample (sort keys %samples) {
		if(defined $sample2gt_filtered_nr{$sample}{$variant}) {
			push @marks,$sample2gt_filtered_nr{$sample}{$variant};
		}
		else {
			push @marks,0;
		}
	}
	print OUT join("\t",@marks),"\n";
}
close OUT;
	

#GT, filtered, anno
print STDERR "Writing $gt_filtered_anno_summary\n";
print LOG "Writing $gt_filtered_anno_summary\n";

open(OUT,">$gt_filtered_anno_summary") || die $!;
print OUT "Variant\t",join("\t",sort keys %samples),"\n";
foreach my $variant (sort keys %varidannos_filtered) {
	print OUT $variant,"\t";
	my @marks;
	
	foreach my $sample (sort keys %samples) {
		if(defined $sample2gt_filtered_anno{$sample}{$variant}) {
			push @marks,$sample2gt_filtered_anno{$sample}{$variant};
		}
		else {
			push @marks,0;
		}
	}
	print OUT join("\t",@marks),"\n";
}
close OUT;


#AD, filtered, NR
print STDERR "Writing $ad_filtered_nr_summary\n";
print LOG "Writing $ad_filtered_nr_summary\n";

open(OUT,">$ad_filtered_nr_summary") || die $!;
print OUT "Variant\t",join("\t",sort keys %samples),"\n";
foreach my $variant (sort keys %varids_filtered) {
	print OUT $variant,"\t";
	my @marks;
	
	foreach my $sample (sort keys %samples) {
		if(defined $sample2ad_filtered_nr{$sample}{$variant}) {
			push @marks,$sample2ad_filtered_nr{$sample}{$variant};
		}
		else {
			push @marks,0;
		}
	}
	print OUT join("\t",@marks),"\n";
}
close OUT;

#AD, filtered, anno
print STDERR "Writing $ad_filtered_anno_summary\n";
print LOG "Writing $ad_filtered_anno_summary\n";

open(OUT,">$ad_filtered_anno_summary") || die $!;
print OUT "Variant\t",join("\t",sort keys %samples),"\n";
foreach my $variant (sort keys %varidannos_filtered) {
	print OUT $variant,"\t";
	my @marks;
	
	foreach my $sample (sort keys %samples) {
		if(defined $sample2ad_filtered_anno{$sample}{$variant}) {
			push @marks,$sample2ad_filtered_anno{$sample}{$variant};
		}
		else {
			push @marks,0;
		}
	}
	print OUT join("\t",@marks),"\n";
}
close OUT;



close LOG;


####
#Functions
####



sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
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






