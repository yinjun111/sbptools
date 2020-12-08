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

#rnaseq-var-summary
#v0.2, summary genotype

#var-summary
#v0.1


my $usage="

var-summary
version: $version
Usage: Summarize variant genotype (GT), allele depth (AD) and impact for dnaseq-process and rnaseq-var

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

my $logfile="$outputfolder/var-summary_run.log";

#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

my $now=current_time();
my $timestamp=build_timestamp($now,"long");

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";



######
#Process input file
######


print STDERR "\nsbptools var-summary $version running ...\n\n" if $verbose;
print LOG "\nsbptools var-summary $version running ...\n\n";


print STDERR "Summarizing snpEff summary files:\n";
print LOG "Summarizing snpEff summary files:\n";

#
my %cate2col=qw(variants_impact_HIGH 4 variants_impact_LOW 5 variants_impact_MODERATE 6 variants_impact_MODIFIER 7);
my %samples;

#infiles,outfile,category
process_snpeff_gene_summary("$inputfolder/*/snpanno/*filtered.snpEff_summary.genes.txt","$outputfolder/var-summary_filtered.snpEff_summary.genes_merged.txt",$cate);
process_snpeff_gene_summary("$inputfolder/*/snpanno/*filtered-cleaned.snpEff_summary.genes.txt","$outputfolder/var-summary_filtered-cleaned.snpEff_summary.genes_merged.txt",$cate);


####
#GT/AD summary
####

print STDERR "\nSummarizing genotyping files:\n";
print LOG "\nSummarizing genotyping files:\n";


#infiles, outfile-prefix
process_snpeff_gtnr("$inputfolder/*/snpanno/*.filtered.snpeff.sift.annotated.edited.txt","$outputfolder/var-summary_filtered");
process_snpeff_gtnr("$inputfolder/*/snpanno/*.filtered-cleaned.snpeff.sift.annotated.edited.txt","$outputfolder/var-summary_filtered-cleaned");



####
#TMB summary
####

#use specific cates...




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



sub process_snpeff_gene_summary {
	my ($infiles,$resultfile,$catesel)=@_;

	#get all summary.genes.txt
	my @snpeff_infiles=glob($infiles);

	my %sample2snp;
	
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
			$sample2snp{$samplename}{$array[0]}=$array[$cate2col{$catesel}];
			$genes{$array[0]}++;
		}
		close IN;
	}

	print STDERR "Writing $resultfile\n";
	print LOG "Writing $resultfile\n";

	open(OUT,">$resultfile") || die $!;
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
}
	

sub process_snpeff_gtnr {
	
	my ($infiles,$outfile_prefix)=@_;
	
	
	#genotype
	my $gt_nr_summary="$outfile_prefix.gt_nr_summary.txt";
	my $gt_anno_summary="$outfile_prefix.gt_anno_summary.txt";

	#allele depth
	my $ad_nr_summary="$outfile_prefix.ad_nr_summary.txt";
	my $ad_anno_summary="$outfile_prefix.ad_anno_summary.txt";
	
	
	#get all summary.genes.txt
	my @infiles=glob($infiles);

	#my %sample2snp;
	#my %samples;
	my %varids;
	my %varidannos;

	my %sample2gt_nr;
	my %sample2gt_anno;

	my %sample2ad_nr;
	my %sample2ad_anno;


	foreach my $file (@infiles) {
		
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
			
			$sample2gt_nr{$samplename}{$varid}=$gt;
			$sample2ad_nr{$samplename}{$varid}=$scores[1];
			
			$varids{$varid}++;
				
			if(defined $array[14]) {		
				my $varid_anno=join(",",@array[0..4,13,11,12]);
				
				$sample2ad_anno{$samplename}{$varid_anno}=$scores[1];
				$sample2gt_anno{$samplename}{$varid_anno}=$gt;
				$varidannos{$varid_anno}++;
			}
		}
		close IN;
	}

	#GT, nofilter, NR

	print STDERR "\nWriting $gt_nr_summary\n";
	print LOG "\nWriting $gt_nr_summary\n";

	open(OUT,">$gt_nr_summary") || die $!;
	print OUT "Variant\t",join("\t",sort keys %samples),"\n";
	foreach my $variant (sort keys %varids) {
		print OUT $variant,"\t";
		my @marks;
		
		foreach my $sample (sort keys %samples) {
			if(defined $sample2gt_nr{$sample}{$variant}) {
				push @marks,$sample2gt_nr{$sample}{$variant};
			}
			else {
				push @marks,0;
			}
		}
		print OUT join("\t",@marks),"\n";
	}
	close OUT;
		

	#GT, nofilter, anno
	print STDERR "Writing $gt_anno_summary\n";
	print LOG "Writing $gt_anno_summary\n";

	open(OUT,">$gt_anno_summary") || die $!;
	print OUT "Variant\t",join("\t",sort keys %samples),"\n";
	foreach my $variant (sort keys %varidannos) {
		print OUT $variant,"\t";
		my @marks;
		
		foreach my $sample (sort keys %samples) {
			if(defined $sample2gt_anno{$sample}{$variant}) {
				push @marks,$sample2gt_anno{$sample}{$variant};
			}
			else {
				push @marks,0;
			}
		}
		print OUT join("\t",@marks),"\n";
	}
	close OUT;


	#AD, nofilter, NR
	print STDERR "Writing $ad_nr_summary\n";
	print LOG "Writing $ad_nr_summary\n";

	open(OUT,">$ad_nr_summary") || die $!;
	print OUT "Variant\t",join("\t",sort keys %samples),"\n";
	foreach my $variant (sort keys %varids) {
		print OUT $variant,"\t";
		my @marks;
		
		foreach my $sample (sort keys %samples) {
			if(defined $sample2ad_nr{$sample}{$variant}) {
				push @marks,$sample2ad_nr{$sample}{$variant};
			}
			else {
				push @marks,0;
			}
		}
		print OUT join("\t",@marks),"\n";
	}
	close OUT;

	#AD, nofilter, anno
	print STDERR "Writing $ad_anno_summary\n";
	print LOG "Writing $ad_anno_summary\n";

	open(OUT,">$ad_anno_summary") || die $!;
	print OUT "Variant\t",join("\t",sort keys %samples),"\n";
	foreach my $variant (sort keys %varidannos) {
		print OUT $variant,"\t";
		my @marks;
		
		foreach my $sample (sort keys %samples) {
			if(defined $sample2ad_anno{$sample}{$variant}) {
				push @marks,$sample2ad_anno{$sample}{$variant};
			}
			else {
				push @marks,0;
			}
		}
		print OUT join("\t",@marks),"\n";
	}
	close OUT;

}




