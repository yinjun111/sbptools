#!/usr/bin/perl
#Andrew Hodges, PhD
#BI Shared Resource, SBP
#11/6/2019
#Update: 4/30/2020
#(C)2019-2020, SBP

#Purpose: generate gct and cls format files for sbptools expression data
#Input format based on sbptools outputs
#

#!/usr/bin/perl -w
#use strict;
#use warnings;
use Getopt::Long;
use List::MoreUtils qw(uniq);
#use String::Util qw(trim);


########
#Prerequisites
########


########
#Interface
########

########
#Update: removing rows with unmatched gene anno; updating cls file counts
########


my $version="0.5.4";


my $usage="
gseagen
version: $version\n
Usage:  perl gsea-gen_caller.pl -e ./data/gene.results.merged.tpm.txt -ga ./data/geneanno.txt -g example.gct -s ./data/configV1.txt -n Group -c example.gct \n

Description: Perl script to generate gct and cls files for GSEA analysis.\n
Parameters:
	--expression|-e      input file of merged expression data
	--sampleAnno|-s      sample annotation 
	--geneAnno|-ga       gene annotation (optional)
	--gctName|-g         output gct file name
	--clsName|-c         output cls file name
	--groupName|-n       name of group column for CLS generation
        --removeUnmatched|-u remove unmatched genes: 1 for remove rows (default/hardcoded)

#--verbose|-v     Verbose\n
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
my $expression;
my $sampleAnno;
my $geneAnno = "";
my $gctName;
my $clsName; ##="F";
my $sampleName="Sample";
#my $dirExprF;
#my $dirSanno; # = "";
#my $dirGanno;
#my $verbose;
my $groupName;
my $unmatched = 1;  #by default we will remove unmatched anno if 


####re-add here
GetOptions(	
    "expression|e=s" => \$expression,
    "sampleAnno|s=s" => \$sampleAnno,
    "geneAnno|ga=s" => \$geneAnno,
    "gctName|g=s" => \$gctName,
    "clsName|c=s" => \$clsName,
    "groupName|n=s" => \$groupName,
    "sampleName|r=s" => \$sampleName, #sample column ID used in config file
    #"removeUnmatched|u=s" => \$unmatched,	  
#    "verbose|v" => \$verbose,
);


#####MAIN functions:
##create a log file based on gct name: updated 4/20/2020
my $tmp = $gctName;
$tmp =~ s/\.gct/_gsea-gen_run.log/;
open LOG, ">".$tmp or die "Error: Unable to open logfile $tmp $! \n";
print LOG "GSEA-Gen Log file: \n";
my $now = localtime();
print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print STDERR "Start time: $now\n\n";
print LOG "Curent version: $version\n\n";
print LOG "\n";


#next get annotation
my %anno;
my $aflag = 0; #aflag is whether or not to use the gene annotation file to label rows
print LOG "Annotation: ";
print STDERR "Annotation: ";
if($geneAnno ne ""){ %anno = getAnno( $geneAnno); $aflag=1; print LOG "$geneAnno";}
else{ print LOG "No gene annotation file selected.";}
print LOG "\n";
my $rowcount;
my $colcount;


#$rowcount = 47643;
my $expr = "wc -l $expression";
$rowcount = `$expr` - 1;
print LOG "Number of starting rows: ".$rowcount."\n";
print STDERR "Number of starting rows: ".$rowcount."\n";

#$colcount = 20;
my $expr2 = "head -n2 $expression | tail -n1 | wc -w";
$colcount = `$expr2` - 1; #assume 1 gene info columns
print LOG "Number of columns: ".$colcount."\n";
print STDERR "Number of columns: ".$colcount."\n";

#Next open main filehandles for data and output gct file
open fout, ">".$gctName or die "Error: Can't write $gctName $! \n";
print LOG "\nOpened GCT file for writing.\n";
print STDERR "\nOpened GCT file for writing.\n";

open expr, $expression or die "Error: Can't read $expression $! \n";
#printf fout "#1.2\n".$rowcount."\t".$colcount."\n";
print LOG "Opened expression file for reading.\n";
print STDERR "Opened expression file for reading.\n";
my $tempgct = "";

my @sampleOrders;


####this can be added into a subroutine later
my $rowcount2 = 0;  ###counter for rows that arent removed during optional gene anno match
my $rowcount1 = 0;
my $y = -1;
while(<expr>){
	$y++; 
	#print $y."\n";
	chomp(my $string = $_);
	#print $string."\n";
	if($string ne ""){
		if($y==0){ #modify header column with "ID" and "Description" tags
			$string =~ s/^Gene/ID\tDescription/;
			my @tempr = split(/\t/,$string);
			@sampleOrders = @tempr;
			shift(@sampleOrders);
			shift(@sampleOrders);
			$tempgct .= $string."\n";
		}
		else{ ###data non-blank rows
			###first get the gene using annotation hash
			my @vals = split(/\t/,$string);
			my $gene2 = $vals[0];
			$gene2 =~ s/^[\s]+//;
			$gene2 =~ s/[\s]+$//;
			my $gene = ""; 
			if($aflag){ #use gene anno
			    			    
			    ####Update 4/24: only write if entry exists
			    if(exists $anno{ $gene2 }){
				#if($anno{ $gene2 } ne ""){
				$gene = $anno{ $gene2 };
				$rowcount2++;
				$tempgct .= $gene."\t".$string."\n";			    
			    }
			    else{next;}			    
			}
			else{ ###don't use gene anno
			    $gene = $gene2;
			    $rowcount1++;
			    $tempgct .= $gene."\t".$string."\n";
			}
			
		}
	}
}

#print "\n\n:::::::::TEST: rowcount $rowcount1 ::::::\n\n";
my $rowcountx = 0;
#print " Test rowcount is $rowcount2 or $rowcount \n";
if($unmatched){ $rowcountx = $rowcount2; } else { $rowcountx = $rowcount; }
if(!$aflag){ $rowcountx = $rowcount1; } #original counts while scanning ensid
print fout "#1.2\n" . $rowcountx ."\t" . $colcount."\n";
print LOG "__GCT final counts: $rowcountx rows, $colcount columns.\n";
print STDERR "__GCT final counts: $rowcountx rows, $colcount columns.\n";
#print "__:".$rowcountx."\t".$colcount."\n\n";
print fout $tempgct;

close expr;
print LOG "__Finished reading expression file.\n\n";
print STDERR "__Finished reading expression file.\n\n";

close fout;
print LOG "Finished generating gct file!\n";
print STDERR "Finished generating gct file!\n";

print LOG "Sample orderings: ".join(" ",@sampleOrders)."\n";
print STDERR "Sample orderings: ".join(" ",@sampleOrders)."\n";


#######Step 2: create CLS file
###note that current version is for categorical CLS and NOT continuous (time series/geneprofile)
print LOG "\nAssembling CLS header, groupings, and orderingso.\n";
print STDERR "\nAssembling CLS header, groupings, and orderings.\n";

my @cats = getUniqCats($sampleAnno, $groupName);
print LOG "__Categories generated: ".join(", ",@cats)."\n";
print STDERR "__Categories generated: ".join(", ",@cats)."\n";

my $clsheader = $colcount ." ". scalar @cats ." 1\n";
$clsheader .= "#".join(" ",@cats)."\n";
my $writestring = getCats($sampleAnno, $groupName, $sampleName, \@sampleOrders);
#print $clsheader.$writestring."\n";

open outcls, ">".$clsName or die "Error: Can't read $clsName $! \n";
printf outcls $clsheader.$writestring."\n";
close outcls;

print LOG "\nFinished generating cls file!\n\n";
print STDERR "\nFinished generating cls file!\n\n";
$now = localtime();
print LOG "Completion time: $now \n";
print STDERR "\nCompletion time: $now \n";

close LOG;


###
#####subroutines
sub getCats {
	#fun to get categories given key and filename
	#my @info = @_;
	#my @info = (@_);
	my $outstring;
	my $fname = $_[0]; ##name of config file
	my $group = $_[1]; ##group name to specify column
	my $sample = $_[2];
	my @sampleOrderx = @{$_[3]};
	my @temp;
	my $samplekey=0;
	my $tempkey = -1;
	my %mapper;
	open fconfig2, $fname or die "Error: Can't read $fname $! \n";
	
	$z = -1;
	my %temphash;

	###Read through configuration file to get mappings
	#print "\n\n___\nTesting sample orderings for config files:\n\n";
	#print join(" ",@sampleOrderx)."\n";
	while(<fconfig2>){
		$z++;
		chomp(my $str = $_);
		$str =~ s/[\r\n]+//;  ###Updated 4/15
		if($str ne ""){	
			my @vs = ();
			@vs = split(/\t/,$str);

			if($z == 0){     		
				#print "VS: ".join(",",@vs)."\n";
				$k=-1;
				foreach $val (values @vs){
					$k++;					
					#print "_next: $k \n";
					my $mp =""; $mp = $vs[$k];
					$temphash{ $mp } = $k; #$vs[$k]
					#print "__Next key: ($vs[$k]) value: $k \n";
				}
				#print join(", ",keys %temphash)."\n\n;;;;\n";
				$tempkey = $temphash{ $group } or die $!;
				#$samplekey = $temphash{ "Sample" } or die $!;
				#$samplekey = 0;
				#print "final key: ".$tempkey."\n";
			}
			else{ #now that tempkey is done, just grab the unique IDs
				#$outstring .= $vs[$tempkey]." ";
				#push(@temp, $vs[$tempkey]);
				#Previous version just assumed orderings were right
				###fix this: create hash to map sample->
				my $key=""; my $val="";
				$key = $vs[ $samplekey ];
				$val = $vs[ $tempkey ];				
				$mapper{ $key } = $val;
			}
		}
	}
	#print join(", ",keys %mapper)."\n";
	#print join(", ",values %mapper)."\n";
	
	foreach my $string (values @sampleOrderx){
		#print $string."\n";
		push(@temp, $mapper{$string});
	}
	$outstring = join("\t",@temp)."\n";
	close fconfig2;
	#chomp($outstring);
	$outstring =~ s/^[\s\t]+//;
	return( $outstring );
}


sub getUniqCats {
	#function to get unique categories based on selected column
	my @info = @_;
	my %outvals;
	my $fname = $info[0]; ##name of config file
	my $group = $info[1]; ##group name to specify column
	#print $fname."\n";
	open fconfig, $fname or die "ERROR: Can't read $fname $!";
	
	$z = -1;
	my %temphash;
        my @tempcatarr;
	my $tempkey = -1;
	while(<fconfig>){
		$z++;
		chomp(my $str = $_);
		$str =~ s/[\r\n]+//;
		if($str ne ""){	
			my @vs = split(/\t/,$str);
			if($z == 0){     		
				foreach my $k (keys @vs){
					$temphash{ $vs[$k] } = $k;
				}
				$tempkey = $temphash{ $group } or die $!;
			}
			else{ #now that tempkey is done, just grab the unique IDs
				#$outvals{$vs[$tempkey]}=1;
				push(@tempcatarr, $vs[$tempkey]);
			}
		}
	}
	close fconfig;
	###new: get list based on original ordering of appearance: modified from O'reilly
	%seen = ();
	my @uniqx = grep{ ! $seen{$_} ++ } @tempcatarr;

	return(@uniqx);
	#return( keys %outvals );
}



sub getAnno {
	my ($anno) = @_;
	$annofile = $anno; #[0];
	open in, $annofile or die $!;
	$x = 0;
	my %genehash={};
	while(<in>){
	    $x++; 
	    #if($x > 20){last;} ###test of missing values via restricting anno genes
		if($x>1){ #ignore header row
			#trim(
			chomp(my $string = $_);
			#);
			$string =~ s/[\r\n]+//;
			if($string ne ""){
				my @vals = split(/\t/,$string);
				my $v1 = ""; my $v2 = "";
				$v1 = $vals[0]; chomp($v2 = $vals[1]);
				$v2 =~ s/\s+//g;
				###trivial case of non-identified ensids reassigned 
				###   (e.g. no gene name)
				if($v2 eq ""){ } #update: don't add unmatched to hash  ####$v2 = $v1; }
				else{ $genehash{$v1} = $v2; }
			}
		}	
	}
	return %genehash;
}
