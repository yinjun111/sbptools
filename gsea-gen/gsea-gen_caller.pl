#!/usr/bin/perl
use strict;
use Cwd qw(abs_path);
use Getopt::Long;
use List::MoreUtils qw(uniq);
use File::Basename qw(basename dirname);


#Andrew Hodges, PhD
#BI Shared Resource, SBP
#11/6/2019
#Update: 4/22/2020
#(C)2019-2020, SBP
#Jun Yin, PhD
#BI Shared Resource, SBP
#7/20/2020
#Update: 4/22/2020
#(C)2019-2020, SBP


########
#Prerequisites
########


########
#Interface
########


my $version="1.1";

#v1.0a, perform GSEA analysis after cls and gct are generated.
#v1.1 versioning

my $usage="
gseagen
version: $version\n
Usage:  sbptools gsea-gen -e ./data/gene.results.merged.tpm.txt -o resultsfolder -s ./data/configV1.txt -n Group -t Human.B38.Ensembl84 -d h.all.v7.1\n

Description: Perl script to generate gct and cls files for GSEA analysis based on gene expression matrix
The script needs --comparisons to run GSEA. If no --comparisons is provided, it will just generate .gct and .cls files.


Parameters:

	--expression|-e   Input file of merged expression data
	--out|-o          Output folder
	
	--sampleAnno|-s   Sample annotation 
	--groupName|-n    Name of group column for CLS generation

	--tx|-t           Transcriptome
                        Currently support Human.B38.Ensembl84, Mouse.B38.Ensembl84
	--chip            Use MSigDB chip file, e.g. Human_ENSEMBL_Gene_MSigDB.7.1 (optional)
	--geneAnno|-ga    Gene annotation (optional)

	--comparisons|-c  Comparisons, one pair comparison is recommended(optional)

	--db|-d           MSigDB database to be used
                       Recommended:
                        h.all.v7.1
                        c5.bp.v7.1
						
                       Others include: 
                       c1.all.v7.1,c2.all.v7.1,c2.cgp.v7.1,c2.cp.biocarta.v7.1,c2.cp.kegg.v7.1,c2.cp.pid.v7.1,c2.cp.reactome.v7.1,c2.cp.v7.1,c3.all.v7.1,c3.mir.mirdb.v7.1,c3.mir.mir_legacy.v7.1,c3.mir.v7.1,c3.tft.gtrd.v7.1,c3.tft.tft_legacy.v7.1,c3.tft.v7.1,c4.all.v7.1,c4.cgn.v7.1,c4.cm.v7.1,c5.all.v7.1,c5.bp.v7.1,c5.cc.v7.1,c5.mf.v7.1,c6.all.v7.1,c7.all.v7.1,h.all.v7.1

    --runmode|-r      Where to run the scripts, local, system, server or none [none]
                                  none, only prints scripts
                                  local, run the scripts using parallel in local workstation(Falco)	
                                  system, run the scripts without paralleling (Falco)
                                  server, run the scripts in cluster (Firefly)

    --verbose|-v      Verbose

";

#	--gctfile|-g      output gct file name
#	--clsfile|-c      output cls file name





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
my $outputfolder;
my $tx;
my $chipfile;
my $db;

my $comparisons;
my $sampleName="Sample";
#my $dirExprF;
#my $dirSanno; # = "";
#my $dirGanno;
my $verbose=0;
my $groupName;

my $jobs=5;
my $runmode="none";

my $dev=0;

####re-add here
GetOptions(	
    "expression|e=s" => \$expression,
	"out|o=s" => \$outputfolder,
    "sampleAnno|s=s" => \$sampleAnno,
    "geneAnno|ga=s" => \$geneAnno,
    "groupName|n=s" => \$groupName,
	"tx|t=s" => \$tx,
	"chip=s" => \$chipfile,	
	"comparisons|c=s" => \$comparisons,
	"db|d=s" => \$db,
    "sampleName|r=s" => \$sampleName, #sample column ID used in config file	  
	"dev" => \$dev,	
	"runmode|r=s" => \$runmode,	
    "verbose|v" => \$verbose,
);


#    "gctfile|g=s" => \$gctfile,
#    "clsfile|c=s" => \$clsfile,

########
#Prerequisites
########

#my $sbptools=locate_cmd("sbptools","/usr/bin/sbptools");
#my $motiffinder="$sbptools motif-finder"; #improve compatibility

my $sbptoolsfolder="/apps/sbptools/";

#Dev version
if($dev) {
	$sbptoolsfolder="/home/jyin/Projects/Pipeline/sbptools/";
}
else {
	#the tools called will be within the same folder of the script
	$sbptoolsfolder=get_parent_folder(abs_path(dirname($0)));
}


my $parallel_job="perl $sbptoolsfolder/parallel-job/parallel-job_caller.pl"; 
my $gsea_cli="/apps/GSEA_Linux_4.0.3/gsea-cli.sh GSEA";




########
#Program begins
########

#Create folders

if(!-e $outputfolder) {
	mkdir($outputfolder);
}

$outputfolder = abs_path($outputfolder);

my $outputfoldername=basename($outputfolder);


my $logfile="$outputfolder/gsea-gen_run.log";

my $gctfile="$outputfolder/$outputfoldername.gct";
my $clsfile="$outputfolder/$outputfoldername.cls";

my $scriptfile1="$outputfolder/gsea-gen_run.sh";
my $scriptlocalrun="$outputfolder/gsea-gen_local_submission.sh";
my $scriptclusterrun="$outputfolder/gsea-gen_cluster_submission.sh";


#write log file
open(LOG, ">$logfile") || die "Error writing into $logfile. $!";

print LOG "gsea-gen log file: \n";
my $now = localtime();
my $timestamp=build_timestamp($now,"long");

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print STDERR "Start time: $now\n\n" if $verbose;
print LOG "Curent version: $version\n\n";
print LOG "\n";



#####
#next get annotation
#####

print STDERR "gsea-gen version $version running.\n\n";
print LOG "gsea-gen version $version running.\n\n";

my %anno;
my $aflag = 0; #aflag is whether or not to use the gene annotation file to label rows
print LOG "Annotation: ";
print STDERR "Annotation: " if $verbose;
if($geneAnno ne ""){ %anno = getAnno( $geneAnno); $aflag=1; print LOG "$geneAnno";}else{ print LOG "No gene annotation file selected.";}
print LOG "\n";
my $rowcount;
my $colcount;


#$rowcount = 47643;
my $expr = "wc -l $expression";
$rowcount = `$expr` - 1;
print LOG "Number of rows: ".$rowcount."\n";
print STDERR "Number of rows: ".$rowcount."\n" if $verbose;

#$colcount = 20;
my $expr2 = "head -n2 $expression | tail -n1 | wc -w";
$colcount = `$expr2` - 1; #assume 1 gene info columns
print LOG "Number of columns: ".$colcount."\n";
print STDERR "Number of columns: ".$colcount."\n" if $verbose;




######
#create gct file
######

#Next open main filehandles for data and output gct file
open fout, ">".$gctfile or die "Error: Can't write $gctfile $! \n";
print LOG "\nOpened GCT file for writing.\n";
print STDERR "\nOpened GCT file for writing.\n" if $verbose;

open expr, $expression or die "Error: Can't read $expression $! \n";
printf fout "#1.2\n".$rowcount."\t".$colcount."\n";
print LOG "Opened expression file for reading.\n";
print STDERR "Opened expression file for reading.\n" if $verbose;


my @sampleOrders;


####this can be added into a subroutine later
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
			printf fout $string."\n";
		}
		else{ ###data non-blank rows
			###first get the gene using annotation hash
			my @vals = split(/\t/,$string);
			my $gene2 = $vals[0];
			my $gene = ""; 
			if($aflag){$gene = $anno{ $gene2 };}
			else{$gene = $gene2;}
			#print $gene."\n";
			###update string
			###write to file
			printf fout $gene."\t".$string."\n";
		}
	}
}
close expr;
print LOG "__Finished reading expression file.\n\n";
print STDERR "__Finished reading expression file.\n\n" if $verbose;

close fout;
print LOG "Finished generating gct file!\n";
print STDERR "Finished generating gct file!\n" if $verbose;

print LOG "Sample orderings: ".join(" ",@sampleOrders)."\n";
print STDERR "Sample orderings: ".join(" ",@sampleOrders)."\n" if $verbose;


######
#create CLS file
######

###note that current version is for categorical CLS and NOT continuous (time series/geneprofile)
print LOG "\nAssembling CLS header, groupings, and orderingso.\n";
print STDERR "\nAssembling CLS header, groupings, and orderings.\n" if $verbose;

my @cats = getUniqCats($sampleAnno, $groupName);
print LOG "__Categories generated: ".join(", ",@cats)."\n";
print STDERR "__Categories generated: ".join(", ",@cats)."\n" if $verbose;

my %cats_hash=map {$_,1} @cats;

my $clsheader = $colcount ." ". scalar @cats ." 1\n";
$clsheader .= "#".join(" ",@cats)."\n";
my $writestring = getCats($sampleAnno, $groupName, $sampleName, \@sampleOrders);
#print $clsheader.$writestring."\n";

open outcls, ">".$clsfile or die "Error: Can't read $clsfile $! \n";
printf outcls $clsheader.$writestring."\n";
close outcls;

print LOG "\nFinished generating cls file!\n\n";
print STDERR "\nFinished generating cls file!\n\n" if $verbose;


######
#Read comparisons
######

my %comparisons_all;
my %comp_groups;

if(defined $comparisons && length($comparisons)>0) {
	foreach my $comp (split(",",$comparisons)) {
		if($comp=~/(.+)_versus_(.+)/) {
			if(!defined $cats_hash{$1}) {
				print STDERR "ERROR:Group $1 not defined.\n\n";
				print LOG "ERROR:Group $1 not defined.\n\n";
				exit;
			}
			if(!defined $cats_hash{$2}) {
				print STDERR "ERROR:Group $1 not defined.\n\n";
				print LOG "ERROR:Group $1 not defined.\n\n";
				exit;
			}
		}
		else {
			print STDERR "ERROR:Comparison $comp is not in the right format. Please use group1_versus_group2 for comparison name.\n\n";
			print LOG "ERROR:Comparison $comp is not in the right format. Please use group1_versus_group2 for comparison name.\n\n";
		}
		$comparisons_all{$comp}++;
	}

	print STDERR scalar(keys %comparisons_all), " comparison(s) are needed:\n";
	print STDERR join("\n",sort keys %comparisons_all),"\n";

	print LOG scalar(keys %comparisons_all), " comparison(s) are needed:\n";
	print LOG join("\n",sort keys %comparisons_all),"\n";
}
else {

	print STDERR "No comparisons are defined. Skip GSEA run. Only .gct and .cls files are generated.\n";
	print LOG "No comparisons are defined. Skip GSEA run. Only .gct and .cls files are generated.\n";
	
	exit; #exit here before generating GSEA scripts
}


######
#Run GSEA
######

#MSigdb, gmt annotation files 
my $db2file="/data/jyin/Databases/GSEA/msigdb_v7.1_files_to_download_locally/msigdb_v7.1_GMTs/"."$db.symbols.gmt";

if(!-e $db2file) {
	print STDERR "ERROR:--db $db not supported. $db2file was not found.\n\n";
	exit;
}

#chip, gene annotation
my $chipfile;
if($tx eq "Human.B38.Ensembl84") {
	$chipfile="/data/jyin/Databases/GSEA/msigdb_v7.1_chip_files_to_download_locally/Human_ENSEMBL_Gene_MSigDB.7.1.chip";
}
elsif ($tx eq "Mouse.B38.Ensembl84") {
	$chipfile="/data/jyin/Databases/GSEA/msigdb_v7.1_chip_files_to_download_locally/Mouse_ENSEMBL_Gene_ID_to_Human_Orthologs_MSigDB.7.1.chip";
}
else {
	print STDERR "ERROR: --tx $tx not recognized.\n\n";
	exit;
}

if(defined $comparisons && length($comparisons)>0) {
	
	open(OUT,">$scriptfile1") || die $!;
	foreach my $comparison (sort keys %comparisons_all) {
		print OUT "$gsea_cli -res $gctfile -cls $clsfile#$comparison -gmx $db2file -chip $chipfile -out $outputfolder -collapse Collapse -mode Max_probe -norm meandiv -nperm 1000 -permute phenotype -rnd_type no_balance -scoring_scheme weighted -rpt_label my_analysis -metric Signal2Noise -sort real -order descending  -create_gcts false -create_svgs false -include_only_symbols true -make_sets true -median false -num 100 -plot_top_x 30 -rnd_seed timestamp -save_rnd_lists false -set_max 500 -set_min 10 -zip_report false\n";
	}
	close OUT;

}



#$now = localtime();
#print LOG "Completion time: $now \n";
#print STDERR "\nCompletion time: $now \n";

#close LOG;



#######
#Run mode
#######

my $jobnumber=0;
my $jobname="gsea-gen-$timestamp";

if($jobs eq "auto") {
	$jobnumber=0;
}
else {
	$jobnumber=$jobs;
}

my $localcommand="screen -S $jobname -dm bash -c \"source ~/.bashrc;cat $scriptfile1 | parallel -j $jobnumber;\"";


if($runmode eq "none") {
	print STDERR "\nTo run locally, in shell type: $localcommand\n\n";
	print LOG "\nTo run locally, in shell type: $localcommand\n\n";
	
	print STDERR "\nTo run in cluster, in shell type:$parallel_job -i $scriptfile1 -o $outputfolder/scripts/ -r\n\n";
	print LOG "\nTo run in cluster, in shell type:$parallel_job -i $scriptfile1 -o $outputfolder/scripts/ -r\n\n";
}
elsif($runmode eq "local") {
	#local mode
	
	#need to replace with "sbptools queuejob" later

	system($localcommand);
	print LOG "$localcommand;\n\n";

	print STDERR "Starting local paralleled processing using $jobnumber tasks. To monitor process, use \"screen -r $jobname\".\n\n";
	print LOG "Starting local paralleled processing using $jobnumber tasks. To monitor process, use \"screen -r $jobname\".\n\n";
	
}
if($runmode eq "system") {
	print STDERR "\nRunning locally without parallel. sh $scriptfile1.\n\n";
	print LOG "\nRunning locally without parallel. sh $scriptfile1.\n\n";
	
	system("sh $scriptfile1");
}
elsif($runmode eq "cluster") {
	#cluster mode	
	print STDERR "\nStart running in cluster using:$parallel_job -i $scriptfile1 -o $outputfolder/scripts/ -r\n\n";
	system("$parallel_job -i $scriptfile1 -o $outputfolder/scripts/ -r");
}


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
	
	my $z = -1;
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
				my $k=-1;
				foreach my $val (values @vs){
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
	chomp($outstring);
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
	
	my $z = -1;
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
	my %seen = ();
	my @uniqx = grep{ ! $seen{$_} ++ } @tempcatarr;

	return(@uniqx);
	#return( keys %outvals );
}



sub getAnno {
	my ($anno) = @_;
	my $annofile = $anno; #[0];
	open in, $annofile or die $!;
	my $x = 0;
	my %genehash;
	while(<in>){
		$x++;
		if($x>1){ #ignore header row
			#trim(
			chomp(my $string = $_);
			#);
			$string =~ s/[\r\n]+//;
			if($string ne ""){
				my @vals = split(/\t/,$string);
				my $v1 = ""; my $v2 = "";
				$v1 = $vals[0]; $v2 = $vals[1];
				###trivial case of non-identified ensids reassigned 
				###   (e.g. no gene name)
				if($v2 eq ""){ $v2 = $v1; }
				$genehash{$v1} = $v2;
			}
		}	
	}
	return %genehash;
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

sub get_parent_folder {
	my $dir=shift @_;
	
	if($dir=~/^(.+\/)[^\/]+\/?/) {
		return $1;
	}
}

