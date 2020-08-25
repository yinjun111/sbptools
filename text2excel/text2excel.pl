#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Excel::Writer::XLSX;

#Andrew P. Hodges, Ph.D., Jun Yin, Ph.D.
#Copyright, Sanford Burnham Prebys Medical Discovery Institute

#Purpose: generate excel tables from input list of txt files
##


########
#Prerequisites
########

#Text to excel file:


########
#Interface
########


my $version="0.3";

#v0.3 by Jun
#supports wildcards. fix a bug for no --names 


#v0.2 by Jun
#add new theme to wrap text, filter and freeze panel
#add log file and timestamp


#v0.92 to be implemented, 
#1) command line in linux (getopt/long included above & below)
#2) accept multiple input text files into different tabs
###    optional: input to name those tabs differently
#3) Prevent automatic conversion (e.g. specify columns that should be txt and not general
#Optional:  
#####first row bold
#####first row has filter function
#####freeze the first row
#####Color theme selection
####example  $text2excel -i file1,file2,file3 -n name1,name2,name3 -o merged.xlsx


my $usage="
text2excel
version: $version\n

Description: Perl script to generate compile xlsx file from individual text files. The script will prevent automatic format changes in Excel for gene names, special texts, and combine different texts in to tabs.\n

Usage: 

    #Simple usage
    perl text2excel.pl -i file1.txt,file2.txt -o result.xlsx

    #Set sheet names, and choose different themes
    perl text2excel.pl -i file1.txt,file2.txt -n ShName1,Shname2 -t 1,2 --theme theme2 -o result.xlsx\n

    #Use wildcards for input files
    perl text2excel.pl -i \"*.txt\" --theme theme2 -o result.xlsx\n


Parameters:
	--in|-i           input file(s) separated by \",\", support wildcards
                           if wildcards are used, --names won't be supported.
	--out|-o          output file

	--names|-n        Sheet names
	--txt|-t          column number starting from 0 that should be txt not general separated by \",\"

	--boldfrow|-bfr   bold first row [F]
	--boldfcol|-bfc   bold first column [F]
	--color|-c        color theme to use

    #themes to be used
	--theme           theme1, by AH 
                      theme2, by JY, adding wrap text, filter etc. Now default option. [theme2]
                      theme0, don't change format

	--delim|-d        default is tab-delimited; use '' for other entries [\\t]
	--verbose|-v      Verbose\n	
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
my $outfile;
my $names="";
my $txt=0;
my $boldfrow="F";
my $boldfcol="F";
my $verbose;
my $color = "";
my $theme="theme2";

GetOptions(
	
	"in|i=s" => \$infiles,
	"out|o=s" => \$outfile,
	"names|n=s" => \$names,
	"theme=s" => \$theme,	
	"txt|t=s" => \$txt,
	"boldfrow|bfr" => \$boldfrow,
	"boldfcol|bfc" => \$boldfcol,
	"color|c=s" => \$color,
	"verbose|v" => \$verbose,
);

my $logfile=$outfile;
$logfile=~s/\.\w+$/_text2excel.log/;

####First check if the file exists & it is xlsx


######write log file
open(LOG, ">$logfile") || die "Error write $logfile. $!";

my $now=current_time();

print LOG "perl $0 $params\n\n";
print LOG "Start time: $now\n\n";
print LOG "Current version: $version\n\n";

print LOG "\n";


#####Parse params

my @ins;
my $out = $outfile;
my @names;

#deal with wildcard support
if($infiles=~/\*/) {
	#double check names for wildcards
	if(defined $names && length($names)>0) {
		print STDERR "ERROR:--names is not supported, when wildcard is used in --infiles $infiles.\n";
		exit;
	}
	
	foreach my $file (split(",",$infiles)) {
		my @files=glob($file);
		push @ins,@files;
	}
}
else {
	@ins = split(/,/,$infiles);
	@names=split(",",$names);
	
	if(@names>0) {
		if(scalar(@ins) != scalar(@names)) {
			print STDERR "ERROR: --in $infiles (",scalar(@ins),") and --names $names (",scalar(@names),") do not match.\n";
			exit;
		}
	}
	
}

####Next create the excel object
print "- Generating excel object: \n";
my $excel = Excel::Writer::XLSX->new( $outfile ) or die $!;
$excel->set_properties(
   #title => $title,
   author => "BI Shared Resource",
   manager => "Andrew P. Hodges, Ph.D.",
   comments => "Auto-generated excel file from script.",
);
#$excel->set_custom_property('Date generated',date(),'number');


my @textcols = split(/,/,$txt);
my %textcols = map { "x_".$_ => 1 } @textcols;
my $bfrow = $boldfrow;
my $bfcol = $boldfcol;
my $col = $color;

####Also create format blocks for the excel components
print "- Setting up excel formatting: \n";

my $formatHeader = $excel->add_format();

#choice to use different theme
if($theme eq "theme1") {
	$formatHeader->set_bold();
	$formatHeader->set_color('red');
	$formatHeader->set_align('center');
}
elsif($theme eq "theme2") {
	$formatHeader->set_bold();
	$formatHeader->set_color('black');
	$formatHeader->set_align('center');
	$formatHeader->set_text_wrap();
}

#);

my $formatColtxt = $excel->add_format(
	type => 'text',
	align => 'left',
);


#####MAIN CODE
###for each file in list,
###   1) get name, remove extension, pass over as the sheet name (if sheet name not specified
###	  2) add worksheets for each file.

my @worksheets;
foreach my $R (keys @ins){   #index
	
	my $filex = $ins[$R];
	
	print "- - Working on $filex\n";
	
	#for now, just open the file & use names for sheets
	my $sname = $filex;
	if(exists($names[$R])){ $sname = $names[$R]; }
	#print length($sname)."\n";
	
	#update: check length & adjust to 30 if too long.
	my $len = length($sname);
	if($len > 31){$sname = substr($sname,0,30);}
	#print "\nUsing: ".$sname."\n";
	$worksheets[$R] = $excel->add_worksheet($sname);
	print "- - Created new sheet: ".$sname."\n";
	#my $tempsheet = $worksheets[$R];
	
	###apply operations to the current worksheet
	my $i = -1; #row index
	my $maxcol=0;
	open IN, $filex or die "ERROR:$filex not found.$!";
	while(<IN>){
		$i++;
		chomp(my $string = $_);
		#first see if data is 
		if($string ne ""){
			my @row = split(/\t/,$string);
			my $count = @row;
			
			#find the max column number
			if($count>$maxcol) {
				$maxcol=$count;
			}
			
			#xl_rowcol_to_cell( 1, 2 )
			if($i == 0){ 
				my $arrayref = \@row;
				$worksheets[$R]->write_row(0,0,$arrayref,$formatHeader); }
			else{  
					#use default formatting here with 
					##look @ each column... if column is in the list, set as txt
					#$formatColtxt
				for my $j (0 .. $count){
					if(exists($textcols{"x_".$j})){ 
						#print "Success in column $j \n";
						$worksheets[$R]->write($i,$j,$row[$j],$formatColtxt);
					}
					else{
						$worksheets[$R]->write($i,$j,$row[$j]);
					}
			
				}
			}
		}
		else{ $worksheets[$R]->write("\n");}
		
	}
	close IN;	
	
	#additional step to change theme
	if($theme eq "theme2") {
			#freeze first row
			$worksheets[$R]->freeze_panes( 1, 0 );
			$worksheets[$R]->autofilter( 0, 0, $i, $maxcol-1 );
	}
	
	###finalize worksheet
	
}



###End program
$excel->close() or die "Error Closing File: $! \n";
print "- Successfully closed & saved excel file.\n";

close LOG;

####Helpful functions:
#use Excel::Writer::XLSX::Utility;
#( $row, $col ) = xl_cell_to_rowcol( 'C2' );    # (1, 2)
#$str           = xl_rowcol_to_cell( 1, 2 ); 



#$worksheet->write_row( 'A1', $array_ref );    # Write a row of data
#$worksheet->write(     'A1', $array_ref );    # Same thing


####Font example:
# my %font = (
    # font  => 'Calibri',
    # size  => 12,
    # color => 'blue',
    # bold  => 1,
# );

# my %shading = (
    # bg_color => 'green',
    # pattern  => 1,
#);

# $fields <- 'A1:A4'
# $excel->conditional_formatting( $fields,
    # {
        # type     => 'text',  ###this is used to bypass the general type for genes etc.
        #criteria => 'containing',
        #value    => 'foo',
        #format   => $format,
    # }
# );


########
#Functions
########

sub current_time {
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
	my $now = sprintf("%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec);
	return $now;
}