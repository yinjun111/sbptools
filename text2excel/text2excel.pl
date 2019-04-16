#!/usr/bin/perl
#Andrew P. Hodges, Ph.D.
#5 April 2019
#Purpose: generate excel tables from input list of txt files
##



#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Excel::Writer::XLSX;

#CutAdapt+FASTQC+RSEM+STAR


########
#Prerequisites
########

#Text to excel file:


########
#Interface
########


my $version="0.1";

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
Usage: perl text2excel.pl -i file1.txt,file2.txt -n ShName1,Shname2 -t 1,2 -bfr -bfc -o result.txt\n
Description: Perl script to generate compile xlsx file from individual text files.\n
Parameters:
	--in|-i           input file(s) separated by \",\"
	--out|-o          output file
	--names|-n        Sheet names
	--txt|-t          column number starting from 0 that should be txt not general separated by \",\"
	--boldfrow|-bfr   bold first row [F]
	--boldfcol|-bfc   bold first column [F]
	--color|-c        color theme to use
	--delim|-d        default is tab-delimited; use '' for other entries
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
my $names;
my $txt;
my $boldfrow="F";
my $boldfcol="F";
my $verbose;
my $color = "";

GetOptions(
	
	"in|i=s" => \$infiles,
	"out|o=s" => \$outfile,
	"names|n=s" => \$names,
	"txt|t=s" => \$txt,
	"boldfrow|bfr" => \$boldfrow,
	"boldfcol|bfc" => \$boldfcol,
	"color|c=s" => \$color,
	"verbose|v" => \$verbose,
);

#my $logfile=$outfile;
#$logfile=~s/\.txt/_mergefiles.log/;

####First check if the file exists & it is xlsx


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

#####Parse params
print "- Setting up excel formatting: \n";
my @ins = split(/,/,$infiles);
my $out = $outfile;
my @names = split(/,/,$names);
my @textcols = split(/,/,$txt);
my %textcols = map { "x_".$_ => 1 } @textcols;
my $bfrow = $boldfrow;
my $bfcol = $boldfcol;
my $col = $color;

####Also create format blocks for the excel components
my $formatHeader = $excel->add_format();
$formatHeader->set_bold();
$formatHeader->set_color('red');
$formatHeader->set_align('center');
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
	open in, $filex or die $!;
	while(<in>){
		$i++;
		chomp(my $string = $_);
		#first see if data is 
		if($string ne ""){
			my @row = split(/\t/,$string);
			my $count = @row;
			
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
	close in;	
	###finalize worksheet
}



###End program
$excel->close() or die "Error Closing File: $! \n";
print "- Successfully closed & saved excel file.\n";



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


