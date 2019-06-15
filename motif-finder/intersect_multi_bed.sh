#!/bin/sh

#cal bedtools intersect to intersect multiple bed files

#now only supports:
#/apps/bedtools2-2.26.0/bin/bedtools intersect -a $infile -b $dbfile -wb > $outfile


#version 0.1
#sh intersect_multi_bed.sh in.bed "db*.bed" out.txt


infile=$1
dbfiles=$2
outfile=$3

echo -e "\nProgram Starts:"
echo "infile:" $infile 
echo "dbfiles:" $dbfiles 
echo -e "outfile:" $outfile "\n"

#filenames
infilename=$(basename $infile)
outfilename=$(basename $outfile)

#make outfile
echo "" > $outfile

#tempfolder
tempfolder=${outfile/.txt/_temp}
mkdir $tempfolder

runfile=$tempfolder/${outfilename/.txt/_intersect_run.sh}


for dbfile in $dbfiles
	do
		#echo $dbfile
		dbfilename=$(basename $dbfile)
		echo "/apps/bedtools2-2.26.0/bin/bedtools intersect -a $infile -b $dbfile -wb > $tempfolder/${infilename/.bed/}_intersect_${dbfilename/.bed/}.txt" >> $runfile
done
		

#use 6 processes for intersect bed
echo -e "Intersect program running.\n\n"
cat $runfile | parallel -j 6

#merge all
cat $tempfolder/${infilename/.bed/}_intersect_* > $outfile

#remove temp files
rm -R $tempfolder
