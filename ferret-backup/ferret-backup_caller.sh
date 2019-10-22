#!/bin/sh


#######
#Author: Jun Yin <jyin@sbpdiscovery.org>
#SBP Bioinformatics Core
#######

version="0.1"


#######
#Usage
#sh ferret-backup inputfolder location_at_ferret
#######

usage="

ferret-backup

version: $version

Usage: sh ferret-backup -i inputfolder -n nameofzippedfile -o location_at_ferret

Description: Zip the folder and backup to a ferret location

Parameters:

	-i      Input file or folder
	-n      Name of zipped file
	-o      Output folder
	
"


if [ $# -eq 0 ]; then 
	printf "$usage";
	exit 1; 
fi


#####
#Functions needed
#####

realpath() {
    path=`eval echo "$1"`
    folder=$(dirname "$path")
    echo $(cd "$folder"; pwd)/$(basename "$path"); 
}



#######
#Input/Output
#######

#receive options
while getopts ":i:o:n:" opt; do
  case ${opt} in
    i )
		input=$(realpath $OPTARG)
      ;;
    n )
		name=$OPTARG
      ;;	  
    o ) 
		output=$OPTARG
      ;;	  	  
    \? ) 
		printf "ERROR: Unknown options.\n\n$usage"
      ;;
	: ) printf "ERROR: Unknown options.\n\n$usage"
      ;;

  esac
done


#######
##start reading
#######


inputbasename=$(basename $input)

if [ -z "$name" ]
	then name=$inputbasename
fi

infolder=$(dirname "$input")


#log
logfile="$infolder/${inputbasename}_ferret-backup.log"

echo "sh " $0 $@ > $logfile
date >> $logfile

printf "sbptools ferret-backup $version running\n\n" 2>&1 | tee -a $logfile

printf "Reading $infolder for $input\n\n" 2>&1 | tee -a $logfile

#######
#Process
#######

if [ ${input: -7} != ".tar.gz" ]
	then 
		printf "Zipping $input into $input.tar.gz\n\n" 2>&1 | tee -a $logfile
		tar -czvf $name.tar.gz $input 2>&1 | tee -a $logfile
		
		printf "\nTransfer $name.tar.gz to ferret:$output\n\n" 2>&1 | tee -a $logfile
		
		echo "Running: echo \"put $name.tar.gz\"  | sftp ferret:$output" 2>&1 | tee -a $logfile
		echo "put $name.tar.gz" | sftp ferret:$output 2>&1 | tee -a $logfile
		
else
	printf "\nTransfer $input to ferret:$output\n\n" 2>&1 | tee -a $logfile
	echo "Running: echo \"put $input\" | sftp ferret:$output" 2>&1 | tee -a $logfile
	echo "put $input" | sftp ferret:$output 2>&1 | tee -a $logfile
fi

printf "\n\nDone.Running log saved in $logfile\n\n" 2>&1 | tee -a $logfile


