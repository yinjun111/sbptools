#!/bin/sh

#######
#Author: Jun Yin <jyin@sbpdiscovery.org>
#SBP Bioinformatics Core
#######

version="2.1"

#version 1.1 add option to control params
#version 1.2 add reformating snpeff and rename snpEff_summary.genes
#version 1.3 remove huge alignment and index files to free up space
#v1.3a, add -f to tabix
#v1.4, compatibility in Firefly
#v1.41 fix samtools bug
#v1.5 start from bam file
#v2, add gnomad common SNP filter and DP filter
#v2.1, add TMB summary. New snpeff. Remove more files

#######
#Usage
#sh gatk3_rnaseq_variant_version.sh yourfastq.fastq.gz outputfolder
#######

usage="

rnaseq-var

version: $version

Usage: sbptools rnaseq-var -i yourfastq.fastq.gz -o outputfolder -s species

Description: RNA-Seq variant calling pipeline built using GATK3. This pipeline is customized to call and annotate variants from RNA-Seq data for human/mouse using B38 annotation. You only need to provde fastq file and an output folder.

Currently SE sequencing only

Parameters:

	-i      Input fastq file, fastq file name must end with fastq.gz
	-b      Start from bam files instead of fastq.
	-o      Output folder
	-s      Species, human or mouse

	-y      Dry run. Only write commands and log files

	Optional
	-d      Whethter to dedup [T]
	-f      DP filter [5.0]
	-c      gnomad v2 common SNP filter 0.01 or 0.05 [0.01]	
	-k      Keep alignment data [F]
	
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

dedup=T
species=human
dpfilter=5.0
csfilter=0.01
keepalignment=F
dryrun=F

#receive options
while getopts "i:b:o:s:d:f:y:" opt; do
  case ${opt} in
    i )
		fastqfile=$(realpath $OPTARG)
		filename=$(basename $fastqfile)
		filename=${filename/.fastq.gz/}
      ;;
    b )
		bamfile=$(realpath $OPTARG)
		filename=$(basename $bamfile)
		filename=${filename/.bam/}
      ;;
    o ) 
		outfolder=$(realpath $OPTARG)
      ;;
    s ) 
		species=$OPTARG
      ;;
    d ) 
		dedup=$OPTARG
      ;;
    f ) 
		dpfilter=$OPTARG
      ;;
    c ) 
		csfilter=$OPTARG
      ;;
    k ) 
		keepalignment=$OPTARG
      ;;
    y ) 
		dryrun=T
      ;;		  
    \? ) 
		printf "ERROR: Unknown options.\n\n$usage"
      ;;
	: ) printf "ERROR: Unknown options.\n\n$usage"
      ;;

  esac
done


#fastqfile=$(realpath $1)
#outfolder=$(realpath $2)
#species=$3 #human or mouse

if [ -f "$fastqfile" ] || [ -f "$bamfile" ];then
	#logfile=$outfolder/$filename.run.log
	logfile=$outfolder/rnaseq-var_run.log
else
	echo "ERROR:No input file was found."
fi
	
#run file for recording command
if [ -f "$fastqfile" ] || [ -f "$bamfile" ];then
	runfile=$outfolder/$filename.run.sh
else
	echo "ERROR:No input file was found."
fi



#log 
mkdir -p $outfolder/

echo "sh " $0 $@ > $logfile
date >> $logfile


#######
#Prerequisites
#######


#databases

#genomedir=/Users/diazmeco/RNASeqVar/db/mm10/Mouse.B38.Ensembl84_STAR
#genomefasta=/Users/diazmeco/RNASeqVar/db/mm10/Mus_musculus.GRCm38.dna.primary_assembly_ucsc.fa
#knownsnp=/Users/diazmeco/RNASeqVar/db/mm10/mgp.v6.merged.norm.snp.indels.sfiltered.short.vcf.gz

if [ "$species" == "human" ]
then
	genomeversion=hg38
	genomedir=/data/jyin/Databases/Genomes/Human/hg38/Human.B38.Ensembl84_STAR
	genomefasta=/data/jyin/Databases/Genomes/Human/hg38/Homo_sapiens.GRCh38.dna.primary_assembly_ucsc.fa
	knownsnp=/data/jyin/Databases/SNP/00-All.vcf.gz
	
	if [ "$csfilter" == "0.01" ]
	then
		gnomadcommonsnp=/data/jyin/Databases/gnomad/fromGATK4/af-only-gnomad.hg38.common001.vcf
	else
		gnomadcommonsnp=/data/jyin/Databases/gnomad/fromGATK4/af-only-gnomad.hg38.common005.vcf
	fi
	
else
	genomeversion=mm10
	genomedir=/data/jyin/Databases/Genomes/Mouse/mm10/Mouse.B38.Ensembl84_STAR
	genomefasta=/data/jyin/Databases/Genomes/Mouse/mm10/Mus_musculus.GRCm38.dna.primary_assembly_ucsc.fa
	knownsnp=/data/jyin/Databases/SNP/mgp.v6.merged.norm.snp.indels.sfiltered.short.vcf.gz
fi	

#programs

#gatk3folder=/home/jyin/Programs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef
#bgzip=/home/jyin/Programs/htslib-1.9/bgzip
#tabix=/home/jyin/Programs/htslib-1.9/tabix
#snpeff=/home/jyin/Programs/snpEff/snpEff.jar
#snpsift=/home/jyin/Programs/snpEff/SnpSift.jar
#reformatsnpeffvcf=/home/jyin/Projects/Pipeline/sbptools/rnaseq-var/reformat_snpeff_vcf.pl
#rnaseq_var_filter=/home/jyin/Projects/Pipeline/sbptools/rnaseq-var/rnaseq-var_filter.pl

#java=/usr/bin/java						  
java=java
star=/apps/STAR-master/bin/Linux_x86_64/STAR
gatk3folder=/apps/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef
bgzip=/apps/htslib-1.9/bgzip
tabix=/apps/htslib-1.9/tabix
snpeff=/apps/snpEff-5/snpEff.jar
snpsift=/apps/snpEff-5/SnpSift.jar
reformatsnpeffvcf=/apps/sbptools/rnaseq-var/reformat_snpeff_vcf.pl
samtools=/apps/bin/samtools
rnaseq_var_filter=/apps/sbptools/rnaseq-var/rnaseq-var_filter.pl

#####
#Program starts
#####

printf "\nrnaseq-var version $version running\n\n" | tee -a $logfile $runfile

#######
#Step 0
#######
#samtools faidx $genomefasta

#java -jar $gatk3folder/picard.jar CreateSequenceDictionary R= $genomefasta O= ${genomefasta/.fa/.dict}


#######
#Step 1, 2-pass algnment
#######


if [ -f "$fastqfile" ];then

	#1st-pass
	printf "#1-pass alignment\n" | tee -a  $logfile $runfile

	if [ "$dryrun" == "F" ];then
		mkdir -p $outfolder/alignment/1pass_alignment
		
		$star --genomeDir $genomedir --readFilesIn $fastqfile --outSAMtype BAM Unsorted --readFilesCommand 'gunzip -c' --outFileNamePrefix $outfolder/alignment/1pass_alignment/$filename >> $logfile 2>&1
	fi
	
	printf "mkdir -p $outfolder/alignment/1pass_alignment;$star --genomeDir $genomedir --readFilesIn $fastqfile --outSAMtype BAM Unsorted --readFilesCommand 'gunzip -c' --outFileNamePrefix $outfolder/alignment/1pass_alignment/$filename >> $logfile 2>&1\n" | tee -a $runfile
	

	#2-pass #need to use a different alignment output folder ...
	printf "\n\n#2-pass alignment, genomeGenerate\n" | tee -a  $logfile $runfile

	
	if [ "$dryrun" == "F" ];then
		mkdir -p $outfolder/alignment/2pass_genomeDir
		$star --runMode genomeGenerate --genomeDir $outfolder/alignment/2pass_genomeDir --outFileNamePrefix $outfolder/alignment/2pass_genomeDir/$filename --genomeFastaFiles $genomefasta --sjdbFileChrStartEnd $outfolder/alignment/1pass_alignment/${filename}SJ.out.tab --sjdbOverhang 75 --runThreadN 4 >> $logfile 2>&1
	fi

	printf "mkdir -p $outfolder/alignment/2pass_genomeDir;$star --runMode genomeGenerate --genomeDir $outfolder/alignment/2pass_genomeDir --outFileNamePrefix $outfolder/alignment/2pass_genomeDir/$filename --genomeFastaFiles $genomefasta --sjdbFileChrStartEnd $outfolder/alignment/1pass_alignment/$filenameSJ.out.tab --sjdbOverhang 75 --runThreadN 4 >> $logfile 2>&1\n" | tee -a $runfile
	
	#2-pass alignment #convert output to bam
	printf "\n\n#2-pass alignment, align\n" | tee -a  $logfile $runfile
	if [ "$dryrun" == "F" ];then
		mkdir -p $outfolder/alignment/2pass_alignment
		$star --genomeDir $outfolder/alignment/2pass_genomeDir --readFilesIn $fastqfile --readFilesCommand 'gunzip -c' --outSAMtype BAM Unsorted --outFileNamePrefix $outfolder/alignment/2pass_alignment/$filename >> $logfile 2>&1
	fi
	
	printf "mkdir -p $outfolder/alignment/2pass_alignment;$star --genomeDir $outfolder/alignment/2pass_genomeDir --readFilesIn $fastqfile --readFilesCommand 'gunzip -c' --outSAMtype BAM Unsorted --outFileNamePrefix $outfolder/alignment/2pass_alignment/$filename >> $logfile 2>&1\n" | tee -a $runfile

	#sam file will be: $outfolder/alignment/2pass_alignment/${filename/fastq.gz/Aligned.out.sam}
else
	printf "#No Fastq file is defined. Skip STAR alignment.\n" | tee -a $logfile $runfile
fi


#need to add this 
#samtools view -Sb DKO1_Normal.Aligned.out.sam > DKO1_Normal.Aligned.out.bam

######
#Step 2
######


#GATK3 Start
#Add read groups, sort, mark duplicates, and create index
printf "\n\n#gatk sort, mark duplicates, and create index\n"  | tee -a  $logfile $runfile

if [ "$dryrun" == "F" ];then
	mkdir -p $outfolder/gatk3
fi

printf "mkdir -p $outfolder/gatk3\n" | tee -a $runfile

if [ -f "$fastqfile" ];then
	if [ "$dryrun" == "F" ];then
		eval $java -jar $gatk3folder/picard.jar AddOrReplaceReadGroups I=$outfolder/alignment/2pass_alignment/${filename}Aligned.out.bam O=$outfolder/gatk3/${filename}Aligned.out_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample USE_JDK_DEFLATER=true USE_JDK_INFLATER=true >> $logfile 2>&1 
	fi
	printf "eval $java -jar $gatk3folder/picard.jar AddOrReplaceReadGroups I=$outfolder/alignment/2pass_alignment/${filename}Aligned.out.bam O=$outfolder/gatk3/${filename}Aligned.out_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample USE_JDK_DEFLATER=true USE_JDK_INFLATER=true >> $logfile 2>&1\n" | tee -a $runfile
	
else
	#existing bam
	if [ "$dryrun" == "F" ];then
		eval $java -jar $gatk3folder/picard.jar AddOrReplaceReadGroups I=$bamfile O=$outfolder/gatk3/${filename}Aligned.out_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample USE_JDK_DEFLATER=true USE_JDK_INFLATER=true >> $logfile 2>&1 
	fi
	printf "eval $java -jar $gatk3folder/picard.jar AddOrReplaceReadGroups I=$bamfile O=$outfolder/gatk3/${filename}Aligned.out_added_sorted.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample USE_JDK_DEFLATER=true USE_JDK_INFLATER=true >> $logfile 2>&1\n" | tee -a $runfile
	
fi


if [ "$dedup" == "T" ]
then
	printf "\n\n#%cd T, run gatk MarkDuplicates\n" "-"  | tee -a  $logfile $runfile
	
	if [ "$dryrun" == "F" ];then
		eval $java -jar $gatk3folder/picard.jar MarkDuplicates I=$outfolder/gatk3/${filename}Aligned.out_added_sorted.bam O=$outfolder/gatk3/${filename}Aligned.out_dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics USE_JDK_DEFLATER=true USE_JDK_INFLATER=true >> $logfile 2>&1
	fi
	
	printf "eval $java -jar $gatk3folder/picard.jar MarkDuplicates I=$outfolder/gatk3/${filename}Aligned.out_added_sorted.bam O=$outfolder/gatk3/${filename}Aligned.out_dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics USE_JDK_DEFLATER=true USE_JDK_INFLATER=true >> $logfile 2>&1\n" | tee -a $runfile
	
else
	printf "\n\n#%cd F, skip gatk MarkDuplicates\n" "-"  | tee -a  $logfile $runfile
	
	if [ "$dryrun" == "F" ];then
		cp $outfolder/gatk3/${filename}Aligned.out_added_sorted.bam $outfolder/gatk3/${filename}Aligned.out_dedupped.bam
		$samtools index $outfolder/gatk3/${filename}Aligned.out_dedupped.bam $outfolder/gatk3/${filename}Aligned.out_dedupped.bai
	fi
	
	printf "cp $outfolder/gatk3/${filename}Aligned.out_added_sorted.bam $outfolder/gatk3/${filename}Aligned.out_dedupped.bam;$samtools index $outfolder/gatk3/${filename}Aligned.out_dedupped.bam $outfolder/gatk3/${filename}Aligned.out_dedupped.bai\n" | tee -a  $runfile
	
fi

######
#Step 3
######

printf "\n\n#gatk split and trim\n"  | tee -a  $logfile $runfile

#Split'N'Trim and reassign mapping qualities
if [ "$dryrun" == "F" ];then
	eval $java -jar $gatk3folder/GenomeAnalysisTK.jar -T SplitNCigarReads -R $genomefasta -I $outfolder/gatk3/${filename}Aligned.out_dedupped.bam -o $outfolder/gatk3/${filename}Aligned.out_split.bam -rf ReassignOneMappingQuality -U ALLOW_N_CIGAR_READS -RMQF 255 -RMQT 60 --use_jdk_deflater --use_jdk_inflater >> $logfile 2>&1
fi

printf "eval $java -jar $gatk3folder/GenomeAnalysisTK.jar -T SplitNCigarReads -R $genomefasta -I $outfolder/gatk3/${filename}Aligned.out_dedupped.bam -o $outfolder/gatk3/${filename}Aligned.out_split.bam -rf ReassignOneMappingQuality -U ALLOW_N_CIGAR_READS -RMQF 255 -RMQT 60 --use_jdk_deflater --use_jdk_inflater >> $logfile 2>&1\n" | tee -a  $runfile


#######
#Step 4
#######

#Indel Realignment (optional)

#######
#Step 5
#######

#Base Recalibration


#######
#Step 6
#######

#Picard ReplaceSamHeader
printf "\n\n#gatk haplotype calling\n"  | tee -a  $logfile $runfile

if [ "$dryrun" == "F" ];then
	eval $java -jar $gatk3folder/GenomeAnalysisTK.jar -T HaplotypeCaller -R $genomefasta -I $outfolder/gatk3/${filename}Aligned.out_split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o $outfolder/gatk3/${filename}output.vcf --use_jdk_deflater --use_jdk_inflater >> $logfile 2>&1
fi

printf "eval $java -jar $gatk3folder/GenomeAnalysisTK.jar -T HaplotypeCaller -R $genomefasta -I $outfolder/gatk3/${filename}Aligned.out_split.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -o $outfolder/gatk3/${filename}output.vcf --use_jdk_deflater --use_jdk_inflater >> $logfile 2>&1\n" | tee -a  $runfile


#######
#Step 7
#######

#Variant filtering
printf "\n\n#gatk variant filtering\n"  | tee -a  $logfile $runfile

if [ "$dryrun" == "F" ];then
	eval $java -jar $gatk3folder/GenomeAnalysisTK.jar -T VariantFiltration -R $genomefasta -V $outfolder/gatk3/${filename}output.vcf -window 35 -cluster 3 -filterName FS -filter \"FS \> 30.0\" -filterName QD -filter \"QD \< 2.0\" -filterName DP -filter \"DP \< $dpfilter\" -o $outfolder/gatk3/${filename}output.filtered.vcf --use_jdk_deflater --use_jdk_inflater >> $logfile 2>&1
fi

printf "eval $java -jar $gatk3folder/GenomeAnalysisTK.jar -T VariantFiltration -R $genomefasta -V $outfolder/gatk3/${filename}output.vcf -window 35 -cluster 3 -filterName FS -filter \"FS \> 30.0\" -filterName QD -filter \"QD \< 2.0\" -filterName DP -filter \"DP \< $dpfilter\" -o $outfolder/gatk3/${filename}output.filtered.vcf --use_jdk_deflater --use_jdk_inflater >> $logfile 2>&1\n" | tee -a $runfile



#######
#Step 8
#######

printf "\n\n#tabix indexing\n"  | tee -a  $logfile $runfile #steps after step 7 are for snpsift/snpEff analysis

#zip vcf file

if [ "$dryrun" == "F" ];then
	$bgzip -f $outfolder/gatk3/${filename}output.filtered.vcf >> $logfile 2>&1

	#index vcf
	$tabix -f -p vcf $outfolder/gatk3/${filename}output.filtered.vcf.gz >> $logfile 2>&1
fi

printf "$bgzip -f $outfolder/gatk3/${filename}output.filtered.vcf >> $logfile 2>&1;$tabix -f -p vcf $outfolder/gatk3/${filename}output.filtered.vcf.gz >> $logfile 2>&1\n" | tee -a  $runfile

#steps after step 7 are for snpsift/snpEff analysis
printf "\n\n#snpsift annotating\n"  | tee -a  $logfile $runfile

#annotate each vcf, may not needed, because there will be an extra step for -cancer annotation

if [ "$dryrun" == "F" ];then
	mkdir $outfolder/snpanno

	#annotate ID only using snpsift #may need MAF info later
	eval $java -jar $snpsift annotate -noInfo -v $knownsnp $outfolder/gatk3/${filename}output.filtered.vcf.gz > $outfolder/snpanno/${filename}output.filtered.sift.annotated.vcf 2>> $logfile
fi

printf "mkdir $outfolder/snpanno;eval $java -jar $snpsift annotate -noInfo -v $knownsnp $outfolder/gatk3/${filename}output.filtered.vcf.gz > $outfolder/snpanno/${filename}output.filtered.sift.annotated.vcf 2>> $logfile\n" | tee -a  $runfile


#filter here by common SNP and DP
if [ "$species" == "human" ]
then
	if [ "$dryrun" == "F" ];then
		#human
		perl $rnaseq_var_filter -i $outfolder/snpanno/${filename}output.filtered.sift.annotated.vcf -c $gnomadcommonsnp -o $outfolder/snpanno/${filename}output.filtered-cleaned.sift.annotated.vcf 2>> $logfile
	fi
	
	printf "perl $rnaseq_var_filter -i $outfolder/snpanno/${filename}output.filtered.sift.annotated.vcf -c $gnomadcommonsnp -o $outfolder/snpanno/${filename}output.filtered-cleaned.sift.annotated.vcf 2>> $logfile" | tee -a  $runfile
	
else
	
	if [ "$dryrun" == "F" ];then
		#mouse
		perl $rnaseq_var_filter -i $outfolder/snpanno/${filename}output.filtered.sift.annotated.vcf -o $outfolder/snpanno/${filename}output.filtered-cleaned.sift.annotated.vcf 2>> $logfile
	fi
	
	printf "perl $rnaseq_var_filter -i $outfolder/snpanno/${filename}output.filtered.sift.annotated.vcf -o $outfolder/snpanno/${filename}output.filtered-cleaned.sift.annotated.vcf 2>> $logfile" | tee -a  $runfile
fi




#annotate using snpEff for effects

printf "\n\n#snpeff annotating\n"  | tee -a  $logfile $runfile


if [ "$dryrun" == "F" ];then
	#####
	#before filter
	#####
	
	eval $java -jar $snpeff -v $genomeversion $outfolder/snpanno/${filename}output.filtered.sift.annotated.vcf -stats $outfolder/snpanno/snpEff_summary.html > $outfolder/snpanno/${filename}output.filtered.snpeff.sift.annotated.vcf 2>> $logfile

	#rename
	mv $outfolder/snpanno/snpEff_summary.genes.txt $outfolder/snpanno/${filename}snpEff_summary.genes.txt
	mv $outfolder/snpanno/snpEff_summary.html $outfolder/snpanno/${filename}snpEff_summary.html

	#reformat snpeff ANN column
	perl $reformatsnpeffvcf -i $outfolder/snpanno/${filename}output.filtered.snpeff.sift.annotated.vcf -o $outfolder/snpanno/${filename}output.filtered.snpeff.sift.annotated.edited.txt


	####
	#after filter
	####
	eval $java -jar $snpeff -v $genomeversion $outfolder/snpanno/${filename}output.filtered-cleaned.sift.annotated.vcf -stats $outfolder/snpanno/snpEff_filtered-cleaned_summary.html > $outfolder/snpanno/${filename}output.filtered-cleaned.snpeff.sift.annotated.vcf 2>> $logfile

	#rename
	mv $outfolder/snpanno/snpEff_filtered-cleaned_summary.genes.txt $outfolder/snpanno/${filename}snpEff_filtered-cleaned_summary.genes.txt
	mv $outfolder/snpanno/snpEff_filtered-cleaned_summary.html $outfolder/snpanno/${filename}snpEff_filtered-cleaned_summary.html

	#reformat snpeff ANN column
	perl $reformatsnpeffvcf -i $outfolder/snpanno/${filename}output.filtered-cleaned.snpeff.sift.annotated.vcf -o $outfolder/snpanno/${filename}output.filtered-cleaned.snpeff.sift.annotated.edited.txt

	
fi

printf "#unfiltered SNPs\n";
printf "eval $java -jar $snpeff -v $genomeversion $outfolder/snpanno/${filename}output.filtered.sift.annotated.vcf -stats $outfolder/snpanno/snpEff_summary.html > $outfolder/snpanno/${filename}output.filtered.snpeff.sift.annotated.vcf 2>> $logfile;mv $outfolder/snpanno/snpEff_summary.genes.txt $outfolder/snpanno/${filename}snpEff_summary.genes.txt;mv $outfolder/snpanno/snpEff_summary.html $outfolder/snpanno/${filename}snpEff_summary.html;perl $reformatsnpeffvcf $outfolder/snpanno/${filename}output.filtered.snpeff.sift.annotated.vcf $outfolder/snpanno/${filename}output.filtered.snpeff.sift.annotated.edited.txt\n" | tee -a  $runfile

printf "#filtered SNPs\n";
printf "eval $java -jar $snpeff -v $genomeversion $outfolder/snpanno/${filename}output.filtered-cleaned.sift.annotated.vcf -stats $outfolder/snpanno/snpEff_filtered-cleaned_summary.html > $outfolder/snpanno/${filename}output.filtered-cleaned.snpeff.sift.annotated.vcf;	mv $outfolder/snpanno/snpEff_filtered-cleaned_summary.genes.txt $outfolder/snpanno/${filename}snpEff_filtered-cleaned_summary.genes.txt;mv $outfolder/snpanno/snpEff_filtered-cleaned_summary.html $outfolder/snpanno/${filename}snpEff_filtered-cleaned_summary.html;perl $reformatsnpeffvcf $outfolder/snpanno/${filename}output.filtered-cleaned.snpeff.sift.annotated.vcf $outfolder/snpanno/${filename}output.filtered-cleaned.snpeff.sift.annotated.edited.txt\n" | tee -a  $runfile


######
#Step 9
######

#clean up to free disk space
#if don't free up, it will take ~40G per sample, after that, about 4G per sample


if [ "$keepalignment" == "F" ] && [ -f "$fastqfile" ]
then
	if [ "$dryrun" == "F" ];then
		rm $outfolder/alignment/*/*.sam
		rm $outfolder/alignment/2pass_genomeDir/SA*
		rm $outfolder/alignment/2pass_genomeDir/Genome
		
		rm $outfolder/gatk3/${filename}Aligned.out_added_sorted.bam
		rm $outfolder/gatk3/${filename}Aligned.out_dedupped.bam
		
	fi
	
	printf "rm $outfolder/alignment/*/*.sam;rm $outfolder/alignment/2pass_genomeDir/SA*;rm $outfolder/alignment/2pass_genomeDir/Genome\n" | tee -a  $runfile
	
fi


printf "\n\ndone\n" | tee -a  $logfile 

date >> $logfile
