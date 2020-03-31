#!/bin/sh

#######
#Author: Jun Yin <jyin@sbpdiscovery.org>
#SBP Bioinformatics Core
#######

version="1.41"

#version 1.1 add option to control params
#version 1.2 add reformating snpeff and rename snpEff_summary.genes
#version 1.3 remove huge alignment and index files to free up space
#v1.3a, add -f to tabix
#v1.4, compatibility in Firefly
#v1.41 fix samtools bug

#######
#Usage
#sh gatk3_rnaseq_variant_macos_dm_v1.sh yourfastq.fastq.gz outputfolder
#######

usage="

rnaseq-var

version: $version

Usage: sbptools rnaseq-var -i yourfastq.fastq.gz -o outputfolder -s species

Description: RNA-Seq variant calling pipeline built using GATK3. This pipeline is customized to call and annotate variants from RNA-Seq data for human/mouse using B38 annotation. You only need to provde fastq file and an output folder.

Currently SE sequencing only

Parameters:

	-i      Input fastq file
	-o      Output folder
	-s      Species, human or mouse

	Optional
	-d      Whethter to dedup [T]
	-f      DP filter [5.0]
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
keepalignment=F

#receive options
while getopts ":i:o:s:d:f:" opt; do
  case ${opt} in
    i )
		fastqfile=$(realpath $OPTARG)
		fastqfilename=$(basename $fastqfile)
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
    k ) 
		keepalignment=$OPTARG
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


logfile=$outfolder/${fastqfilename/fastq.gz/run.log}


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


#java=/usr/bin/java						  
java=java
star=/apps/STAR-master/bin/Linux_x86_64/STAR
gatk3folder=/apps/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef
bgzip=/apps/htslib-1.9/bgzip
tabix=/apps/htslib-1.9/tabix
snpeff=/apps/snpEff/snpEff.jar
snpsift=/apps/snpEff/SnpSift.jar
reformatsnpeffvcf=/apps/sbptools/rnaseq-var/reformat_snpeff_vcf.pl
samtools=/apps/bin/samtools

#######
#Step 0
#######
#samtools faidx $genomefasta

#java -jar $gatk3folder/picard.jar CreateSequenceDictionary R= $genomefasta O= ${genomefasta/.fa/.dict}


#######
#Step 1, 2-pass algnment
#######


#1st-pass
printf "1-pass alignment\n" | tee -a  $logfile

mkdir -p $outfolder/alignment/1pass_alignment

$star --genomeDir $genomedir --readFilesIn $fastqfile --readFilesCommand 'gunzip -c' --outFileNamePrefix $outfolder/alignment/1pass_alignment/${fastqfilename/fastq.gz/} >> $logfile 2>&1

#2-pass #need to use a different alignment output folder ...
printf "2-pass alignment\n" | tee -a  $logfile

mkdir -p $outfolder/alignment/2pass_genomeDir

$star --runMode genomeGenerate --genomeDir $outfolder/alignment/2pass_genomeDir --outFileNamePrefix $outfolder/alignment/2pass_genomeDir/${fastqfilename/fastq.gz/} --genomeFastaFiles $genomefasta --sjdbFileChrStartEnd $outfolder/alignment/1pass_alignment/${fastqfilename/fastq.gz/SJ.out.tab} --sjdbOverhang 75 --runThreadN 4 >> $logfile 2>&1

#2-pass alignment
mkdir -p $outfolder/alignment/2pass_alignment
$star --genomeDir $outfolder/alignment/2pass_genomeDir --readFilesIn $fastqfile --readFilesCommand 'gunzip -c' --outFileNamePrefix $outfolder/alignment/2pass_alignment/${fastqfilename/fastq.gz/} >> $logfile 2>&1

#sam file will be: $outfolder/alignment/2pass_alignment/${fastqfilename/fastq.gz/Aligned.out.sam}

#need to add this 
#samtools view -Sb DKO1_Normal.Aligned.out.sam > DKO1_Normal.Aligned.out.bam

######
#Step 2
######


#GATK3 Start
#Add read groups, sort, mark duplicates, and create index
printf "gatk sort, mark duplicates, and create index\n"  | tee -a  $logfile

mkdir -p $outfolder/gatk3

eval $java -jar $gatk3folder/picard.jar AddOrReplaceReadGroups I=$outfolder/alignment/2pass_alignment/${fastqfilename/fastq.gz/Aligned.out.sam} O=$outfolder/gatk3/${fastqfilename/fastq.gz/Aligned.out_added_sorted.bam} SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample USE_JDK_DEFLATER=true USE_JDK_INFLATER=true >> $logfile 2>&1 


if [ "$dedup" == "T" ]
then
	printf "%cd T, run gatk MarkDuplicates\n" "-"  | tee -a  $logfile
	eval $java -jar $gatk3folder/picard.jar MarkDuplicates I=$outfolder/gatk3/${fastqfilename/fastq.gz/Aligned.out_added_sorted.bam} O=$outfolder/gatk3/${fastqfilename/fastq.gz/Aligned.out_dedupped.bam}  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=output.metrics USE_JDK_DEFLATER=true USE_JDK_INFLATER=true >> $logfile 2>&1
else
	printf "%cd F, skip gatk MarkDuplicates\n" "-"  | tee -a  $logfile
	cp $outfolder/gatk3/${fastqfilename/fastq.gz/Aligned.out_added_sorted.bam} $outfolder/gatk3/${fastqfilename/fastq.gz/Aligned.out_dedupped.bam}
	$samtools index $outfolder/gatk3/${fastqfilename/fastq.gz/Aligned.out_dedupped.bam} $outfolder/gatk3/${fastqfilename/fastq.gz/Aligned.out_dedupped.bai}
fi

######
#Step 3
######

printf "gatk split and trim\n"  | tee -a  $logfile

#Split'N'Trim and reassign mapping qualities
eval $java -jar $gatk3folder/GenomeAnalysisTK.jar -T SplitNCigarReads -R $genomefasta -I $outfolder/gatk3/${fastqfilename/fastq.gz/Aligned.out_dedupped.bam} -o $outfolder/gatk3/${fastqfilename/fastq.gz/Aligned.out_split.bam} -rf ReassignOneMappingQuality -U ALLOW_N_CIGAR_READS -RMQF 255 -RMQT 60 --use_jdk_deflater --use_jdk_inflater >> $logfile 2>&1


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
printf "gatk haplotype calling\n"  | tee -a  $logfile

eval $java -jar $gatk3folder/GenomeAnalysisTK.jar -T HaplotypeCaller -R $genomefasta -I $outfolder/gatk3/${fastqfilename/fastq.gz/Aligned.out_split.bam} -dontUseSoftClippedBases -stand_call_conf 20.0 -o $outfolder/gatk3/${fastqfilename/fastq.gz/output.vcf} --use_jdk_deflater --use_jdk_inflater >> $logfile 2>&1

#######
#Step 7
#######

#Variant filtering
eval $java -jar $gatk3folder/GenomeAnalysisTK.jar -T VariantFiltration -R $genomefasta -V $outfolder/gatk3/${fastqfilename/fastq.gz/output.vcf} -window 35 -cluster 3 -filterName FS -filter \"FS \> 30.0\" -filterName QD -filter \"QD \< 2.0\" -filterName DP -filter \"DP \< $dpfilter\" -o $outfolder/gatk3/${fastqfilename/fastq.gz/output.filtered.vcf} --use_jdk_deflater --use_jdk_inflater >> $logfile 2>&1



#######
#Step 8
#######

printf "tabix indexing\n"  | tee -a  $logfile #steps after step 7 are for snpsift/snpEff analysis

#zip vcf file
$bgzip -f $outfolder/gatk3/${fastqfilename/fastq.gz/output.filtered.vcf} >> $logfile 2>&1

#index vcf
$tabix -f -p vcf $outfolder/gatk3/${fastqfilename/fastq.gz/output.filtered.vcf.gz} >> $logfile 2>&1

#steps after step 7 are for snpsift/snpEff analysis
printf "snpsift/snpeff annotating\n"  | tee -a  $logfile 

#annotate each vcf, may not needed, because there will be an extra step for -cancer annotation
mkdir $outfolder/snpanno

#annotate ID only using snpsift #may need MAF info later
eval $java -jar $snpsift annotate -noInfo -v $knownsnp $outfolder/gatk3/${fastqfilename/fastq.gz/output.filtered.vcf.gz} > $outfolder/snpanno/${fastqfilename/fastq.gz/output.filtered.sift.annotated.vcf} 2>> $logfile

#annotate using snpEff for effects
eval $java -jar $snpeff -v $genomeversion $outfolder/snpanno/${fastqfilename/fastq.gz/output.filtered.sift.annotated.vcf} -stats $outfolder/snpanno/snpEff_summary.html > $outfolder/snpanno/${fastqfilename/fastq.gz/output.filtered.snpeff.sift.annotated.vcf} 2>> $logfile

#rename
mv $outfolder/snpanno/snpEff_summary.genes.txt $outfolder/snpanno/${fastqfilename/fastq.gz/snpEff_summary.genes.txt}
mv $outfolder/snpanno/snpEff_summary.html $outfolder/snpanno/${fastqfilename/fastq.gz/snpEff_summary.html}

#reformat snpeff ANN column
perl $reformatsnpeffvcf $outfolder/snpanno/${fastqfilename/fastq.gz/output.filtered.snpeff.sift.annotated.vcf} $outfolder/snpanno/${fastqfilename/fastq.gz/output.filtered.snpeff.sift.annotated.edited.txt}


######
#Step 9
######

#clean up to free disk space
#if don't free up, it will take ~40G per sample, after that, about 4G per sample

if [ "$keepalignment" == "F" ]
then
	rm $outfolder/alignment/*/*.sam
	rm $outfolder/alignment/2pass_genomeDir/SA*
	rm $outfolder/alignment/2pass_genomeDir/Genome
fi



printf "done\n" | tee -a  $logfile 

date >> $logfile
