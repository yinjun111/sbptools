#!/bin/sh
# *******************************************
# Script to perform DNA seq variant calling
# using a single sample with fastq files
# named 1.fastq.gz and 2.fastq.gz
# *******************************************



#######
#Usage
#######

version="0.2"

usage="

dnaseq-process-sentieon

version: $version

Usage: sh sentieon_quickstart_caller.sh -f1 f1.fastq.gz  -f2 f2.fastq.gz -o Outputfolder

Description: Use sentieon workflow to call variants and CNV...

Currently PE sequencing only

Parameters:

	-f      Fastq_R1.fastq.gz
	-r      Fastq_R2.fastq.gz
	-o      Output folder
	
"



if [ $# -eq 0 ]; then 
	printf "$usage";
	exit 1; 
fi




#f1
#f2
#workdir

#data_dir=/home/jyin/Data/Pipeline_Test/dnaseq-process-sentieon/example1
#fastq_1=$data_dir/UACC1093-2_R1.fastq.gz
#fastq_2=$data_dir/UACC1093-2_R2.fastq.gz #If using Illumina paired data


#receive options
while getopts "f:r:o:" opt; do
  case ${opt} in
    f )
		fastq_1=$(realpath $OPTARG)
      ;;
    r )
		fastq_2=$(realpath $OPTARG)
      ;;
    o ) 
		workdir=$(realpath $OPTARG)
      ;;
    \? ) 
		printf "ERROR: Unknown options.\n\n$usage"
      ;;
	: ) printf "ERROR: Unknown options.\n\n$usage"
      ;;

  esac
done



## double comments are commands changed

# Update with the fullpath location of your sample fastq
set -x
#data_dir="$( cd -P "$( dirname "$0" )" && pwd )"
#data_dir=/home/jyin/Data/Pipeline_Test/dnaseq-process-sentieon/example1
#fastq_1=$data_dir/UACC1093-2_R1.fastq.gz
#fastq_2=$data_dir/UACC1093-2_R2.fastq.gz #If using Illumina paired data

# Update with the location of the reference data files
#fasta=$data_dir/reference/ucsc.hg19_chr22.fasta
fasta=/data/jyin/Pipeline_Dev/dnaseq-process-sentieon/hg38/Homo_sapiens.GRCh38.dna.primary_assembly_ucsc.fa

##dbsnp=$data_dir/reference/dbsnp_135.hg19_chr22.vcf
dbsnp=/data/jyin/Databases/gnomad/fromGATK4/af-only-gnomad.hg38.common001.vcf.gz

##known_1000G_indels=$data_dir/reference/1000G_phase1.snps.high_confidence.hg19_chr22.sites.vcf
##known_Mills_indels=$data_dir/reference/Mills_and_1000G_gold_standard.indels.hg19_chr22.sites.vcf

# Set SENTIEON_LICENSE if it is not set in the environment
## export SENTIEON_LICENSE=10.1.1.1:8990
export SENTIEON_LICENSE=10.1.66.46:8990


# Update with the location of the Sentieon software package
##SENTIEON_INSTALL_DIR=/home/release/sentieon-genomics-201808.05
SENTIEON_INSTALL_DIR=/home/apps/sentieon/

#7 Update with the location of temporary fast storage and uncomment
SENTIEON_TMPDIR=/tmp

# It is important to assign meaningful names in actual cases.
## It is particularly important to assign different read group names.
sample="sample_name"
group="read_group_name"
platform="ILLUMINA" 

# Other settings
nt=16 #number of threads to use in computation

# ******************************************
# 0. Setup
# ******************************************
#workdir=$data_dir/result4

mkdir -p $workdir
logfile=$workdir/run.log
exec >$logfile 2>&1
cd $workdir


echo "sh " $0 $@ > $logfile
date >> $logfile

cat "Fastq1:",$fastq_1
cat "Fastq2:",$fastq_1
cat "Output:",$workdir

#Sentieon proprietary compression
bam_option="--bam_compression 1"

# ******************************************
# 1. Mapping reads with BWA-MEM, sorting
# ******************************************
#The results of this call are dependent on the number of threads used. To have number of threads independent results, add chunk size option -K 10000000 

# speed up memory allocation malloc in bwa
export LD_PRELOAD=$SENTIEON_INSTALL_DIR/lib/libjemalloc.so
export MALLOC_CONF=lg_dirty_mult:-1

( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -M -R "@RG\tID:$group\tSM:$sample\tPL:$platform" -t $nt -K 10000000 $fasta $fastq_1 $fastq_2 || echo -n 'error' ) | $SENTIEON_INSTALL_DIR/bin/sentieon util sort $bam_option -r $fasta -o sorted.bam -t $nt --sam2bam -i -

# ******************************************
# 2. Metrics
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i sorted.bam --algo MeanQualityByCycle mq_metrics.txt --algo QualDistribution qd_metrics.txt --algo GCBias --summary gc_summary.txt gc_metrics.txt --algo AlignmentStat --adapter_seq '' aln_metrics.txt --algo InsertSizeMetricAlgo is_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot GCBias -o gc-report.pdf gc_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualDistribution -o qd-report.pdf qd_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot MeanQualityByCycle -o mq-report.pdf mq_metrics.txt
$SENTIEON_INSTALL_DIR/bin/sentieon plot InsertSizeMetricAlgo -o is-report.pdf is_metrics.txt

# ******************************************
# 3. Remove Duplicate Reads
# To mark duplicate reads only without removing them, remove "--rmdup" in the second command
# ******************************************
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i sorted.bam --algo LocusCollector --fun score_info score.txt
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -i sorted.bam --algo Dedup --rmdup --score_info score.txt --metrics dedup_metrics.txt $bam_option deduped.bam 

# ******************************************
# 4. Indel realigner
# This step is optional for haplotyper-based
# caller like HC, but necessary for any
# pile-up based caller. If you want to use
# this step, you need to update the rest of
# the commands to use realigned.bam instead
# of deduped.bam
# ******************************************
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam --algo Realigner -k $known_Mills_indels -k $known_1000G_indels $bam_option realigned.bam

# ******************************************
# 5. Base recalibration
# ******************************************

# Perform recalibration
##$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam --algo QualCal -k $dbsnp recal_data.table

# Perform post-calibration check (optional)
##$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam -q recal_data.table --algo QualCal -k $dbsnp -k $known_Mills_indels -k $known_1000G_indels recal_data.table.post
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam -q recal_data.table --algo QualCal -k $dbsnp recal_data.table.post
$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt --algo QualCal --plot --before recal_data.table --after recal_data.table.post recal.csv   
$SENTIEON_INSTALL_DIR/bin/sentieon plot QualCal -o recal_plots.pdf recal.csv

# ******************************************
# 5b. ReadWriter to output recalibrated bam
# This stage is optional as variant callers
# can perform the recalibration on the fly
# using the before recalibration bam plus
# the recalibration table
# ******************************************
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam -q recal_data.table --algo ReadWriter recaled.bam


# ******************************************
# 6. HC Variant caller
# Note: Sentieon default setting matches versions before GATK 3.7.
# Starting GATK v3.7, the default settings have been updated multiple times. 
# Below shows commands to match GATK v3.7 - 4.1
# Please change according to your desired behavior.
# ******************************************

# Matching GATK 3.7, 3.8, 4.0
##$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam -q recal_data.table --algo Haplotyper -d $dbsnp --emit_conf=10 --call_conf=10 output-hc.vcf.gz


# Matching GATK 4.1
#$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam -q recal_data.table --algo Haplotyper -d $dbsnp --genotype_model multinomial --emit_conf 30 --call_conf 30 output-hc.vcf.gz

$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $fasta -t $nt -i deduped.bam -q recal_data.table --algo Haplotyper -d $dbsnp --genotype_model multinomial --emit_conf 30 --call_conf 30 output-hc.vcf.gz


# ******************************************
# 7. CNV caller
# ******************************************

$SENTIEON_INSTALL_DIR/bin/sentieon driver -t $nt -r $fasta -i deduped.bam --algo CNV output-cnv.vcf.gz
