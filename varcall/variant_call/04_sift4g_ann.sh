#!/bin/bash

### Use this script to annotate deleterious variants in vcf file with sift-4g 

### USAGE:
# batch 04_sift4g_ann.sh REF SAMPLE
		#REF= reference genome (e.g. LA1589)
		#SAMPLE= variant call file (e.g. 82_acc_Spim0.1)
		
# define reference and variant file
prefix1=`echo $1|awk -F"/" '{print $(NF)}'`
prefix2=`echo $2|awk -F"/" '{print $(NF)}'`
ref=$1
sample=$2

# set date and name of output file
mkdir log
date=(date +%Y%m%d)
exec 1> log/sift4g_ann_"$ref"_"$sample"_`date '+%Y%m%d_%H%M%S'`.out 2>&1

# load modules
module load bcftools

# Make tmp and results directory
mkdir sift4g_results
results_dir=sift4g_results #Path to your results folder

# Path to sift4g annotator executable
sift4g_ann=/FILE-PATH-TO/SIFT4G_Annotator.jar

# Path to input vcf file
vcf=vcf

# Path to SIFT4G database directory
sift_db_dir=/FILE-PATH-TO/sift4g_lib/

# run sift4g annotator
java -jar "$sift4g_ann" -c -i "$vcf"/"$sample"_filtered.vcf -d "$sift_db_dir"/"$ref"/"$ref" -r "$results_dir" -t

# write table format for R analysis
bcftools query -f '[%CHROM\t%POS\t%REF\t%ALT\t%SAMPLE\t%GT\t%INFO/SIFTINFO\n]' "$results_dir"/"$sample"_filtered_SIFTpredictions.vcf | grep CDS > "$results_dir"/"$sample"_filtered_SIFTpredictions.list





