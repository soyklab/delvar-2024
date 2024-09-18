#!/bin/bash

### Use this file to call and filter variants


### USAGE:
# batch variantcall_region.sbatch REGION
        # REGION = region in bam file (e.g. Spim0.1ch02; Spim0.1ch02:10000-20000)


### Requirements

# load modules
module add samtools
module add bcftools
module add vcftools

### parameters
echo "reading parameters"; date

#########################
# name experiment for file naming
experiment=82_acc

# select reference genome
ref=SpimLA1589_v0.1_chromosomes.fa
#########################

# read parameters
region=$1

# set date and name of output file
mkdir log
date=(date +%Y%m%d)
exec 1> log/variants_run_"$region"_`date '+%Y%m%d_%H%M%S'`.out 2>&1
mkdir 82_acc
mkdir 82_acc/vcf
tmp_vcf=82_acc/vcf
mkdir vcf
vcf=vcf

# write list of paths to selected bam files
bamdir=82_acc/bwa_results
ls $bamdir/*.bam | grep -f sample.list | sort  > bam.list

# write sample list for new samplenames in vcf header; extract from  bam.list (& double-check order!!)
ls $bamdir/*.bam | grep -f sample.list | rev |  cut -c 19- | rev | cut -c 39- | sort  > newsample.names
newsamplenames=newsample.names

echo "done reading parameters"; date



### Calls variants

echo "Start running mpileup"; date

bcftools mpileup --region "$region" --threads 16 --no-BAQ --ignore-RG -d 1000000 -Q0 --annotate FORMAT/AD,FORMAT/DP --fasta-ref "$ref" --bam-list bam.list | bcftools call --threads 16 -mvO z -o "$tmp_vcf"/"$experiment"_"$region"_prelim.vcf.gz

echo "Done running mpileup"; date



### Rename bcfheader
# use simple sample names
# list new sample names in the "newsamplenames" file with one new sample name per line in same order as listed in the vcf file header)

echo "Start changing vcf file header and index"; date

cat "$newsamplenames"
bcftools reheader -s "$newsamplenames" -o "$tmp_vcf"/"$experiment"_"$region".vcf.gz  "$tmp_vcf"/"$experiment"_"$region"_prelim.vcf.gz
bcftools index -f  "$tmp_vcf"/"$experiment"_"$region".vcf.gz

echo "Done changing vcf file header and index finished"; date



### Filter variants by type and quality
# Biallelic variants
# Quality > 30
# mac = 2
# D 5-to-100

echo "filter variants"; date

vcftools --gzvcf "$tmp_vcf"/"$experiment"_"$region".vcf.gz --min-alleles 2 --max-alleles 2 --minQ 30 --minDP 5 --maxDP 100 --mac 2 --recode --recode-INFO-all --stdout > "$tmp_vcf"/"$experiment"_"$region"_filtered.vcf
bcftools index -f "$tmp_vcf"/"$experiment"_"$region"_filtered.vcf

echo "filter variants finished"; date

echo "all done"; date

