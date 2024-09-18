#!/bin/bash

# Use this script to concatenate vcf files from separate chromosomes


# define experiment

experiment=82_acc

# set date and name of output file
mkdir log
date=(date +%Y%m%d)
exec 1> log/concat_vcf_`date '+%Y%m%d_%H%M%S'`.out 2>&1

# load modules
module load bcftools
module load vcftools

# directories
mkdir 82_acc/vcf
tmp_vcf=82_acc/vcf
mkdir vcf
vcf=vcf


echo "start concatenating vcf files" ; date

# concatenate files
bcftools concat -o "$tmp_vcf"/"$experiment"_Spim0.1_filtered.vcf \
"$tmp_vcf"/"$experiment"_Spim0.1ch00_filtered.vcf \
"$tmp_vcf"/"$experiment"_Spim0.1ch01_filtered.vcf \
"$tmp_vcf"/"$experiment"_Spim0.1ch02_filtered.vcf \
"$tmp_vcf"/"$experiment"_Spim0.1ch03_filtered.vcf \
"$tmp_vcf"/"$experiment"_Spim0.1ch04_filtered.vcf \
"$tmp_vcf"/"$experiment"_Spim0.1ch05_filtered.vcf \
"$tmp_vcf"/"$experiment"_Spim0.1ch06_filtered.vcf \
"$tmp_vcf"/"$experiment"_Spim0.1ch07_filtered.vcf \
"$tmp_vcf"/"$experiment"_Spim0.1ch08_filtered.vcf \
"$tmp_vcf"/"$experiment"_Spim0.1ch09_filtered.vcf \
"$tmp_vcf"/"$experiment"_Spim0.1ch10_filtered.vcf \
"$tmp_vcf"/"$experiment"_Spim0.1ch11_filtered.vcf \
"$tmp_vcf"/"$experiment"_Spim0.1ch12_filtered.vcf 

## index file
bcftools index "$vcf"/"$experiment"_Spim0.1_filtered.vcf.gz

echo "done" ; date


