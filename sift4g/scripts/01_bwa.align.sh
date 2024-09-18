#!/bin/bash

### Use this script to align reads to reference genome


# define list of samples
samples=acc.list # provide list of genome samples (accessions)

# set date and name of output file
mkdir log
date=(date +%Y%m%d)
exec 1> log/bwa.align_run_"$samples"_`date '+%Y%m%d_%H%M%S'`.out 2>&1

# destinations
fastqdir=/fastq # directory with fastq files
mkdir /82_acc
mkdir /82_acc/bwa_logs
mkdir /82_acc/bwa_results
bwa_logs=/82_acc/bwa_logs
bwa_results=/82_acc/bwa_results

# load modules
module add bwa #aligner
module add picard #mark duplicates
module add samtools #sorting and indexing

# load refs (make sure that bwa index has been created in same directory)
ref=SpimLA1589_v0.1_chromosomes.fa

###################################################
# Creating read alignments
##################################################

echo "### start creating read alignments ###"; date

for SAMPLE in $(cut -f1 "$samples");
do

# Aligning reads with BWA-MEM

	echo "align reads for sample ${SAMPLE}"; date
	bwa mem -M -t 16 \
		$ref \
		"$fastqdir"/${SAMPLE}_R1.fastq "$fastqdir"/${SAMPLE}_R2.fastq \
		2> "$bwa_logs"/${SAMPLE}_bwa.err \
		> "$bwa_results"/${SAMPLE}.sam
	echo "### aligning reads for sample ${SAMPLE} finished  ###"; date

# Sorting SAM by coordinates

	echo "sort sam for sample ${SAMPLE}"; date
	picard SortSam \
		INPUT="$bwa_results"/${SAMPLE}.sam \
		OUTPUT="$bwa_results"/${SAMPLE}_sorted.sam \
		SORT_ORDER=coordinate \
		VALIDATION_STRINGENCY=SILENT
	echo "### sorting sams for sample ${SAMPLE} finished  ###"; date

# Marking duplicates

        echo "mark duplicates for sample ${SAMPLE}"; date
		picard  MarkDuplicates \
		INPUT="$bwa_results"/${SAMPLE}_sorted.sam \
		OUTPUT="$bwa_results"/${SAMPLE}_sorted_marked.bam \
		METRICS_FILE="$bwa_results"/${SAMPLE}_metrics.txt \
		ASSUME_SORTED=true \
		VALIDATION_STRINGENCY=SILENT
	echo "### marking duplicates for sample ${SAMPLE} finished  ###"; date

# Creating index for BAM file

        echo "create bam index for sample ${SAMPLE}"; date
		samtools index "$bwa_results"/${SAMPLE}_sorted_marked.bam
	echo "### creating bam indexes for sample ${SAMPLE} finished  ###"; date

# Clear sam files
	
	rm "$bwa_results"/${SAMPLE}*.sam

done
echo "### finished creating read alignments ###"; date


