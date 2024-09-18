### prep bed file with peak coordinates

cat peak_table_fdr_0.01.bed | sort -uk1,1 -k2,2n > peak_table_fdr_0.01_unique_sorted.bed

echo "nrow peak_table_fdr_0.01.bed"
cat peak_table_fdr_0.01.bed | wc -l

echo "peak_table_fdr_0.01_unique_sorted.bed"
cat peak_table_fdr_0.01_unique_sorted.bed | wc -l


### extract read coverage

module load singularityce
export SINGULARITY_BINDPATH="/work,/users,/scratch"
singularity run /dcsrsoft/singularity/containers/megadepth-1.2.0.sif --help

# select sample:
sample=SlycSSP2
#sample=SpimSSP2

bwdir=../../bw
bed=peak_table_fdr_0.01_unique_sorted.bed

singularity run /dcsrsoft/singularity/containers/megadepth-1.2.0.sif "$bwdir"/"$sample"_1_binned_filtered_normalized.bw --annotation "$bed" > "$sample"_1_cov.bed
singularity run /dcsrsoft/singularity/containers/megadepth-1.2.0.sif "$bwdir"/"$sample"_2_binned_filtered_normalized.bw --annotation "$bed" > "$sample"_2_cov.bed




