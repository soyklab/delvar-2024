#-------------------------------------------------------------------------------------------
#
# PURPOSE: Differential binding analysis of ChIP-seq data using the csaw tool
#
# VERSION: 1
#
# AUTHOR : Giovanna Ambrosini
#          giovanna.ambrosini@epfl.ch 
#
#
# INPUT :  List of bam files and metadata file describing the samples and experimental design
# OUTPUT:  Differentially bound regions in different formats (BED, narrowPeak, tsv and 
#          annotated tsv) and data traces in bigWig format for visualisation
#
#
# REMARK: This script is based on the csaw User's Guide
#         https://bioconductor.org/books/release/csawBook/
#
#---------------------------------------------------------------------------------------------
#
library(csaw)
library(edgeR)
library(Rsamtools)
library(writexl)
library(naturalsort)
library(ggplot2)
library(BiocFileCache)
library(rtracklayer)
library(BiocParallel)
library(parallel)
library(Gviz)
library(clusterProfiler)
library(GenomicFeatures)  # makeTxDbFromGFF function

### Bioconductor version 3.13 (BiocManager 1.30.18), R 4.1.1 (2021-08-10)

#
### csaw pipeline
#
# It's a fairly modular process:

# 1. Counting reads into windows
# 2. Computing normalization factors
# 3. Filtering out low-abundance windows
# 4. Testing for differential binding
# 5. Aggregating windows into clusters
# 6. Visualization, annotation and other stuff

# The definition of these variables depends on the structure 
# of the working directory (can be changed)
setwd("/Users/giovanna/tomato_bZIP_sebastien_soyk/analysis")
inputdir <- "/Users/giovanna/tomato_bZIP_sebastien_soyk/data/bams"
outputdir <- "."
metadatafile <- "/Users/giovanna/tomato_bZIP_sebastien_soyk/data/metadata.txt"
getwd()

### Set bam files path                                                            ####
metadata <-  read.table(metadatafile,sep="\t",stringsAsFactors=FALSE,header=TRUE)
metadata
metadata$filename=normalizePath(paste0(inputdir,"/",metadata$sample,".bam"))
metadata

if(!all(file.exists(metadata$filename))) {
    cat("file not found:\n")
    cat(paste0(" ",metadata[!file.exists(metadata$filename),"filename"],"\n",collapse=""))
    cat("\n")
    stop("Aborting")
}
rownames(metadata)=metadata$filename
metadata

#  Pre-processing checks                                                          ####
#  Checking some mapping statistics for the SSP/SSP2 dataset                      ####
diagnostics <- list()
for (b in seq_along(metadata$filename)) {
    bam <- metadata$filename[b]
    total <- countBam(bam)$records
    mapped <- countBam(bam, param=ScanBamParam(
    flag=scanBamFlag(isUnmapped=FALSE)))$records
    diagnostics[[b]] <- c(Total=total, Mapped=mapped, Marked=0)
}

diag.stats <- data.frame(do.call(rbind, diagnostics))
diag.stats

rownames(diag.stats) <- metadata$sample
diag.stats

write.table(diag.stats, file = "bam_diagnostics.txt", quote = F, col.names = T, row.names = T, sep = "\t")

# Definition of readParam to standardise parameter settings for this analysis     ####
# 
# Obtaining the ENCODE blacklist                                                  ####
# Tomato genome S100_v2.1.0                                                       ####
standard.chr <- paste0("chr", c(1:12))
f <- "/Users/giovanna/tomato_bZIP_sebastien_soyk/results/genome_files/S100_v2.1.0_blacklist.bed"
# Import blacklist                                                                ####
blacklist=GRanges()
cat("loading regions to discard from",f,"\n")
tmp=import.bed(f)
blacklist=reduce(c(blacklist, tmp))
length(blacklist)

#
# All non-input files together
cat("evaluating cross correlation (all files)\n")
bams <- metadata$filename[-c(1,2)]
control <-  metadata$filename[c(1,2)]
bams

# Setting up extraction parameters                                               ####
cores <- as.numeric(min(8,detectCores()))
multicoreParam <- MulticoreParam(workers = cores)
# Paired-end
param <- readParam(minq=20, discard=blacklist, restrict=standard.chr, pe="first")
# Ignore duplicate reads
param=reform(param, dedup=TRUE)

# Cross correlation for Computing the average fragment length                    ####
x <- correlateReads(metadata$filename, param=param, BPPARAM=multicoreParam)      ### All bams
saveRDS(file='correlateReads.rds', x)
### Load read correlation file                                                   #### 
x <- readRDS(file='correlateReads.rds')
frag.len <- maximizeCcf(x)
frag.len

pdf(file = "Frag_len.pdf", width = 8, height = 7, useDingbats = F)
par(mar=c(5,9,2,3))
plot(1:length(x)-1, x, xlab="Delay (bp)", ylab="CCF", type="l", cex.lab=3, cex.axis=2,
     mgp = c(3.5, 1, 0))
abline(v=frag.len, col="red")
text(x=frag.len, y=min(x), paste(frag.len, "bp"), pos=4, col="red", cex=2)
dev.off()

### Compute the length of the sequenced fragment for each read pair in paired-end tag (PE) data.  ###
getPESizes(bams[1], param=readParam(pe="both"))

### Counting reads into windows                                                   ####
#
# Reads are then counted into sliding windows using csaw (Lun and Smyth 2016). 
# For TF data analyses, smaller windows are necessary to capture sharp binding sites. 
# A large window size will be suboptimal as the count for a particular site will be 
# “contaminated” by non-specific background in the neighbouring regions. 
# In this case, a window size of 10 bp is used.
#
# The windowCounts function returns a RangedSummarizedExperiment object where the matrix of counts 
# is stored as the first assay. Each row corresponds to a genomic window while each column corresponds to a library.
# The coordinates of each window are stored in the rowRanges. The total number of reads in each library 
# (also referred to as the library size) is stored as totals in the colData.
#
# In this experiment, BAM files haven't been marked for duplication (so dedup=TRUE is not used).
#
# We use all BAMS as later we want to compare SSP (IP) vs control. 
#
# The csaw pipeline can also be applied to search for DB between ChIP libraries and control libraries. 
# The ChIP and control libraries can be treated as separate groups, in which most DB events are expected
# to be enriched in the ChIP samples. 
#
# Define/Read param file                                                         ####
param <- readParam(minq=20, discard=blacklist, restrict=standard.chr, pe="both", max.frag = 1000)
#param <- readRDS('readParams.rds')
param

window.width <- 10
win.data <- windowCounts(metadata$filename, param=param, width=window.width, ext=frag.len)
win.data
dim(win.data)

# From previous analysis                                                         ####
#
# Preview the counts:
head(assay(win.data))
# Preview the genomic coordinates: 
head(rowRanges(win.data))
# Preview the totals
win.data$totals

# Save win.data                                                                  ####
saveRDS(file='win_data.rds', win.data)

# Filtering of low-abundance windows                                             ####
#
# Filtering by global abundance                                                  ####
#
# We remove low-abundance windows by computing the coverage in each window relative
# to a global  estimate of background enrichment (Section 4.4 of csaw Book).
# The majority of windows in background regions are filtered out upon applying a 
# modest fold-change threshold. 
# This leaves a small set of relevant windows for further analysis.
# logFC=3
# logFC=4
# logFC=5
#
# Binning is necessary here to increase the size of the counts when examining low-density background regions.
# This ensures that precision is maintained when estimating the background abundance.
#
# The median of the average abundances across all bins is computed and used as a global estimate of the background coverage. 
# This global background is then compared to the window-based abundances. This determines whether a window 
# is driven by background enrichment, and thus, unlikely to be interesting.
#
# However, some care is required as the sizes of the regions used for read counting are different between bins and windows.
# The average abundance of each bin must be scaled down to be comparable to those of the windows.
#
# The filterWindowsGlobal() function returns the increase in the abundance of each window over the 
# global background. Windows are filtered by setting some minimum threshold on this increase. 
#
# The aim is to eliminate the majority of uninteresting windows prior to further analysis.
# Here, a fold change of 3 or 4 is necessary for a window to be considered as containing a binding site. 
#
# Note that the 10 kbp bins are used here for filtering.
#
#
# Estimate background.bins                                                        ####
background.bin.width=10000
bins <- windowCounts(metadata$filename, bin=TRUE, width=background.bin.width, param=param)

saveRDS(file='bkg_bins.rds', bins)
bins <- readRDS('bkg_bins.rds')

filter.stat <- filterWindowsGlobal(win.data, bins)

# Filtering options                                                               ####
#filtering_fc <- 3
filtering_fc <- 4
#filtering_fc <- 5
filtering.min.fc <- filtering_fc
filtering.min.fc
keep <- filter.stat$filter > log2(filtering.min.fc)
summary(keep)

tt=table(factor(keep,levels=c(TRUE,FALSE)))
cat(" Nb windows  kept:",tt["TRUE"]," (",signif(100*tt["TRUE"]/sum(tt),2),"%), filtered out:",tt["FALSE"],"(",signif(100*tt["FALSE"]/sum(tt),2),"%)\n")
#Nb windows  kept: 391063  ( 4.2 %), filtered out: 8875933 ( 96 %)

### Plot
#
# Distribution of the log-increase in coverage over the global background for each window. 
# The red line denotes the chosen threshold for filtering.
#
# We can visualize the effect of filtering to confirm that the bulk of windows - 
# presumably in background regions - are indeed discarded upon filtering. 
#
# One might hope to see a bimodal distribution due to windows containing genuine binding sites,
# but this is usually not visible due to the dominance of background regions in the genome.
#
# For TF data, a large cut-off works well as narrow binding sites will have high read densities
# and are unlikely to be lost during filtering. 
# Smaller minimum fold changes are recommended for diffuse marks where the difference from background is less obvious.
#
getwd()
pdf(file = "Filtering_stats_fc4.pdf", 10,10)
par(mar=c(5,7,2,3))
hist(filter.stat$filter, main="", xlab="Log-fold change from global background", 
     breaks=100, col="grey80", xlim=c(0, 5), cex.lab=2.5, cex.axis=2,
     mgp = c(3.5, 1, 0))
abline(v=log2(3), col="red", lwd=2)
abline(v=log2(4), col="red", lwd=2)
dev.off()

filtered.data <- win.data[keep,]
saveRDS(file='filtered_data.rds', filtered.data)

# Read Filtered data                                                                ####
filtered.data <- readRDS('filtered_data.rds')
filtered.data

# Normalizing for technical biases                                                  ####
#
# This includes composition biases, efficiency biases and trended biases.
#
# Normalization for composition biases                                              ####
#
# Composition biases are formed when there are differences in the composition of sequences across libraries.
# Highly enriched regions consume more sequencing resources and thereby suppress the representation of other regions. 
# Differences in the magnitude of suppression between libraries can lead to spurious DB calls. 
# Scaling by library size fails to correct for this as composition biases can still occur in libraries of the same size.
#
# We expect unbalanced DB in this dataset as SPP2 function should be compromised in the domesticated Tomato 
#
# Using the TMM method on binned counts
#
# This uses the trimmed mean of M-values (TMM) method (Robinson and Oshlack 2010) to correct 
# for any systematic fold change in the coverage of the bins.
#
# To remove this bias, we assign reads to large genomic bins and assume that most bins represent 
# non-DB background regions (Lun and Smyth 2014).
# Any systematic differences in the coverage of those bins is attributed to composition bias and
# is normalized out. Specifically, the trimmed mean of M-values (TMM) method (Robinson and Oshlack 2010) 
# is applied to compute normalization factors from the bin counts.
#
# These factors are stored in win.data so that they will be applied during the DB analysis 
# with the window counts.
#
filtered.data <- normFactors(bins, se.out=filtered.data)
(normfacs <- filtered.data$norm.factors)

saveRDS(file='filtered_data.rds', filtered.data)
### Load normalized filtered data                                                 ####
filtered.data <- readRDS('filtered_data.rds')
filtered.data

# Visualization                                                                   ####
#
# We visualize the effect of normalization with mean-difference plots between pairs of samples.
# The red line represents the log-ratio of the normalization factors between samples.
# The dense cloud in each plot represents the majority of bins in the genome. 
# These are assumed to mostly contain background regions.
# A non-zero log-fold change for these bins indicates that composition bias is present between samples.
# The red line represents the log-ratio of normalization factors and passes through the centre 
# of the cloud in each plot, indicating that the bias has been successfully identified and removed.
#
bin.ab <- scaledAverage(bins)
adjc <- calculateCPM(bins, use.norm.factors=FALSE)

par(cex.lab=1.5, mfrow=c(3,3))
smoothScatter(bin.ab, adjc[,1]-adjc[,3], ylim=c(-6, 6),
              xlab="Average abundance", ylab="Log-ratio (1 vs 3)")
abline(h=log2(normfacs[1]/normfacs[3]), col="red")

smoothScatter(bin.ab, adjc[,2]-adjc[,3], ylim=c(-6, 6),
              xlab="Average abundance", ylab="Log-ratio (2 vs 3)")
abline(h=log2(normfacs[2]/normfacs[3]), col="red")

smoothScatter(bin.ab, adjc[,3]-adjc[,3], ylim=c(-6, 6),
              xlab="Average abundance", ylab="Log-ratio (3 vs 3)")
abline(h=log2(normfacs[3]/normfacs[3]), col="red")

smoothScatter(bin.ab, adjc[,4]-adjc[,3], ylim=c(-6, 6),
              xlab="Average abundance", ylab="Log-ratio (4 vs 3)")
abline(h=log2(normfacs[3]/normfacs[3]), col="red")

smoothScatter(bin.ab, adjc[,5]-adjc[,3], ylim=c(-6, 6),
              xlab="Average abundance", ylab="Log-ratio (5 vs 3)")
abline(h=log2(normfacs[3]/normfacs[3]), col="red")

smoothScatter(bin.ab, adjc[,6]-adjc[,3], ylim=c(-6, 6),
              xlab="Average abundance", ylab="Log-ratio (6 vs 3)")
abline(h=log2(normfacs[3]/normfacs[3]), col="red")

smoothScatter(bin.ab, adjc[,7]-adjc[,3], ylim=c(-6, 6),
              xlab="Average abundance", ylab="Log-ratio (7 vs 3)")
abline(h=log2(normfacs[3]/normfacs[3]), col="red")

smoothScatter(bin.ab, adjc[,8]-adjc[,3], ylim=c(-6, 6),
              xlab="Average abundance", ylab="Log-ratio (8 vs 3)")
abline(h=log2(normfacs[3]/normfacs[3]), col="red")

### Statistical modelling  of biological variability                    ####
#
# We model counts for each window using edgeR (McCarthy, Chen, and Smyth 2012; Robinson, McCarthy, and Smyth 20103). 
# First, we convert our RangedSummarizedExperiment object into a DGEList
#
y <- asDGEList(filtered.data)
summary(y)

# Design matrix for our experimental design                                      ####
#
# We have a simple one-way layout with (three + Input) groups of two replicates
#
antibody <- metadata$antibody
antibody

antibody <- factor(antibody)
antibody
design <- model.matrix(~0+antibody)
colnames(design) <- levels(antibody)
design
#
# Estimating the dispersions                                                      ####
# We use the negative binomial (NB) and quasi-likelihood (QL) dispersions 
# for each window. 
#
# NB
# This accounts for low, discrete counts that exhibit over dispersion between biological replicates.
# Specifically, variability between replicates is modelled using the NB dispersion parameter, 
# as estimated through the estimateDisp function.
#
y <- estimateDisp(y, design)
summary(y$trended.dispersion)
# Plot the NB dispersions 
# Only the trended NB dispersions will be used later. We usually see a decreasing trend
# with abundance, possibly plateauing off at very large abundances.
# Little variability among replicates
#
par(mfrow=c(1,1))
plotBCV(y)

# QL
#
# We can augment this model with quasi-likelihood (QL) methods. 
# The QL dispersions account for window-specific variability
#
# We fit a generalized linear model (GLM) to the counts for each window, 
# where the QL dispersion for that window is estimated from the GLM deviance.
#
fit <- glmQLFit(y, design, robust=TRUE)
summary(fit$var.post)
summary(fit$df.prior)

# The effect of EB (Empirical Bayes) stabilisation/shrinkages can be visualized by examining the
# biological coefficient of variation (for the NB dispersion) and the quarter-root deviance 
# (for the QL dispersion).
#
par(mfrow=c(1,2))
o <- order(y$AveLogCPM)
plot(y$AveLogCPM[o], sqrt(y$trended.dispersion[o]), type="l", lwd=2,
     ylim=c(0, 1), xlab=expression("Ave."~Log[2]~"CPM"),
     ylab=("Biological coefficient of variation"))
plotQLDisp(fit)

#
# Examining replicate similarity with MDS plots                                   ####
#
# As a quality control measure, the window counts can be used to examine the similarity of 
# replicates through multi-dimensional scaling (MDS) plots 
# https://www.statisticshowto.com/multidimensional-scaling/
#
# Multidimensional scaling is a visual representation of distances or dissimilarities between sets of objects.
# The term scaling comes from psychometrics, where abstract concepts (“objects”) are assigned numbers according to a rule (Trochim, 2006). 
# You can also think of “scaling” as the fact that you’re essentially scaling down the data (i.e. making it simpler by creating
# lower-dimensional data). Data that is scaled down in dimension keeps similar properties. For example, two data points that are
# close together in high-dimensional space will also be close together in low-dimensional space (Martinez, 2005). 
# PCA is another similar tool, but while MDS uses a similarity matrix to plot the graph, PCA uses the original data.
#
# The distance between each pair of libraries is computed as the square root of the mean squared 
# log-fold change across the top set of bins with the highest absolute log-fold changes. 
#
# A small top set visualizes the most extreme differences whereas
# a large set visualizes overall differences. 
# Checking a range of top values may be useful when the scope of DB is unknown.
#
### MDS plot with two dimensions for all samples in the SSP data set
#
# NO batch effect between replicates (GOOD)
par(mfrow=c(1,1))
pdf(file = "MDS_filt_fc3.pdf", 10,10)
pdf(file = "MDS_filt_fc4.pdf", 10,10)
par(mar=c(5,7,2,3))
plotMDS(cpm(y, log=TRUE), top=10000, labels=antibody,
col=c("grey37", "green", "red", "blue")[as.integer(antibody)], xlim=c(-4,4), cex=2, cex.axis=2,
     cex.lab=3)
dev.off()

#--------------------------------------------------------------------------------------
# Testing for DB                                                                  ####
#
# Using the QL F-test.
#
#
fdr <- 0.05
fdr <- 0.01
flag_ignore_down <- 1 # (for comparisons with Input)
design

### 1. SSP vs Input
contrast <- makeContrasts(SSP-Input, levels=design)
comparison <- paste0("SSP_vs_Input_filt_fc_",filtering.min.fc,"_minq.20_tol_100_fdr_",fdr)
comparison
flag_ignore_down <- 1 
cond <- strsplit(gsub('.{35}$', '', comparison), split ="_vs_")[[1]][1]
ctrl <- strsplit(gsub('.{35}$', '', comparison), split ="_vs_")[[1]][2]
cond
ctrl

### 2. SlycSSP2 vs Input
contrast <- makeContrasts(SlycSSP2-Input, levels=design)
comparison <- paste0("SlycSSP2_vs_Input_filt_fc_",filtering.min.fc,"_minq.20_tol_100_fdr_",fdr)
comparison
flag_ignore_down <- 1 
cond <- strsplit(gsub('.{35}$', '', comparison), split ="_vs_")[[1]][1]
ctrl <- strsplit(gsub('.{35}$', '', comparison), split ="_vs_")[[1]][2]
cond
ctrl

### 3. SpimSSP2 vs Input
contrast <- makeContrasts(SpimSSP2-Input, levels=design)
comparison <- paste0("SpimSSP2_vs_Input_filt_fc_",filtering.min.fc,"_minq.20_tol_100_fdr_",fdr)
comparison
flag_ignore_down <- 1 
cond <- strsplit(gsub('.{35}$', '', comparison), split ="_vs_")[[1]][1]
ctrl <- strsplit(gsub('.{35}$', '', comparison), split ="_vs_")[[1]][2]
cond
ctrl

### 4. SSP vs SlycSSP2
contrast <- makeContrasts(SSP-SlycSSP2, levels=design)
comparison <- paste0("SSP_vs_SlycSSP2_filt_fc_",filtering.min.fc,"_minq.20_tol_100_fdr_",fdr)
comparison
flag_ignore_down <- 0
cond <- strsplit(gsub('.{35}$', '', comparison), split ="_vs_")[[1]][1]
ctrl <- strsplit(gsub('.{35}$', '', comparison), split ="_vs_")[[1]][2]
cond
ctrl

### 5. SSP vs SpimSSP2
contrast <- makeContrasts(SSP-SpimSSP2, levels=design)
comparison <- paste0("SSP_vs_SpimSSP2_filt_fc_",filtering.min.fc,"_minq.20_tol_100_fdr_",fdr)
comparison
flag_ignore_down <- 0
cond <- strsplit(gsub('.{35}$', '', comparison), split ="_vs_")[[1]][1]
ctrl <- strsplit(gsub('.{35}$', '', comparison), split ="_vs_")[[1]][2]
cond
ctrl

### 6. SpimSSP2 vs SlycSSP2
contrast <- makeContrasts(SpimSSP2-SlycSSP2, levels=design)
comparison <- paste0("SpimSSP2_vs_SlycSSP2_filt_fc_",filtering.min.fc,"_minq.20_tol_100_fdr_",fdr)
comparison
flag_ignore_down <- 0
cond <- strsplit(gsub('.{35}$', '', comparison), split ="_vs_")[[1]][1]
ctrl <- strsplit(gsub('.{35}$', '', comparison), split ="_vs_")[[1]][2]
cond
ctrl

#
contrast
#
res <- glmQLFTest(fit, contrast=contrast)
dim(res$table)
head(res$table)
dim(res$table)

# Windows less than 100 bp apart are clustered into regions (Lun and Smyth 2014) with a maximum cluster width of 5 kbp. 
# We then control the region-level FDR by combining per-window  p-values using Simes’ method (Simes 1986).
#
# This computes a single combined p-value for each cluster/region, representing the evidence against the global null.
#
tol.res <- 100
max.region <- 5000
merged <- mergeResults(filtered.data, res$table, tol=tol.res, 
                       merge.args=list(max.width=max.region))
merged$regions

# We also get some other statistics, like the total number of windows in each cluster as well as
# the number that are changing substantially in each direction.
tabcom <- merged$combined
head(tabcom)
dim(tabcom)

# Significant DB regions                                                          ####
#
is.sig <- tabcom$FDR <= fdr
summary(is.sig)

# All significant regions have increased SSP binding 
table(tabcom$direction[is.sig])

# Direction according the best window in each cluster
tabbest <- merged$best
log2fc.threshold <- 0
is.sig.pos <- (tabbest$rep.logFC > log2fc.threshold)[is.sig]
summary(is.sig.pos)

#
# Save the results to file in the form of a serialized R object                    #### 
out.ranges <- merged$regions
mcols(out.ranges) <- DataFrame(tabcom,
                               best.pos=mid(ranges(rowRanges(filtered.data[tabbest$rep.test]))),
                               best.logFC=tabbest$rep.logFC)
out.ranges
length(out.ranges)
outputfile=paste0(outputdir,"/", comparison,"_results.rds")
outputfile
saveRDS(file=outputfile, out.ranges)

# Save the significant regions in BED format                                      ####
assembly <- "S100_v2"
simplified <- out.ranges[is.sig]
simplified$score <- simplified$best.logFC   #-10*log10(simplified$FDR)
## remove "invalid" regions with input>IP
if(flag_ignore_down) {
    to.remove=(simplified$direction%in%c("down","mixed"))|(simplified$best.logFC<0)
    tt=table(factor(to.remove,levels=c(TRUE,FALSE)))
    cat("Ignoring invalid regions (with input>IP). Nb invalid regions:",tt["TRUE"]," (",signif(100*tt["TRUE"]/sum(tt),2),"%), nb valid regions :",tt["FALSE"],"(",signif(100*tt["FALSE"]/sum(tt),2),"%)\n")
    simplified=simplified[!to.remove]
}
length(simplified)     # nb. of pos.only/or both pos and neg. fc peaks
mcols(simplified)
range(simplified)
outputfile=paste0(outputdir,"/", comparison,"_results.bed")
cat("creating",outputfile,"\n")
trackline=new("BasicTrackLine",db=assembly,name=comparison,description=comparison)
trackline
### Export BED file                                                              ####
export(con=outputfile, object=simplified, trackLine=trackline)

# Prepare narrowPeak format                                                      ####
simplified_df <- data.frame(simplified)
head(simplified_df)
dim(simplified_df)
simplified_df$strand = '.'
head(simplified_df)
simplified_df$start <- simplified_df$start -1
simplified_df$name <- gsub('.{35}$', '', comparison)
simplified_df$score <- simplified$best.logFC
simplified_df$signalValue <- simplified$best.logFC*10 
simplified_df$pValue <- -log10(simplified$PValue)
simplified_df$qValue <- -10*log10(simplified$FDR)
simplified_df$peak <- (simplified_df$end-simplified_df$start)/2
head(simplified_df)
simplified_df <- simplified_df[, -c(4,6:15)]
head(simplified_df)
simplified_df <- simplified_df[, c(1:3,6,5,4,7:10)] 
#if(flag_ignore_down) {
#    simplified_df <- simplified_df[, c(1:4,6,5,7:10)]
#}
head(simplified_df)
colnames(simplified_df)[1] <- "chrom"
head(simplified_df)

trackline=paste0("track ", "type=narrowPeak ", "db=",assembly," name=",comparison," description=",comparison,"\n")
trackline
# Export narrowPeak BED                                                           ####
outputfile=paste0(outputdir,"/", comparison,"_results_narrowPeak.bed")
cat("creating",outputfile,"\n")
cat(NULL,file=outputfile)
cat(trackline, file=outputfile, append=TRUE)
write.table(simplified_df, file=outputfile, quote = F, sep = "\t", col.names = F, row.names = F, append=TRUE)

# Save table result from glmQLFTest                                               ####
outputfile=paste0(outputdir,"/",comparison,"_filtered_data_stats.rds")
cat("creating",outputfile,"\n")
saveRDS(file=outputfile, res$table)

# Save full filtered.data object (with res$table in rowData)                      ####
filtered.data.tmp <- filtered.data
rowData(filtered.data.tmp) <- cbind(rowData(filtered.data.tmp), res$table)

outputfile=paste0(outputdir,"/",comparison,"_filtered_data.rds")
cat("creating",outputfile,"\n")
saveRDS(file=outputfile, filtered.data.tmp)

# Save significant DB regions as tsv                                              ####
simplified <- out.ranges[is.sig]
simplified=as.data.frame(simplified)
colnames(simplified)=gsub("^seqnames","chr",colnames(simplified))
simplified=simplified[,!colnames(simplified)%in%c("strand","mean.logFC")]
dim(simplified)
outputfile=paste0(outputdir,"/",comparison,"_results.tsv")
cat("creating",outputfile,"\n")
write.table(simplified,sep="\t",quote=FALSE,row.names=FALSE,file=outputfile)

# Save significant DB regions as xlsx                                             ####
outputfile=paste0(outputdir,"/",comparison,"_results.xlsx")
cat("creating",outputfile,"\n")
write_xlsx(simplified,path=outputfile)

### Save general settings                         ####
outputfile=paste0(outputdir,"/","readParams.rds")
cat("creating",outputfile,"\n")
saveRDS(file=outputfile, param)

window.spacing = 50   # default
settings=list(
    outputdir=outputdir,
    inputdir=inputdir,
    assembly=assembly,
    standard.chr=standard.chr,
    flag_paired=1,
    flag_nodup=0,
    fdr.threshold=fdr,
    log2fc.threshold=log2fc.threshold,
    window.width=window.width,
    window.spacing=window.spacing,
    filtering.min.fc=filtering.min.fc,
    metadatafile=metadatafile,
    fragment.length=frag.len,
    max.fragment.length=1000,
    normalization.method.option="composition_biases"
    #comparison_condition1=comparison_condition1,
    #comparison_antibody1=comparison_antibody1,
    #comparison_condition2=comparison_condition2,
    #comparison_antibody2=comparison_antibody2,
    #comparison=comparison
)
outputfile=paste0(outputdir,"/","settings.rds")
cat("creating",outputfile,"\n")
saveRDS(file=outputfile,settings)


# Export normalized data to bigWig                                                ####
#
# DO only once
normalization_method <- "composition_biases"
if(normalization_method=="composition_biases")
{
    output.normalized.mat <- calculateCPM(filtered.data,use.norm.factors=TRUE,log=FALSE)
}
if(normalization_method=="trended_biases")
{
    output.normalized.mat <- calculateCPM(filtered.data,use.offsets=TRUE,log=FALSE)
}


outputdir.tmp=paste0(outputdir,"/","binned_filtered_fc4_normalized_bigwig/")
outputdir.tmp
dir.create(outputdir.tmp, showWarnings = FALSE,recursive=TRUE)

# Create out.data with ranges from rowRanges(filtered.data) but reducing window size to spacing
out.data=rowRanges(filtered.data)
ranges(out.data)=IRanges(start=start(out.data)+(window.width-window.spacing)/2,width=window.spacing)
##Note:
## metadata(filtered.data)$width==window.width
## metadata(filtered.data)$spacing==window.spacing

output.normalized=list()
for(i in seq_along(colData(filtered.data)[,"bam.files"])) {
    sample=colData(filtered.data)[,"bam.files"][i]
    sample.name=paste(metadata[sample,c("condition","replicate")],collapse="_")
    filename=paste0(outputdir.tmp,"/",sample.name,"_binned_filtered_normalized.bw")
    
    mcols(out.data)=data.frame(score=output.normalized.mat[,i])
    output.normalized[[i]]=out.data
    cat("creating",filename,"\n")
    export.bw(out.data,con=filename)
}

# Done 
cat("Done 1\n")

#--------------------------------------------------------------------------------------
# Annotation and visualization                                                     ####
#--------------------------------------------------------------------------------------
#
# The makeTxDbFromGFF function allows the user to make a TxDb object from transcript 
# annotations available as a GFF3 or GTF file.

## Load gff file for assembly of interest
gffFile <- "/Users/giovanna/tomato_bZIP_sebastien_soyk/data/genome/S100_v2.1.0/SollycSweet-100_genes_v2.1.1.gff3"
gffFile <- "/Users/giovanna/tomato_bZIP_sebastien_soyk/data/genome/S100_v2.1.0/SollycSweet-100_genes_v2.1.1.new.gff3"

## Generate annotation object
txdb <- makeTxDbFromGFF(file=gffFile,
                        dataSource="Alonge et al., 2021",
                        organism="Solanum lycopersicum", dbxrefTag="Alias")

keytypes(txdb)
columns(txdb)

# Install custom library
# install.packages("/Users/giovanna/tomato_bZIP_sebastien_soyk/data/GO_term/org.SlycopersicumNatalia.eg.db", repos = NULL, type = "source")
library(org.SlycopersicumNatalia.eg.db)
org.SlycopersicumNatalia.eg.db

keytypes(org.SlycopersicumNatalia.eg.db)
columns(org.SlycopersicumNatalia.eg.db)

egid <- head(keys(org.SlycopersicumNatalia.eg.db, "GID"))
egid <- keys(org.SlycopersicumNatalia.eg.db, "GID")
select(org.SlycopersicumNatalia.eg.db, egid, c("SYMBOL", "GENENAME"), "GID")
select(org.SlycopersicumNatalia.eg.db, egid, c("CHROMOSOME", "SYMBOL", "GENENAME", "GO"), "GID")

gid <- head(keys(txdb, "GENEID"))
gid <- keys(txdb, "GENEID")
select(txdb, gid, c("TXNAME", "TXID"), "GENEID")
select(txdb, gid, c("CDSID", "CDSNAME", "EXONID", "EXONNAME", "TXID",  "TXNAME"), "GENEID")

# Adding gene-based annotation                                                    ####
#
# Here, the promoter region of each gene is defined as some interval 3 kbp up- and 1 kbp downstream of the TSS for that gene.
# Any exonic features within dist on the left or right side of each supplied region will also be reported.
#
differential.binding.regions <- out.ranges[is.sig]
### Remove down and mixed peaks (if DB is against Input)
if(flag_ignore_down) {
    to.remove=(differential.binding.regions$direction%in%c("down","mixed"))|(differential.binding.regions$best.logFC<0)
    tt=table(factor(to.remove,levels=c(TRUE,FALSE)))
    cat("Ignoring invalid regions (with input>IP). Nb invalid regions:",tt["TRUE"]," (",signif(100*tt["TRUE"]/sum(tt),2),"%), nb valid regions :",tt["FALSE"],"(",signif(100*tt["FALSE"]/sum(tt),2),"%)\n")
    differential.binding.regions=differential.binding.regions[!to.remove]
}

length(differential.binding.regions)

anno <- detailRanges(differential.binding.regions, txdb=txdb,
                     orgdb=org.SlycopersicumNatalia.eg.db, key.field="GID", promoter=c(3000, 2000))
#
# Character vectors of compact string representations are provided to summarize the features overlapped 
# by each supplied region. Each pattern contains GENE|STRAND|TYPE to describe the strand and overlapped 
# features of that gene. Exons are labelled as E, promoters are P and introns are I. 
# For left and right, TYPE is replaced by DISTANCE. This indicates the gap (in base pairs) between the
# supplied region and the closest non-overlapping exon of the annotated feature. 
# All of this annotation can be stored in the metadata of the GRanges object for later use.
head(anno$overlap)

differential.binding.regions$overlap <- anno$overlap

# While the string representation saves space in the output, it is not easy to work with. 
# If the annotation needs to manipulated directly, users can obtain it from the detailRanges() command
# by not specifying the regions of interest. This can then be used for interactive manipulation, e.g., 
# to identify all genes where the promoter contains DB sites.
anno.ranges <- detailRanges(txdb=txdb, 
                            orgdb=org.SlycopersicumNatalia.eg.db, key.field="GID")
anno.ranges

length(anno.ranges)
length(differential.binding.regions)
head(differential.binding.regions)

# Save binding sites Annotation                                                    ####
outputfile=paste0(outputdir,"/",comparison,"_db_anno_results.rds")
cat("creating",outputfile,"\n")
saveRDS(file=outputfile, differential.binding.regions)

# Save binding Annotation as tsv                                                   ####
anno=as.data.frame(differential.binding.regions)
head(anno)
colnames(anno)=gsub("^seqnames","chr",colnames(anno))
anno=anno[,!colnames(anno)%in%c("strand","mean.logFC")]
dim(anno)
head(anno)
outputfile=paste0(outputdir,"/",comparison,"_db_anno_results.tsv")
cat("creating",outputfile,"\n")
write.table(anno,sep="\t",quote=FALSE,row.names=FALSE,file=outputfile)

### Simple visualization of genomic coverage                                       ####
#
# for SSP - sample 3
#
# Here, coverage is visualized as the number of reads covering each base pair in the interval of interest. 
# Specifically, the reads-per-million is shown to allow comparisons between libraries of different size. 
# The plots themselves are constructed using methods from the Gviz package. The blue and red tracks represent
# the coverage on the forward and reverse strands, respectively. Strong strand bimodality is consistent 
# with a genuine TF binding site. For paired-end data, coverage can be similarly plotted for fragments, 
# i.e., proper read pairs.
#
cur.region <- GRanges(c("chr3:64800-64856", "chr6:38000000-38151000", "chr11:57500-57550"))
cur.region <- GRanges("chr3", IRanges(64840000, 64855000))
cur.region <- GRanges("chr6", IRanges(38000000, 38151000))
cur.region <- GRanges("chr11", IRanges(57500000, 57700000))
#library(Gviz)
gax <- GenomeAxisTrack(col="black", fontsize=15, size=7)  # size increases the top margin

options(ucscChromosomeNames=FALSE)
greg <- GeneRegionTrack(txdb, showId=TRUE, ucscChromosomeNames=FALSE,
                        geneSymbol=TRUE, name="", background.title="transparent")

# mapIds gets the mapped ids (column) for a set of keys that are of a particular keytype.
# Usually returned as a named character vector.
symbols <- unlist(mapIds(org.SlycopersicumNatalia.eg.db, gene(greg), "GID", "GID",
                         multiVals = "first"))
symbol(greg) <- symbols[gene(greg)]

### Print coverage plot of top-50/100 diff-binding peaks ordered by logFC             ####
###
### Read binding regions input file for plots                                         ####
inputfile=paste0(outputdir,"/",comparison,"_db_anno_results.rds")
cat("reading",inputfile,"\n")
differential.binding.regions <- readRDS(file=inputfile)
#
### Read filtered data for calculating normalized counts                              ####
filtered.data <- readRDS('filtered_data.rds')
data.normalized.all.mat <- calculateCPM(filtered.data,use.norm.factors=TRUE,log=FALSE)
# create list
data.normalized.all=list()
out.data=rowRanges(filtered.data)
for(i in 1:ncol(filtered.data)) {
    mcols(out.data)=data.frame(score=data.normalized.all.mat[,i])
    data.normalized.all[[i]]=out.data
}

# Order diff-binding peaks by logFC                                                ####
o <-order(abs(differential.binding.regions$best.logFC), decreasing=TRUE)
#o <- order(differential.binding.regions$PValue)    
o
length(o)
# Plot Sig DB track , raw tracks, and normalized tracks                            ####
signif_DB_track <- AnnotationTrack(differential.binding.regions, name="DB", col.axis="black", col.title="black",fill="yellow",cex.title=0.5,rotation.title=0)
signif_DB_track
### Open PDF file                                                                      ####
outputfile=paste0(outputdir,"/",comparison,"_results.pdf")
cat("creating",outputfile,"\n")
cairo_pdf(outputfile,10,2*nrow(metadata),onefile =TRUE) 
max.plots <- 100
for (i in seq(max.plots)) {
    cur.region <-differential.binding.regions[o[i]] 
    ### Extend bs region (-5000, + 3000)
    start(cur.region) <- start(cur.region) - 5000
    end(cur.region) <- end(cur.region) + 3000
    ### raw data
    collected <- list()
    for (j in seq_along(metadata$filename)) { 
        reads <- extractReads(metadata$filename[j], cur.region, param=param, ext=list(init.ext=rep(NA,length(metadata$filename[j])), final.ext=settings[["fragment.length"]]))
        cov <- as(coverage(reads)/lib.sizes[j], "GRanges")
        collected[[j]] <- DataTrack(cov, type="histogram", lwd=0, ylim=c(0,50), # 15000 for chr3
                                    name=metadata$sample[j], col.axis="black", col.title="black",
                                    fill="darkgray", col.histogram=NA, genome=assembly)
    }
    ##normalized
    collected.normalized=list()
    data.normalized=lapply(1:nrow(metadata),function(j){
            sample=metadata$filename[j]
            sample.name=metadata$sample[j]
            subsetByOverlaps(data.normalized.all[[j]], cur.region) # data normalized is a list of 8 Granges norm objects
    })
    maxnorm=max(sapply(unlist(data.normalized,recursive=FALSE),function(x){max(x$score)}))
    if(!is.finite(maxnorm))maxnorm=1
    for (j in 1:nrow(metadata)) {
        window.width=window.width
        window.spacing=window.spacing
        sample=metadata$filename[j]
        sample.name=metadata$sample[j]
        collected.normalized=c(collected.normalized, DataTrack(data.normalized[[j]], type="histogram", lwd=0, ylim=c(0,maxnorm),
                                  name=metadata$sample[j], col.axis="black", col.title="black",
                                  fill="blue", col.histogram=NA,genome=assembly))
    }
    #plotTracks(c(gax, signif_DB_track, collected, collected.normalized, greg), chromosome=as.character(seqnames(cur.region)), transcriptAnnotation="symbol",
    #           from=start(cur.region)-5000, to=end(cur.region)+2000, margin = 6, innerMargin = 3)
    plotTracks(c(gax, signif_DB_track, collected, collected.normalized, greg), chromosome=as.character(seqnames(cur.region)), transcriptAnnotation="symbol",
               from=start(cur.region), to=end(cur.region), margin = 6, innerMargin = 3)
}
dev.off()
# Done
cat("Done 2\n")