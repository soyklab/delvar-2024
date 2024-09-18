
## Experiment description
# summary of meristems/stages that were collected
# - WT and mutants in M82 background
# - meristems were collected at transition stage (TM) of meristem maturation
# 
# sample	name	n_meristems
# S01	WT_TM_1	16
# S02	WT_TM_2	15
# S03	WT_TM_3	15
# S04	ssp_TM_1	15
# S05	ssp_TM_2	14
# S06	ssp_TM_3	14
# S07	ssp2_TM_1	15
# S08	ssp2_TM_2	15
# S09	ssp2_TM_3	14
# S10	sspssp2_TM_1	22
# S11	sspssp2_TM_2	22
# S12	sspssp2_TM_3	22
# 
# Brief summary of raw read analysis
# - reads were aligned to M82_v1.1.0 reference using STAR
# - reads were counted using htseq-count and the SollycM82_genes_v1.1.1.gff3 annotation



## Set environment

# set directories
setwd("./")
home <- getwd()
htseqdir <- paste(home,"htseq/",sep="/")

# create dated directory for DESeq2 file exports
date <- format(Sys.time(), "%Y%m%d")
dir.create(file.path(home,paste0("DESeq2_",date)), showWarnings = FALSE)
dir.create(file.path(home,paste0("DESeq2_",date,"/plots")), showWarnings = FALSE)
dir.create(file.path(home,paste0("DESeq2_",date,"/tables")), showWarnings = FALSE)


### Load packages

library("DEGreport")
library("tidyverse")
library("DESeq2")
library("apeglm")
library("vsn")
library("dplyr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("EnhancedVolcano")
library("viridis")
library("DEGreport")
library("ggfortify")
library("PCAtools")
require(vidger)

# Load description and position file (M82_v1.0)
descdat <- read.delim(file="/FILEPATHTO/SollycM82_v1.0.fasta.descriptions", header=FALSE, sep='\t')
colnames(descdat) <- c("GeneName", "M82_v1.0_geneid", "M82_v1.0_chr", "M82_v1.0_start" ,"M82_v1.0_end" ,"M82_v1.0_strand", "M82_v1.0_description")
desc <- descdat[ -c(1, 3:6) ]
colnames(desc) <- c("geneid", "description")
head(desc)
desc_tb <- as_tibble(desc)

## Data import and preparation for DESeq2
### Process HTSeq input data

# read sample names
samples <- read.table(file.path(htseqdir, "samples"), sep = '\t', header = FALSE)
colnames(samples) <- c("sample_id", "sample_name")
samples

# clean up htseq output (remove "mRNA:", "gene-", "rna-", ^NC_" part and "__" entries at the tail of the HTseq count file)

for (F in c(grep(".htseq$", list.files(htseqdir),value=TRUE))) { # the $ anchor matches the end of the line
  print(paste("process sample", F))
  counts <- read.table(file.path(htseqdir, F), header = FALSE, row.names = 1)
  print(paste("lines before processing: ", nrow(counts)))
  rownames(counts) <- gsub("mRNA:","",rownames(counts)) # remove mRNA: part of gene ids
  rownames(counts) <- gsub("gene-","",rownames(counts)) # remove gene- part of gene ids
  rownames(counts) <- gsub("rna-","",rownames(counts)) # remove rna- part of gene ids
  rownames(counts) <- gsub("^NC_","",rownames(counts)) # remove NC_ part of gene ids
  counts_processed <- subset(counts, !grepl("^__", rownames(counts))) # delete lines with "__" (i.e. last lines of HTseq files that contain alignment parameters)
  print(paste("lines after processing: ", nrow(counts_processed)))
  write.table(counts_processed, (paste0(htseqdir,F,".processed")), append = FALSE, sep = " ", row.names = TRUE, col.names = FALSE, quote = FALSE)
}

### Import HTSeq input data

# Import individual HTseq data files
sampleFiles <- grep(".htseq.processed",list.files(htseqdir),value=TRUE)
sampleName <- samples$sample_name
sampleCondition <- gsub(".*_(.*)_.*..*", "\\1", sampleName) # by meristem stage (string between first and second underscore in filename)
sampleStrain <- gsub("(.*)_.*_.*..*", "\\1", sampleName) # by genotype (string before first underscore in filename)
sampleGroup <- factor(paste0(sampleStrain,"_",sampleCondition))


### Build SampleTable
sampleTable <- data.frame(sampleName = sampleName,
                          fileName = sampleFiles,
                          condition = sampleCondition,
                          strain = sampleStrain,
                          group = sampleGroup)

sampleTable$condition <- factor(sampleTable$condition)
sampleTable$strain <- factor(sampleTable$strain)
sampleTable$group <- factor(sampleTable$group)
row.names(sampleTable) <- sampleTable$sampleName # important for "degPatterns" function when extracting clusters from LRT analysis
sampleTable # make sure that sampleName and fileName match/fit

### Build the DESeqDataSet
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = htseqdir,
                                  design= ~ strain)
dds

### DESeqDataSet Design
design(dds)

## Factor levels to set order for later comparisons
dds$strain <- factor(dds$strain, levels = c("WT","ssp", "ssp2", "sspssp2"))
dds$strain

### Prefiltering DESeqDataSet
# Prefilter by removing rows with very few reads (BUT messes with combining CPM and DESEQ files)
keep <- rowSums(counts(dds)) >= 5 # removes rows with equal or less than 1 reads 
dds <- dds[keep,]
dds

## Conduct differential expression analysis
dds <- DESeq(dds)

## Check the size factors
sizeFactors(dds)
## Total number of raw counts per sample
colSums(counts(dds))
## Total number of normalized counts per sample
colSums(counts(dds, normalized=TRUE))



# Check replicates reproducibility

# Convert normalized_counts to a data frame 
ncounts <- counts(dds, normalized=TRUE)

# check replicates reproducibility
c <- cor(as.matrix(ncounts), method="spearman")
d_c <- as.dist(1-c) 
# hierarchical clustering
fit_h <- hclust(d_c, method="ward.D") 
# display dendogram
plot(fit_h) +
  abline(h=0.2, col = "red")


# Heatmap of the sample-to-sample distances

# extract transformed values using vst function
vsd <- vst(dds, blind=FALSE)


#  Apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances.
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$strain, vsd$condition, sep=":")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
x <- pheatmap(sampleDistMatrix,
              clustering_distance_rows=sampleDists,
              clustering_distance_cols=sampleDists,
              col=colors)

save_pheatmap_pdf <- function(x, filename, width=5, height=4) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(x, (file.path(paste0("./DESeq2_",date,"/plots/sample_distances.pdf"))))


# Principal component plot of the samples
# plot transformed values in PCA
ntop=500

plotPCA(vsd, intgroup=c("strain", "condition"), ntop = ntop) +
  theme_light() +
  ggtitle("n(top genes) =", ntop)

ggsave((file.path(paste0("./DESeq2_",date,"/plots/PCA_top",ntop,".pdf"))), width = 5, height = 4, useDingbats=FALSE)


# Additional PCA plots with all genes
vst <- assay(vst(dds))
p <- pca(vst, metadata=sampleTable, removeVar = 0.1)
n_pca <- nrow(vst)

screeplot(p, axisLabSize = 18, titleLabSize = 22)
ggsave((file.path(paste0("./DESeq2_",date,"/plots/PCA_all_screeplot.pdf"))), width = 5, height = 4, useDingbats=FALSE)

biplot(p, showLoadings = FALSE,
       labSize = 3, pointSize = 5, sizeLoadingsNames = 5, title=paste("n(genes) =", n_pca))
ggsave((file.path(paste0("./DESeq2_",date,"/plots/PCA_all_biplot.pdf"))), width = 5, height = 4, useDingbats=FALSE)


# Check which names and comparisons can be selected for the contrasts
resultsNames(dds) 


# Define contrasts
contrast_ssp_WT <- c("strain", "ssp", "WT") # contrast <- c("condition", "level_to_compare", "base_level")
contrast_ssp2_WT <- c("strain", "ssp2", "WT") # contrast <- c("condition", "level_to_compare", "base_level")
contrast_sspssp2_WT <- c("strain", "sspssp2", "WT") # contrast <- c("condition", "level_to_compare", "base_level")

## Extract results for mutants vs control
# Perform independent filtering based on the mean of normalized counts for each gene, optimizing the number of genes which will have an adjusted p value below a given FDR cutoff, "alpha"

resultsNames(dds)

# for design= ~ strain)
res_ssp_vs_WT <- results(dds, contrast = contrast_ssp_WT , alpha =0.1)
res_ssp2_vs_WT <- results(dds, contrast = contrast_ssp2_WT , alpha =0.1)
res_sspssp2_vs_WT <- results(dds, contrast = contrast_sspssp2_WT , alpha =0.1)

#print summary for results
summary(res_ssp_vs_WT)
summary(res_ssp2_vs_WT)
summary(res_sspssp2_vs_WT)

# MA plots
plotMA(res_ssp_vs_WT, ylim=c(-4,4))
plotMA(res_ssp2_vs_WT, ylim=c(-4,4))
plotMA(res_sspssp2_vs_WT, ylim=c(-4,4))


## DE gene selection from contrasts
### Set thresholds for analyses

# Create a tibble of results
res_ssp_vs_WT_tb <- res_ssp_vs_WT %>%
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  left_join(desc_tb, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE)

res_ssp2_vs_WT_tb <- res_ssp2_vs_WT %>%
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  left_join(desc_tb, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE)

res_sspssp2_vs_WT_tb <- res_sspssp2_vs_WT %>%
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  left_join(desc_tb, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE)


# Print results to file
write.csv(as.data.frame(res_ssp_vs_WT_tb), file=(paste0("DESeq2_",date,"/tables/res_ssp_vs_WT.csv")), row.names = FALSE)
write.csv(as.data.frame(res_ssp2_vs_WT_tb), file=(paste0("DESeq2_",date,"/tables/res_ssp2_vs_WT.csv")), row.names = FALSE)
write.csv(as.data.frame(res_sspssp2_vs_WT_tb), file=(paste0("DESeq2_",date,"/tables/res_sspssp2_vs_WT.csv")), row.names = FALSE)


###########################
# Set cutoffs
padj.cutoff <- 0.05
padj.cutoff
lfc.cutoff <- 0.58
lfc.cutoff
###########################


sum(res_ssp_vs_WT$padj < padj.cutoff, na.rm=TRUE)
sum(res_ssp2_vs_WT$padj < padj.cutoff, na.rm=TRUE)
sum(res_sspssp2_vs_WT$padj < padj.cutoff, na.rm=TRUE)


# Filter by padj and lfc
sig_ssp_vs_WT <- res_ssp_vs_WT_tb %>%
  filter(padj < padj.cutoff) %>%
  filter(!between(log2FoldChange, -lfc.cutoff, lfc.cutoff))
sig_ssp2_vs_WT <- res_ssp2_vs_WT_tb %>%
  filter(padj < padj.cutoff) %>%
  filter(!between(log2FoldChange, -lfc.cutoff, lfc.cutoff))
sig_sspssp2_vs_WT <- res_sspssp2_vs_WT_tb %>%
  filter(padj < padj.cutoff) %>%
  filter(!between(log2FoldChange, -lfc.cutoff, lfc.cutoff))


### Data Visualization
# Volcano plots

library(gridExtra)
library(grid)

p1 <- EnhancedVolcano(res_ssp_vs_WT,
                      # lab = rownames(res_ssp_vs_WT),
                      lab = NA,
                      title = 'ssp vs WT',
                      subtitle = paste0("lfc.cutoff = ", lfc.cutoff, "; padj.cutoff = ", padj.cutoff ),
                      caption = paste0("significant = ", nrow(sig_ssp_vs_WT), "; total = ", nrow(res_ssp_vs_WT)),
                      x = "log2FoldChange",
                      y = "padj",
                      # selectLab = genes,
                      drawConnectors = TRUE,
                      FCcutoff = lfc.cutoff,
                      pCutoff = padj.cutoff)

p2 <- EnhancedVolcano(res_ssp2_vs_WT,
                      # lab = rownames(res_ssp2_vs_WT),
                      lab = NA,
                      title = 'ssp2 vs WT',
                      subtitle = paste0("lfc.cutoff = ", lfc.cutoff, "; padj.cutoff = ", padj.cutoff ),
                      caption = paste0("significant = ", nrow(sig_ssp2_vs_WT), "; total = ", nrow(res_ssp2_vs_WT)),
                      x = "log2FoldChange",
                      y = "padj",
                      # selectLab = genes,
                      drawConnectors = TRUE,
                      FCcutoff = lfc.cutoff,
                      pCutoff = padj.cutoff)

p3 <- EnhancedVolcano(res_sspssp2_vs_WT,
                      # lab = rownames(res_sspssp2_vs_WT),
                      lab = NA,
                      title = 'sspssp2 vs WT',
                      subtitle = paste0("lfc.cutoff = ", lfc.cutoff, "; padj.cutoff = ", padj.cutoff ),
                      caption = paste0("significant = ", nrow(sig_sspssp2_vs_WT), "; total = ", nrow(res_sspssp2_vs_WT)),
                      x = "log2FoldChange",
                      y = "padj",
                      # selectLab = genes,
                      drawConnectors = TRUE,
                      FCcutoff = lfc.cutoff,
                      pCutoff = padj.cutoff)



# Print as grid
grid.arrange(p1, p2, p3, ncol=3)

# run in console:
pdf((file.path(paste0("./DESeq2_",date,"/plots/volcano_padj",padj.cutoff,"_lfc",lfc.cutoff,".pdf"))), width = 14, height = 6, useDingbats=FALSE)
grid.arrange(p1, p2, p3, ncol=3)
dev.off() # Close the file


# Print significant DEGs to file
write.csv(as.data.frame(sig_ssp_vs_WT), file=(paste0("DESeq2_",date,"/tables/sig_ssp_vs_WT_padj",padj.cutoff,"_lfc", lfc.cutoff,".csv")), row.names = FALSE)
write.csv(as.data.frame(sig_ssp2_vs_WT), file=(paste0("DESeq2_",date,"/tables/sig_ssp2_vs_WT_padj",padj.cutoff,"_lfc", lfc.cutoff,".csv")), row.names = FALSE)
write.csv(as.data.frame(sig_sspssp2_vs_WT), file=(paste0("DESeq2_",date,"/tables/sig_sspssp2_vs_WT_padj",padj.cutoff,"_lfc", lfc.cutoff,".csv")), row.names = FALSE)

### Normalized counts

# Convert normalized_counts to a data frame 
ncounts <- counts(dds, normalized=TRUE) %>% #DESeq2 creates a matrix when you use the counts() function
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble()
head(ncounts)

write.csv(as.data.frame(ncounts), file=(paste0("DESeq2_",date,"/tables/ncounts.cpm.csv")), row.names = FALSE)

# Filter for significantly DEGs
padj.cutoff
lfc.cutoff

ncounts_sig <- ncounts %>%
  filter(ncounts$geneid %in% sig_ssp_vs_WT$geneid | ncounts$geneid %in% sig_ssp2_vs_WT$geneid | ncounts$geneid %in% sig_sspssp2_vs_WT$geneid)
nrow(ncounts_sig)

write.csv(as.data.frame(ncounts_sig), file=(paste0("DESeq2_",date,"/tables/sigDEGs_norm_padj",padj.cutoff,"_lfc",lfc.cutoff,".csv")), row.names = FALSE)

# collaps replicates of the normalised count data
ncounts_sig_cps <- ncounts_sig

ncounts_sig_cps$WT_TM=rowMeans(ncounts_sig_cps[,c('WT_TM_1', 'WT_TM_2', 'WT_TM_3')])
ncounts_sig_cps$ssp_TM=rowMeans(ncounts_sig_cps[,c('ssp_TM_1', 'ssp_TM_2', 'ssp_TM_3')])
ncounts_sig_cps$ssp2_TM=rowMeans(ncounts_sig_cps[,c('ssp2_TM_1', 'ssp2_TM_2', 'ssp2_TM_3')])
ncounts_sig_cps$sspssp2_TM=rowMeans(ncounts_sig_cps[,c('sspssp2_TM_1', 'sspssp2_TM_2', 'sspssp2_TM_3')])

ncounts_sig_cps_tb <- ncounts_sig_cps %>%
  as_tibble() %>%
  dplyr::select("geneid", "WT_TM", "ssp_TM", "ssp2_TM", "sspssp2_TM" )
ncounts_sig_cps_tb$geneid <- gsub("\\..*","",ncounts_sig_cps_tb$geneid)




### Heatmap of significant DE genes

# # Set a color palette
# heat_colors <- brewer.pal(6, "YlOrRd")

# Run pheatmap using the metadata data frame for the annotation
xx<-pheatmap(ncounts_sig_cps_tb[2:5], # specify columns used for clustering; only use data columns and exclude geneid column, which results in NA error)
             color = viridis(10), 
             cluster_rows = TRUE, 
             cluster_cols = TRUE, 
             show_rownames = FALSE,
             cutree_rows = 1, # n number of kmeans cluster
             # annotation = meta, 
             border_color = NA, 
             fontsize = 10, 
             scale = "row", 
             fontsize_row = 10, 
             height = 20,
             main = paste0("n= ", nrow(ncounts_sig_cps_tb)," genes (padj=",padj.cutoff,")" ))

save_pheatmap_pdf <- function(x, filename, width=2, height=6) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(xx, (file.path(paste0("./DESeq2_",date,"/plots/heatmap_ncounts_padj",padj.cutoff,"_lfc",lfc.cutoff,".pdf"))))


