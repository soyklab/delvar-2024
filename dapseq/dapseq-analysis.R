################################################################################
# Analysis of DAP-seq results for Glaus et al., 2024
################################################################################

### Experiment description

##### - DAP-seq libraries:
###### - SSP-1:SSP as bait, replicate 1 
###### - SSP-2: SSP as bait, replicate 2
###### - SlycSSP2-1: S. lycopersicum SSP2 as bait, replicate 1
###### - SlycSSP2-2: S. lycopersicum SSP2 as bait, replicate 1
###### - SpimSSP2-1: S. pimpinellifolium SSP2 as bait, replicate 1
###### - SpimSSP2-2: S. pimpinellifolium SSP2 as bait, replicate 2

##### - Peaks were generated using csaw and prefiltered by padj (<=0.01)

################################################################################
# 0.1 Set environment

# set experiment name and date
experiment = "DAP-seq_bZIP" ### change name according to data
# set directories
setwd("./")
home <- getwd()
peakdir <- ("../peaks/") # directory which the csaw result table "peak_table_fdr_0.01.tsv"
degdir <-  ("../degs/")
countdir <- ("../counts/")
# load packages
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(eulerr))
suppressMessages(library(VennDiagram))
suppressMessages(library(pheatmap))
suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(RColorBrewer))
suppressMessages(library(viridis))
suppressMessages(library("gridExtra"))
suppressMessages(library("grid"))
suppressMessages(library("EnhancedVolcano"))
suppressMessages(library(reshape2))
suppressMessages(library(multcompView))
library("gridExtra")
library("grid")
# Create dated parent directory for file exports
date <- format(Sys.time(), "%Y%m%d")
dir.create(file.path(home,paste0(date,"_",experiment)), showWarnings = FALSE)
# Create output directories
dir.create(file.path(paste0(date, "_", experiment, "/peaks/")), showWarnings = FALSE)
dir.create(file.path(paste0(date, "_", experiment, "/genes/")), showWarnings = FALSE)
dir.create(file.path(paste0(date, "_", experiment, "/Venn/")), showWarnings = FALSE)
dir.create(file.path(paste0(date, "_", experiment, "/DAP-RNAseq_integration/")), showWarnings = FALSE)
dir.create(file.path(paste0(date, "_", experiment, "/plots/")), showWarnings = FALSE)
dir.create(file.path(home,paste0(date,"_",experiment, "/statistics")), showWarnings = FALSE)
dir.create(file.path(home, paste0(date, "_", experiment, "/ClusterProfiler")), showWarnings = FALSE)

################################################################################
# 0.2 Import data

# import gene annotation file (S100v2.0.0)
desc <- read_delim(file="/FILEPATHTO/SollycSweet-100_v2.0.fasta.descriptions", col_names = F) %>%
  dplyr::select(X1, X6) %>%
  dplyr::rename(geneid=X1, description=X6)
desc$geneid <- gsub("\\..*","",desc$geneid)

# import DAP-seq peak table
dat <- read_tsv(paste0(peakdir, "peak_table_fdr_0.01.tsv"), col_names = TRUE)
# Add peak IDs (Chr:start)
dat_tb <- dat %>%
  unite("peak_id", chr:start, remove = FALSE)
# Check data format
dat_tb
# Print total number of rows (features)
nrow(dat_tb)

# Write bed files for Megadepth / coverage analysis
dat.bed <- dat_tb %>%
  dplyr::select(chr, start, end)
write.table(dat.bed, (file.path(paste0(date, "_", experiment, "/peaks/peak_table_fdr_0.01.bed"))),row.names = F,col.names = F, sep="\t", quote=FALSE)
# this bed file is used to extract coverage using Megadepth on curnagl (at 49,170 peaks)

################################################################################



################################################################################
# 1. IDENTIFICATION OF SIGNIFICANT PEAKS AND ASSOCIATED GENES

# Set fold change (fs) cutoff
fc<-3
fc

#################################
# 1.1 SSP

# 1.1.1 SSP - peaks
# subset peaks by bait and fc cutoff
dat_SSP_tb <- subset(dat_tb, logFC.SSP  >= fc)
# Total number of peaks
nrow(dat_tb)
# Number of peaks >fc cutoff
nrow(dat_SSP_tb)
# Transform peaks into vector
SSP_peaks <- c(dat_SSP_tb$peak_id)
# Print number of peaks in vector
length(SSP_peaks)  #12358
# Write to file
write(SSP_peaks, (file.path(paste0(date, "_", experiment, "/peaks/SSP_peaks_fc",fc,".txt"))))

# 1.1.2 SSP - genes
# Split genes into separate columns
dat_SSP_tb <- separate(data = dat_SSP_tb, col = gene, into = c("gene_1", "gene_2", "gene_3", "gene_4"), sep = ",")
# Combine genes into vector and remove NAs and duplicated genes (genes with multiple peaks)
SSP_genes <- c(dat_SSP_tb$gene_1, dat_SSP_tb$gene_2, dat_SSP_tb$gene_3, dat_SSP_tb$gene_4)
SSP_genes <- unique(SSP_genes[!is.na(SSP_genes)])
SSP_genes <- gsub("\\..*","",SSP_genes)
# Print number of genes in vector
length(SSP_genes) #7688
# Write to file
write(SSP_genes, (file.path(paste0(date, "_", experiment, "/genes/SSP_genes_fc",fc,".txt"))))
# Add description and write to .csv file
SSP_genes_tb <- as.tibble(SSP_genes) %>%
  dplyr::rename(geneid=value) %>% # rename columns
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  distinct(geneid, .keep_all = TRUE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/genes/SSP_genes_fc",fc,".csv")))
length(SSP_genes)
nrow(SSP_genes_tb)

#################################
# 1.2 SlycSSP2

# 1.2.1 SlycSSP2 - peaks
# subset peaks by bait and fc cutoff
dat_SlycSSP2_tb <- subset(dat_tb, logFC.SlycSSP2  >= fc)
# Total number of peaks
nrow(dat_tb)
# Number of peaks >fc cutoff
nrow(dat_SlycSSP2_tb)
# Transform peaks into vector
SlycSSP2_peaks <- c(dat_SlycSSP2_tb$peak_id)
# Print number of peaks in vector
length(SlycSSP2_peaks) #2433
# Write to file
write(SlycSSP2_peaks, (file.path(paste0(date, "_", experiment, "/peaks/SlycSSP2_peaks_fc",fc,".txt"))))

# 1.2.2 SlycSSP2 - genes
# Split genes into separate columns
dat_SlycSSP2_tb <- separate(data = dat_SlycSSP2_tb, col = gene, into = c("gene_1", "gene_2", "gene_3", "gene_4"), sep = ",")
# Combine genes into vector and remove NAs and duplicated genes (genes with multiple peaks)
SlycSSP2_genes <- c(dat_SlycSSP2_tb$gene_1, dat_SlycSSP2_tb$gene_2, dat_SlycSSP2_tb$gene_3)
SlycSSP2_genes <- unique(SlycSSP2_genes[!is.na(SlycSSP2_genes)])
SlycSSP2_genes <- gsub("\\..*","",SlycSSP2_genes)
# Print number of genes in vector
length(SlycSSP2_genes) #1616
# Write to file
write(SlycSSP2_genes, (file.path(paste0(date, "_", experiment, "/genes/SlycSSP2_genes_fc",fc,".txt"))))
# Add description and write to .csv file
SlycSSP2_genes_tb <- as.tibble(SlycSSP2_genes) %>%
  dplyr::rename(geneid=value) %>% # rename columns
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  distinct(geneid, .keep_all = TRUE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/genes/SlycSSP2_genes_fc",fc,".csv")))
length(SlycSSP2_genes)
nrow(SlycSSP2_genes_tb)

#################################
# 1.3 SpimSSP2

# 1.3.1 SpimSSP2 - peaks
# subset peaks by bait and fc cutoff
dat_SpimSSP2_tb <- subset(dat_tb, logFC.SpimSSP2  >= fc)
# Total number of peaks
nrow(dat_tb)
# Number of peaks >fc cutoff
nrow(dat_SpimSSP2_tb) #8130
# Transform peaks into vector
SpimSSP2_peaks <- c(dat_SpimSSP2_tb$peak_id)
# Print number of peaks in vector
length(SpimSSP2_peaks)
# Write to .txt file
write(SpimSSP2_peaks, (file.path(paste0(date, "_", experiment, "/peaks/SpimSSP2_peaks_fc",fc,".txt"))))

# 1.3.2 SpimSSP2 - genes
# Split genes into separate columns
dat_SpimSSP2_tb <- separate(data = dat_SpimSSP2_tb, col = gene, into = c("gene_1", "gene_2", "gene_3", "gene_4"), sep = ",")
# Combine genes into vector and remove NAs and duplicated genes (genes with multiple peaks)
SpimSSP2_genes <- c(dat_SpimSSP2_tb$gene_1, dat_SpimSSP2_tb$gene_2, dat_SpimSSP2_tb$gene_3)
SpimSSP2_genes <- unique(SpimSSP2_genes[!is.na(SpimSSP2_genes)])
SpimSSP2_genes <- gsub("\\..*","",SpimSSP2_genes)
# Print number of genes in vector
length(SpimSSP2_genes) #5010
# Write to .txt file
write(SpimSSP2_genes, (file.path(paste0(date, "_", experiment, "/genes/SpimSSP2_genes_fc",fc,".txt"))))
# Add description and write to .csv file
SpimSSP2_genes_tb <- as.tibble(SpimSSP2_genes) %>%
  dplyr::rename(geneid=value) %>% # rename columns
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  distinct(geneid, .keep_all = TRUE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/genes/SpimSSP2_genes_fc",fc,".csv")))
length(SpimSSP2_genes)
nrow(SpimSSP2_genes_tb)



################################################################################
# 2. VISUALISATION OF SIGNIFICANT PEAKS AND ASSOCIATED GENES

# 2.1 Compare gene and peak sets in Venn-Euler diagrams to determine the intersection of peaks/genes bound by the different bait proteins. Genes of different sectors are written with gene descriptions into .csv files.

# 2.1.1 Gene comparisons
type="genes"
# Set data
SSP <- SSP_genes
SlycSSP2 <- SlycSSP2_genes
SpimSSP2 <- SpimSSP2_genes

# Calculate total number of unique features
length(SSP)
length(SlycSSP2)
length(SpimSSP2)
combined <- c(SSP,SlycSSP2,SpimSSP2)
length(combined)
length(unique(combined))
total <- length(unique(combined))

# Calculate overlap between vectors
overlap <- calculate.overlap(
  x = list(
    "SSP" = SSP,
    "SlycSSP2" = SlycSSP2,
    "SpimSSP2" = SpimSSP2
  )
);

# rename the overlap lists (for a 3-way overlap)
names(overlap) <- c("a123", "a12", "a13", "a23", "a1", "a2", "a3")

# calculate length of overlaps
a1=length(overlap$a1)
a2=length(overlap$a2)
a3=length(overlap$a3)
a12=length(overlap$a12)
a13=length(overlap$a13)
a23=length(overlap$a23)
a123=length(overlap$a123)

# write overlaps with gene annotation to file
overlap_a1 <- as_data_frame(overlap$a1) %>%
  rename(geneid=value) %>%
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SSP_fc",fc,".csv")))

overlap_a2 <- as_data_frame(overlap$a2) %>%
  rename(geneid=value) %>%
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SlycSSP2_fc",fc,".csv")))

overlap_a3 <- as_data_frame(overlap$a3) %>%
  rename(geneid=value) %>%
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SpimSSP2_fc",fc,".csv")))

overlap_a12 <- as_data_frame(overlap$a12) %>%
  rename(geneid=value) %>%
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SSP&SlycSSP2_fc",fc,".csv")))

overlap_a13 <- as_data_frame(overlap$a13) %>%
  rename(geneid=value) %>%
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SSP&SpimSSP2_fc",fc,".csv")))

overlap_a23 <- as_data_frame(overlap$a23) %>%
  rename(geneid=value) %>%
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SlycSSP2&SpimSSP2_fc",fc,".csv")))

overlap_a123 <- as_data_frame(overlap$a123) %>%
  rename(geneid=value) %>%
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SSP&SlycSSP2&SpimSSP2_fc",fc,".csv")))

# Make Euler object
euler <- euler(c("SSP" = a1,
                 "SlycSSP2" = a2,
                 "SpimSSP2" = a3,
                 "SSP&SlycSSP2" = a12,
                 "SSP&SpimSSP2" = a13,
                 "SlycSSP2&SpimSSP2" = a23,
                 "SSP&SlycSSP2&SpimSSP2" = a123))

# Select colors 
myCol <- c("#FCC427", "#35B779", "#31688E") # Custom for SSP-SlycSSP2-SpimSSP2

length(overlap$a1)
SSP_genes_unique <- overlap$a1
length(overlap$a2)
SlycSSP2_genes_unique <- overlap$a2
length(overlap$a3)
SpimSSP2_genes_unique <- overlap$a3
length(overlap$a13)
SSP_SpimSSP2_genes_shared <- overlap$a13
length(overlap$a13)
SSP_SpimSSP2_genes_shared <- overlap$a13
length(overlap$a123)
SSP_SlSpSSP2_genes_shared <- overlap$a123

# Plot Euler object
euler.plot <- plot(euler,
                   quantities = list(type = c("counts")),
                   fills = myCol,
                   alpha=0.5, # change transparency depending on the used color palette
                   edges=F,
                   main=paste0("DAP-seq ", type, " (fdr<0.01, fc>",fc,", total=", total,")" ))
euler.plot

# Save plot to file
# 1.Open a pdf file
pdf((file.path(paste0(date, "_", experiment, "/Venn/Venn_", type, "_fc", fc, ".pdf"))), width = 5, height = 4, useDingbats=FALSE)
# 2. Create a plot
euler.plot
# Close the pdf file
dev.off() 


# 2.1.2 Peak comparisons
type="peaks"
# Set data
SSP <- SSP_peaks
SlycSSP2 <- SlycSSP2_peaks
SpimSSP2 <- SpimSSP2_peaks

# Calculate total number of unique features
length(SSP)
length(SlycSSP2)
length(SpimSSP2)
combined <- c(SSP,SlycSSP2,SpimSSP2)
length(combined)
length(unique(combined))
total <- length(unique(combined))
# Calculate overlap between vectors
overlap <- calculate.overlap(
  x = list(
    "SSP" = SSP,
    "SlycSSP2" = SlycSSP2,
    "SpimSSP2" = SpimSSP2
  )
);

# rename the overlap lists (for a 3-way overlap)
names(overlap) <- c("a123", "a12", "a13", "a23", "a1", "a2", "a3")
# calculate length of overlaps
a1=length(overlap$a1)
a2=length(overlap$a2)
a3=length(overlap$a3)
a12=length(overlap$a12)
a13=length(overlap$a13)
a23=length(overlap$a23)
a123=length(overlap$a123)

# write overlaps with gene annotation to file
overlap_a1 <- as_data_frame(overlap$a1) %>%
  rename(peak_id=value) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SSP_fc",fc,".csv")))
overlap_a2 <- as_data_frame(overlap$a2) %>%
  rename(peak_id=value) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SlycSSP2_fc",fc,".csv")))
overlap_a3 <- as_data_frame(overlap$a3) %>%
  rename(peak_id=value) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SpimSSP2_fc",fc,".csv")))
overlap_a12 <- as_data_frame(overlap$a12) %>%
  rename(peak_id=value) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SSP&SlycSSP2_fc",fc,".csv")))
overlap_a13 <- as_data_frame(overlap$a13) %>%
  rename(peak_id=value) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SSP&SpimSSP2_fc",fc,".csv")))
overlap_a23 <- as_data_frame(overlap$a23) %>%
  rename(peak_id=value) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SlycSSP2&SpimSSP2_fc",fc,".csv")))
overlap_a123 <- as_data_frame(overlap$a123) %>%
  rename(peak_id=value) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/Venn/", type,"_overlap_SSP&SlycSSP2&SpimSSP2_fc",fc,".csv")))

# Make Euler object
euler <- euler(c("SSP" = a1,
                 "SlycSSP2" = a2,
                 "SpimSSP2" = a3,
                 "SSP&SlycSSP2" = a12,
                 "SSP&SpimSSP2" = a13,
                 "SlycSSP2&SpimSSP2" = a23,
                 "SSP&SlycSSP2&SpimSSP2" = a123))

# Select colors 
myCol <- c("#FCC427", "#35B779", "#31688E") # Custom for SSP-SlycSSP2-SpimSSP2

# Plot Euler object
euler.plot <- plot(euler,
                   quantities = list(type = c("counts")),
                   fills = myCol,
                   alpha=0.5, # change transparency depending on the used color palette
                   edges=F,
                   main=paste0("DAP-seq ", type, " (fdr<0.01, fc>",fc,", total=", total,")" )
)
euler.plot

# Save plot to file
# 1.Open a pdf file
pdf((file.path(paste0(date, "_", experiment, "/Venn/Venn_", type, "_fc", fc, ".pdf"))), width = 5, height = 4, useDingbats=FALSE)
# 2. Create a plot
euler.plot
# Close the pdf file
dev.off()
################################################################################


################################################################################
# 3. Integration of RNA-seq data

# Set DAP-seq data
dap_SSP <- SSP_genes
dap_SlycSSP2 <- SlycSSP2_genes
dap_SpimSSP2 <- SpimSSP2_genes

# Combine DAP-seq genes in one vector (define selection and datasets to combine)
combi <- "SSP_SlycSSP2_SpimSSP2" #all

# Keep unique (remove duplicates)
dapseq <- unique(c(dap_SSP, dap_SlycSSP2, dap_SpimSSP2)) #all
dapseq <- unique(c(dap_SSP, dap_SlycSSP2, dap_SpimSSP2)) #all

# Print number of genes in analysis
length(dapseq) # 8395


# Import RNAseq data

# Import log2FC and p-values (from DESeq2 analysis)
deg_ssp <- read_csv(paste0(degdir, "res_ssp_vs_WT.csv"), col_names = T)
deg_ssp$geneid <- gsub("\\..*","",deg_ssp$geneid)
deg_ssp2 <- read_csv(paste0(degdir, "res_ssp2_vs_WT.csv"), col_names = T)
deg_ssp2$geneid <- gsub("\\..*","",deg_ssp2$geneid)
deg_sspssp2 <- read_csv(paste0(degdir, "res_sspssp2_vs_WT.csv"), col_names = T)
deg_sspssp2$geneid <- gsub("\\..*","",deg_sspssp2$geneid)

# Import normalized counts (from DESeq2 analysis)
counts <- read_csv(paste0(countdir, "ncounts.cpm.csv"), col_names = T) %>%
  rowwise() %>% dplyr::mutate(WT = mean(c(WT_TM_1, WT_TM_2, WT_TM_3))) %>%
  rowwise() %>% dplyr::mutate(ssp = mean(c(ssp_TM_1, ssp_TM_2, ssp_TM_3))) %>%
  rowwise() %>% dplyr::mutate(ssp2 = mean(c(ssp2_TM_1, ssp2_TM_2, ssp2_TM_3))) %>%
  rowwise() %>% dplyr::mutate(sspssp2 = mean(c(sspssp2_TM_1, sspssp2_TM_2, sspssp2_TM_3)))
counts$geneid <- gsub("\\..*","",counts$geneid)
counts

####################################
### Set cutoff for DEGs
####################################
lfc=0.58
p=0.05
####################################

deg_ssp_lfc <- deg_ssp %>%
  subset(log2FoldChange >= lfc | log2FoldChange <= -lfc) %>%
  subset(padj <= p) %>%
  drop_na(padj)
nrow(deg_ssp_lfc)

deg_ssp2_lfc <- deg_ssp2 %>%
  subset(log2FoldChange >= lfc | log2FoldChange <= -lfc) %>%
  subset(padj <= p) %>%
  drop_na(padj)
nrow(deg_ssp2_lfc)

deg_sspssp2_lfc <- deg_sspssp2 %>%
  subset(log2FoldChange >= lfc | log2FoldChange <= -lfc) %>%
  subset(padj <= p) %>%
  drop_na(padj)

# Combine significant DEGs in one vector
rnaseq <- unique(c(deg_ssp_lfc$geneid, deg_ssp2_lfc$geneid, deg_sspssp2_lfc$geneid))
length(rnaseq) #1816

length(rnaseq)
length(deg_ssp_lfc$geneid)
length(deg_ssp2_lfc$geneid)
length(deg_sspssp2_lfc$geneid)

# Compare DAP-seq and RNA-seq datasets in Venn-Euler diagrams to determine intersection of genes bound by the bait protein(s) and the genes differentially expressed genes (DEGs) in the mutant(s)
# calculate overlap between vectors
overlap <- calculate.overlap(
  x = list(
    "deg" = rnaseq,
    "dap" = dapseq
  )
);

n.rnaseq <- length(rnaseq)
n.dapseq <- length(dapseq)
str(overlap)

length(overlap$a1) # 1816
length(overlap$a2) # 8395
length(overlap$a3) # 591
direct_deg <- overlap$a3 # make vector for ClusterProfiler

# rename the overlap lists
# for a 2-way overlap
names(overlap) <- c("a1", "a2", "a12")

# calculate length of overlaps
a1=length(overlap$a1)
a2=length(overlap$a2)
a12=length(overlap$a12)

# Make Euler object
euler <- euler(c("rnaseq" = a1,
                 "dapseq" = a2,
                 "rnaseq&dapseq" = a12))

# Select colors
myCol <- c("#7fc97f", "#beaed4", "#fdc086") #Option 2: select specific colors, see https://colorbrewer2.org

# plot euler
plot.euler <- plot(euler,
                   quantities = list(type = c("counts")),
                   fills = myCol,
                   alpha=0.5, # change transparency depending on the used color palette
                   edges=F,
                   main=paste0(combi, " (lfc>", lfc, "/ p<", p, "), n(rnaseq)=",n.rnaseq,", n(dapseq)=", n.dapseq)
)
plot.euler

# save plot to file euler
# 1.Open a pdf file
pdf((file.path(paste0(date, "_", experiment, "/DAP-RNAseq_integration/Venn_",combi, "_targets_dap_fc",fc,"_rna_lfc", lfc,"_padj", p,".pdf"))), width = 10, height = 8, useDingbats=FALSE)
# 2. Create a plot
plot.euler
# Close the pdf file
dev.off()

# SIMPLE VENN
# mycolors <- c("#7fc97f", "#beaed4", "#fdc086")
mycolors <- c("#7fc97f", "#beaed4")

# Generate sets
set1 <- rnaseq
set2 <- dapseq

# Chart
venn <- venn.diagram(
  x = list(set1, set2),
  category.names = c("rnaseq" , "dapseq"),
  filename = file.path(paste0(date, "_", experiment, "/DAP-RNAseq_integration/Venn_",combi, "_targets_dap_fc",fc,"_rna_lfc", lfc,"_padj", p,".svg")),
  output=FALSE,
  
  # Output features
  imagetype="svg" ,
  height = 4 ,
  width = 4 ,
  # resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = mycolors
)

# Plot overlap between datasets in heatmaps (normalized counts for all differentially expressed DAP-Seq genes) to determine potential direct targets with similar expression patterns between mutant(s)

# Filter normalized counts for genes that are differentially expressed and contain a DAP-seq binding peak
counts_deg_dap <- unique(counts) %>%
  filter(geneid %in% dapseq) %>%
  filter(geneid %in% rnaseq) %>%
  select(geneid, WT, ssp, ssp2, sspssp2) %>%
  column_to_rownames(var = "geneid")
nrow(counts_deg_dap) # 591
length(rnaseq) # 1816
length(dapseq) # 8395

total <- nrow(counts_deg_dap)
nrow(counts)
counts$geneid

# Run pheatmap
pheatmap.plot <- pheatmap(counts_deg_dap, # specify columns used for clustering; only use data columns and exclude geneid column, which results in NA error)
                          color = viridis(100), #inferno / viridis / plasma
                          cluster_rows = T,
                          cluster_cols = F,
                          show_rownames = T,
                          cutree_rows = 2, # set n number of kmeans cluster
                          border_color = NA,
                          fontsize = 6,
                          scale = "row",
                          fontsize_row = 2,
                          fontsize_col = 6,
                          height = 80,
                          treeheight_col = 6,
                          treeheight_row = 20,
                          main = paste0("DAP-seq targets fc=", fc,"\n",combi, "\n (logFC>=", lfc, "/ padj<=", p, ")\n total=", total)
)

# Save pheatmap (change plot height according to number of genes)
save_pheatmap_pdf <- function(x, filename, width=2, height=12) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(pheatmap.plot, (file.path(paste0(date, "_", experiment, "/DAP-RNAseq_integration/pheatmap_",combi, "_targets_dap_fc",fc,"_rna_lfc", lfc,"_padj", p,".pdf"))))

# Export heatmaps clusters: Genes from individual clusters in the heatmap are exported together with their gene description into a .csv file. This file can be searched for genes of special interest.

# Export pheatmap data and combine with logFC, padj, and gene description
pheat.cluster <- counts_deg_dap %>%
  cbind(cluster = cutree(pheatmap.plot$tree_row, k=2)) %>% # change number of cluster according to heatmap
  rownames_to_column(var = "geneid") %>%
  as_tibble %>%
  dplyr::rename(cpm_ssp = ssp, cpm_ssp2 = ssp2, cpm_sspssp2 = sspssp2) %>%
  left_join(deg_ssp) %>%
  dplyr::rename(ssp_log2FoldChange = log2FoldChange, ssp_padj = padj) %>%
  select(!c(baseMean, lfcSE, stat, pvalue, description)) %>%
  left_join(deg_ssp2) %>%
  dplyr::rename(ssp2_log2FoldChange = log2FoldChange, ssp2_padj = padj) %>%
  select(!c(baseMean, lfcSE, stat, pvalue, description)) %>%
  left_join(deg_sspssp2) %>%
  dplyr::rename(sspssp2_log2FoldChange = log2FoldChange, sspssp2_padj = padj) %>%
  select(!c(baseMean, lfcSE, stat, pvalue))

# Write to .csv
write.csv(as.data.frame(pheat.cluster), file=(paste0(date, "_", experiment, "/DAP-RNAseq_integration/pheatmap_",combi, "_targets_dap_fc",fc,"_rna_lfc", lfc,"_padj", p,"_2.csv")), row.names = FALSE)


# Lineplots of selected, putative direct targets

# define order of genotypes in vector
level_order <- c('WT', 'ssp', 'ssp2', "sspssp2")

# Select gene list
# goi <- c("Solyc03g007230", "Solyc06g051940", "Solyc05g052980") # PP2C
# goi <- c("Solyc06g050500", "Solyc03g095780", "Solyc10g076410") # PYL
# goi <- c("Solyc02g065730", "Solyc12g087830", "Solyc12g056460") # MADS
# goi <- c("Solyc04g071990", "Solyc12g056650", "Solyc01g005300") # GIGANTEA & FKF1
# goi <- c("Solyc04g016430", "Solyc04g080820", "Solyc10g084150") # CKX & LOG
# goi <- c("Solyc05g014260", "Solyc03g115770", "Solyc03g081240") # ARR
goi <- c("Solyc06g069430", "Solyc03g114830") #FUL1/2

plot.tpm <- counts %>%
  select(-WT, -ssp, -ssp2, -sspssp2) %>% # remove columns with mean values
  filter(geneid %in% goi) %>% # filter for GOI
  pivot_longer(!geneid, names_to = "sample", values_to = "counts") %>% # transpose
  mutate(genotype = gsub("(.*)_.*_.*..*", "\\1", sample)) %>% # by genotype (string before first underscore in filename)
  mutate(tag = str_remove(geneid, "Solyc")) %>% # shorten geneid to genetag
  ggplot(aes(x = factor(genotype, level = level_order), y = counts)) +
  geom_point(size = 0.5, alpha = 1, position = position_jitter(0.1, seed = 666)) +
  stat_summary(fun = mean, geom = "point", aes(group = tag, color = genotype), size = 1 ) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(group = tag, color = genotype)) +
  facet_wrap(~tag, scales = "free_y", ncol=3) +
  theme_minimal(base_size = 6) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 8)) +
  theme(legend.position = "none") +
  labs(x = "", y = "TPM") +
  scale_colour_brewer(palette = "Set1")
plot.tpm

# change file name based on gene list
ggsave((paste0(date, "_", experiment, "/DAP-RNAseq_integration/tmp_FUL12_tall.pdf")), height = 1.4, width = 1.5, bg = "white")


################################################################################
# 4. Identification of private peaks and associated genes

# 4.1 Extract differential binding data (SpimSSP2_vs_SlycSSP2) from peak table and add peak IDs (Chr:start)

# Make list of peaks that are bound by Slyc/Spim vs input by fc>3
length(SlycSSP2_peaks) #2433
length(SpimSSP2_peaks) #8130
SlycSSP2_SlycSSP2_peaks <- unique(c(SlycSSP2_peaks, SpimSSP2_peaks))
length(SlycSSP2_SlycSSP2_peaks) #9199

res_SpimSSP2_vs_SlycSSP2_tb <- dat %>%
  unite("peak_id", chr:start, remove = FALSE) %>%
  filter(peak_id %in% SlycSSP2_SlycSSP2_peaks) %>%
  # filter(T_flag_SpimSSP2_vs_SlycSSP2 == 1) %>% # set by Giovanna, all peaks that are found in either Slyc or SpimSSP2 with FDR<0.01
  dplyr::select(peak_id, chr, start, end, logFC.SpimSSP2_vs_SlycSSP2, PValue.SpimSSP2_vs_SlycSSP2, gene)
nrow(res_SpimSSP2_vs_SlycSSP2_tb) # 9199

# Set fold change (fs) cutoff
fc<-3
fc

# Filter by fc cutoff (both positive (private SpimSSP2 binding) and negative (private SlycSSP2 binding)) and write filtered peaks and genes into .txt files
sig_SpimSSP2_vs_SlycSSP2_tb <- subset(res_SpimSSP2_vs_SlycSSP2_tb,
                                      logFC.SpimSSP2_vs_SlycSSP2 >= fc | 
                                        logFC.SpimSSP2_vs_SlycSSP2 <= -fc) 
sig_SpimSSP2_vs_SlycSSP2_up_tb <- subset(res_SpimSSP2_vs_SlycSSP2_tb,
                                         logFC.SpimSSP2_vs_SlycSSP2 >= fc) # "up" peaks = private SpimSSP2
sig_SpimSSP2_vs_SlycSSP2_dn_tb <- subset(res_SpimSSP2_vs_SlycSSP2_tb,
                                         logFC.SpimSSP2_vs_SlycSSP2 <= -fc)  # "down" peaks = private SlycSSP2

# Number of peaks >fc cutoff
nrow(sig_SpimSSP2_vs_SlycSSP2_tb) #5063
nrow(sig_SpimSSP2_vs_SlycSSP2_up_tb) #4637
nrow(sig_SpimSSP2_vs_SlycSSP2_dn_tb) #426

# Write to file
write.csv(as.data.frame(sig_SpimSSP2_vs_SlycSSP2_tb), file=(paste0(date, "_", experiment, "/peaks/SpimSSP2_vs_SlycSSP2_peaks_diff_fc",fc,".csv")), row.names = FALSE)
write.csv(as.data.frame(sig_SpimSSP2_vs_SlycSSP2_up_tb), file=(paste0(date, "_", experiment, "/peaks/SpimSSP2_peaks_private_fc",fc,".csv")), row.names = FALSE)
write.csv(as.data.frame(sig_SpimSSP2_vs_SlycSSP2_dn_tb), file=(paste0(date, "_", experiment, "/peaks/SlycSSP2_peaks_private_fc",fc,".csv")), row.names = FALSE)


# 4.2 plot data

# 4.2.1 Plot differentially-bound peaks in Volcano Plots

# str(dat_SpimSSP2_vs_SlycSSP2_tb)
str(sig_SpimSSP2_vs_SlycSSP2_tb)

p1 <- EnhancedVolcano(res_SpimSSP2_vs_SlycSSP2_tb,
                      # lab = rownames(sig_SpimSSP2_vs_SlycSSP2_tb),
                      lab = NA,
                      title = 'SpimSSP2 vs SlycSSP2',
                      # subtitle = paste0("lfc.cutoff = ", lfc.cutoff, "; padj.cutoff = ", padj.cutoff ),
                      subtitle = paste0("logfc > ", fc),
                      caption = paste0("total = ", nrow(res_SpimSSP2_vs_SlycSSP2_tb), "; sig = ", nrow(sig_SpimSSP2_vs_SlycSSP2_tb), "\n down = ", nrow(sig_SpimSSP2_vs_SlycSSP2_dn_tb), "; up = ", nrow(sig_SpimSSP2_vs_SlycSSP2_up_tb)),
                      x = "logFC.SpimSSP2_vs_SlycSSP2",
                      y = "PValue.SpimSSP2_vs_SlycSSP2",
                      col = c('grey', 'green1', 'grey', 'purple1'),
                      colAlpha = 1,
                      # selectLab = genes,
                      drawConnectors = TRUE,
                      FCcutoff = fc,
                      pCutoff = 0.01
)
p1

pdf((file.path(paste0(date, "_", experiment, "/plots/volcano_lfc",fc,".pdf"))), width = 3, height = 5, useDingbats=FALSE)
p1
dev.off() # Close the file


# 4.2.2 Plot Venn diagrams between private SlycSSP2 peaks from DB and Venn analysis
ssp2_priv_venn.df <- read_csv(file.path(paste0(date, "_", experiment, "/Venn/peaks_overlap_SlycSSP2_fc3.csv")))
nrow(ssp2_priv_venn.df) #991
ssp2_priv_venn <- ssp2_priv_venn.df$peak_id
length(ssp2_priv_venn) #991

ssp2_priv_db.df <- read_csv(file.path(paste0(date, "_", experiment, "/peaks/SlycSSP2_peaks_private_fc3.csv")))
nrow(ssp2_priv_db.df) #527
ssp2_priv_db <- ssp2_priv_db.df$peak_id
length(ssp2_priv_db) #527

# mycolors <- c("#7fc97f", "#beaed4", "#fdc086")
mycolors <- c("#7fc97f", "#beaed4")

# Generate sets
set1 <- ssp2_priv_venn
set2 <- ssp2_priv_db

# Chart svg
venn <- venn.diagram(
  x = list(set1, set2),
  category.names = c("Venn" , "DB"),
  filename = file.path(paste0(date, "_", experiment, "/Venn/SlycSSP2_priv_peaks_Venn_vs_DB.svg")),
  output=TRUE,
  
  # Output features
  imagetype="svg" ,
  height = 4 ,
  width = 4 ,
  # resolution = 300,
  # compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = mycolors
)


# 4.3 Write bed files

# 4.3.1  for MEME / motif analysis
SpimSSP2_private_peaks.bed <- sig_SpimSSP2_vs_SlycSSP2_up_tb %>% 
  select(chr, start, end)
write.table(SpimSSP2_private_peaks.bed, (file.path(paste0(date, "_", experiment, "/peaks/SpimSSP2_peaks_private_fc",fc,".bed"))),row.names = F,col.names = F, sep="\t", quote=FALSE)

SlycSSP2_private_peaks.bed <- sig_SpimSSP2_vs_SlycSSP2_dn_tb %>% 
  select(chr, start, end)
write.table(SlycSSP2_private_peaks.bed, (file.path(paste0(date, "_", experiment, "/peaks/SlycSSP2_peaks_private_fc",fc,".bed"))),row.names = F,col.names = F, sep="\t", quote=FALSE)

# select top 1000 peaks (with lowest p-value) (optimal number for MEME to use) (DO NOT do this for lists with less than 1000)
SpimSSP2_private_top1000 <- sig_SpimSSP2_vs_SlycSSP2_up_tb %>%
  slice_min(PValue.SpimSSP2_vs_SlycSSP2, n = 1000)

# write bed files for 200 bp intervals around the peak center
SlycSSP2_private_peaks_200bp.bed <- sig_SpimSSP2_vs_SlycSSP2_dn_tb %>%
  mutate(length = end-start) %>%
  mutate(center = round(start+length/2)) %>%
  mutate(start_200 = center-100) %>%
  mutate(end_200 = center+100) %>%
  select(chr, start_200, end_200)
write.table(SlycSSP2_private_peaks_200bp.bed, (file.path(paste0(date, "_", experiment, "/peaks/SlycSSP2_peaks_private_fc",fc,"_200bp.bed"))),row.names = F,col.names = F, sep="\t", quote=FALSE)
# this bed file is used for motif enrichment analysis using MEME (SlycSSP2, n=426)

SpimSSP2_private_peaks_200bp.bed <- SpimSSP2_private_top1000 %>%
  mutate(length = end-start) %>%
  mutate(center = round(start+length/2)) %>%
  mutate(start_200 = center-100) %>%
  mutate(end_200 = center+100) %>%
  select(chr, start_200, end_200)
write.table(SpimSSP2_private_peaks_200bp.bed, (file.path(paste0(date, "_", experiment, "/peaks/SpimSSP2_peaks_private_fc",fc,"_200bp.bed"))),row.names = F,col.names = F, sep="\t", quote=FALSE)
# this bed file is used for motif enrichment analysis using MEME  (SpimSSP2, n=1000)


# 4.4 Plot normalized read coverage at private peaks

# load coverage data (generated using bed files above and bigwig files with Megadepth)
SlycSSP2_1 <- read_delim(file="../bed/SlycSSP2_1_cov.bed", col_names = F) %>%
  dplyr::rename(chr=X1, start=X2, end=X3, SlycSSP2_1_cov=X4)
SlycSSP2_2 <- read_delim(file="../bed/SlycSSP2_2_cov.bed", col_names = F) %>%
  dplyr::rename(chr=X1, start=X2, end=X3, SlycSSP2_2_cov=X4) %>%
  select(SlycSSP2_2_cov)
SpimSSP2_1 <- read_delim(file="../bed/SpimSSP2_1_cov.bed", col_names = F) %>%
  dplyr::rename(chr=X1, start=X2, end=X3, SpimSSP2_1_cov=X4) %>%
  select(SpimSSP2_1_cov)
SpimSSP2_2 <- read_delim(file="../bed/SpimSSP2_2_cov.bed", col_names = F) %>%
  dplyr::rename(chr=X1, start=X2, end=X3, SpimSSP2_2_cov=X4) %>%
  select(SpimSSP2_2_cov)

# bind data
dat_cov <- SlycSSP2_1 %>% 
  bind_cols(SlycSSP2_2, SpimSSP2_1, SpimSSP2_2) %>%
  unite("peak_id", chr:start, remove = FALSE) # add peak id

# calculate means
dat_cov <- dat_cov %>%
  dplyr::mutate(SlycSSP2_mean = rowMeans(select(dat_cov, starts_with("SlycSSP2")), na.rm = TRUE)) %>% # calculate coverage mean
  dplyr::mutate(SpimSSP2_mean = rowMeans(select(dat_cov, starts_with("SpimSSP2")), na.rm = TRUE)) # calculate coverage mean
nrow(dat_cov) # 8311
dat_cov

# set private peaks
SlycSSP2_private_peaks <- sig_SpimSSP2_vs_SlycSSP2_dn_tb %>%
  mutate(SlycSSP2_private_peaks = if_else(!is.na(peak_id), "SlycSSP2_private", "NA")) %>%
  dplyr::select(peak_id, SlycSSP2_private_peaks)
nrow(SlycSSP2_private_peaks)

SpimSSP2_private_peaks <- sig_SpimSSP2_vs_SlycSSP2_up_tb %>%
  mutate(SpimSSP2_private_peaks = if_else(!is.na(peak_id), "SpimSSP2_private", "NA")) %>%
  dplyr::select(peak_id, SpimSSP2_private_peaks)
nrow(SpimSSP2_private_peaks)

# join coverage data with private peaks
dat_cov_sig_private <- dat_cov %>%
  left_join(SlycSSP2_private_peaks, by=NULL, suffix = c("peak_id", "peak_id"), copy = FALSE,  keep = FALSE) %>% # join private peak data
  left_join(SpimSSP2_private_peaks, by=NULL, suffix = c("peak_id", "peak_id"), copy = FALSE,  keep = FALSE) # join private peak data
nrow(dat_cov_sig_private) # 8311
dat_cov_sig_private

# assign private peaks
dat_cov_sig_private_sig <- mutate(dat_cov_sig_private, priv_class = case_when(
  SlycSSP2_private_peaks == "SlycSSP2_private" ~ "SlycSSP2_private",
  SpimSSP2_private_peaks == "SpimSSP2_private" ~ "SpimSSP2_private",
  TRUE ~ NA_character_ ))
dat_cov_sig_private_sig_plot <- dat_cov_sig_private_sig %>%
  filter(peak_id %in% SlycSSP2_SlycSSP2_peaks) %>% # plot all SSP2 peaks (n=9199)
  dplyr::select(peak_id, SlycSSP2_mean, SpimSSP2_mean, priv_class) %>%
  replace_na(list(priv_class = "common")) # rename non-private peaks
nrow(dat_cov_sig_private) #49170 (all peaks in dataset peak file)
nrow(dat_cov_sig_private_sig) #49170  (all peaks in dataset peak file)
nrow(dat_cov_sig_private_sig_plot) #9199 (all SSP2 peaks only)

# plots
ggplot(dat_cov_sig_private_sig_plot, aes(x=SlycSSP2_mean, y=SpimSSP2_mean, color=priv_class)) +
  geom_point(alpha=1, size=1) +
  scale_x_log10(limits=c(1e-01, 1e+07)) +
  scale_y_log10(limits=c(1e-01, 1e+07)) +
  geom_abline(slope = 1.0,
              intercept = 0,
              color="red",
              linetype = "dashed") +
  scale_color_manual(values=c("#BEBEBE", "#35B779", "#31688E")) +
  theme_light() + 
  ggtitle(paste0("peak coverage with fc fc>", fc, ";n(Slyc)=", nrow(SlycSSP2_private_peaks), " ;n(Spim)=", nrow(SpimSSP2_private_peaks)))

ggsave((file.path(paste0(date, "_", experiment, "/plots/private_cov_distr_fc",fc,".pdf"))), width = 4, height = 2.5, useDingbats=FALSE)


# 4.4.2 plot read coverage distribution in violin plots

# melt data frame from wide to long format
dat_plot <- melt(dat_cov_sig_private_sig_plot, id=c('peak_id', "priv_class"))
dat_plot
#rename columns
names(dat_plot) <- c('peak_id', 'priv_class', 'sample',"mean_cov")
head(dat_plot)

# plot as violin
ggplot(dat_plot, aes(x=sample, y=mean_cov, fill = sample)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_log10(limits=c(1e-01, 1e+07)) +
  ggtitle(paste0("norm coverage with fc>", fc, ";n(Slyc)=", nrow(SlycSSP2_private_peaks), " ;n(Spim)=", nrow(SpimSSP2_private_peaks))) +
  facet_wrap(~priv_class) +
  scale_fill_manual(values=c("#35B779", "#31688E")) +
  theme(legend.position="none")
ggsave((file.path(paste0(date, "_", experiment, "/plots/private_cov_fc",fc,".pdf"))), width = 5, height = 4.5, useDingbats=FALSE)


# 4.4.3 Statistics for violin plot

# load scripts
source(file.path(scriptdir,"/FILEPATHTO/summarySE.R"))

# bind data
dat_sig <- dat_plot %>% 
  unite("sample_class", priv_class:sample, remove = FALSE) # add sample_class
head(dat_sig)
tail(dat_sig)
str(dat_sig)

# calculate summarySE
dat_sum <- summarySE(dat_sig, measurevar="mean_cov", groupvars=c("sample", "priv_class"), na.rm=TRUE)
head(dat_sum) ### prints head of the dataset to console 
# write summarySE to file
write.csv(dat_sum, file = (paste0(date, "_", experiment, "/statistics/summarySE_CPM.csv")))

# ANOVA and Tukey post-hoc on normalized read distribution
model=lm( dat_sig$mean_cov ~ dat_sig$sample_class )
ANOVA=aov(model)
# Tukey test to study each pair of treatment :
TUKEY <- TukeyHSD(x=ANOVA, 'dat_sig$sample_class', conf.level=0.95)
# Tuckey test representation :
plot(TUKEY , las=1 , col="brown" )
# I need to group the treatments that are not different each other together.
generate_label_df <- function(TUKEY, variable){
  # Extract labels and factor levels from Tukey post-hoc
  Tukey.levels <- TUKEY[[variable]][,4]
  Tukey.labels <- data.frame(multcompLetters(Tukey.levels)['Letters'])
  #I need to put the labels in the same order as in the boxplot :
  Tukey.labels$sample_class=rownames(Tukey.labels)
  Tukey.labels=Tukey.labels[order(Tukey.labels$sample_class) , ]
  return(Tukey.labels)
}

# Apply the function on my dataset
GROUPS=generate_label_df(TUKEY , "dat_sig$sample_class")
colnames(GROUPS) <- c("mean_cov", "sample_class")
GROUPS
# Extract results
TK_data<-as.data.frame(TUKEY[1])
# Write to file
write.csv(TK_data, file = (paste0(date, "_", experiment, "/statistics/Tukey_CPM_TUKEY.csv")))
write.csv(GROUPS, file = (paste0(date, "_", experiment, "/statistics/Tukey_CPM_GROUPS.csv")))


# TTEST on coverage
dat_cov_sig_private_sig_plot

dat_cov_SlycSSP2_private <- dat_cov_sig_private_sig_plot %>%
  select(peak_id, SlycSSP2_mean, priv_class) %>%
  filter(priv_class == "SlycSSP2_private") %>%
  dplyr::rename(CPM=SlycSSP2_mean)
nrow(dat_cov_SlycSSP2_private)
dat_cov_SpimSSP2_private <- dat_cov_sig_private_sig_plot %>%
  select(peak_id, SpimSSP2_mean, priv_class) %>%
  filter(priv_class == "SpimSSP2_private") %>%
  dplyr::rename(CPM=SpimSSP2_mean)
nrow(dat_cov_SpimSSP2_private)

# Compute t-test from df
# bind data
dat_cov_private <- dat_cov_SlycSSP2_private %>%
  bind_rows(dat_cov_SpimSSP2_private, .id=NULL)
nrow(dat_cov_private)

# conduct F-test to determine if variances are equal (or not)
res.ftest <- var.test(CPM ~ priv_class, data = dat_cov_private)
res.ftest
# conduct T-test
res <- t.test(CPM ~ priv_class, data = dat_cov_private, var.equal = FALSE)
res
# print the p-value
res$p.value
# print the mean
res$estimate
# print the confidence interval
res$conf.int
# write to file
write.csv(res$p.value, file = (paste0(date, "_", experiment, "/statistics/TTEST_CPM_pval.csv")))


# 4.5 Extract genes associated with private peaks

# Number of peaks >fc cutoff
nrow(sig_SpimSSP2_vs_SlycSSP2_up_tb) #SpimSSP2 private ("up") #7784
nrow(sig_SpimSSP2_vs_SlycSSP2_dn_tb) #SlycSSP2 private ("down") #527

# Split genes into separate columns
SpimSSP2_private_tb <- separate(data = sig_SpimSSP2_vs_SlycSSP2_up_tb, col = gene, into = c("gene_1", "gene_2", "gene_3", "gene_4"), sep = ",")
# Combine genes into vector and remove NAs and duplicated genes (genes with multiple peaks)
SpimSSP2_private_genes <- c(SpimSSP2_private_tb$gene_1, SpimSSP2_private_tb$gene_2, SpimSSP2_private_tb$gene_3, SpimSSP2_private_tb$gene_4)
SpimSSP2_private_genes <- unique(SpimSSP2_private_genes[!is.na(SpimSSP2_private_genes)])
SpimSSP2_private_genes <- gsub("\\..*","",SpimSSP2_private_genes)
# Print number of genes in vector
length(SpimSSP2_private_genes) # 3046
# Write to txt file
write(SpimSSP2_private_genes, (file.path(paste0(date, "_", experiment, "/genes/SpimSSP2_genes_private_fc",fc,".txt"))))
# Add description and write to .csv file
SpimSSP2_private_genes_tb <- as.tibble(SpimSSP2_private_genes) %>%
  dplyr::rename(geneid=value) %>% # rename columns
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  distinct(geneid, .keep_all = TRUE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/genes/SpimSSP2_genes_private_fc",fc,".csv")))
length(SpimSSP2_private_genes) #3046
nrow(SpimSSP2_private_genes_tb) #3046

# Split genes into separate columns
SlycSSP2_private_tb <- separate(data = sig_SpimSSP2_vs_SlycSSP2_dn_tb, col = gene, into = c("gene_1", "gene_2", "gene_3", "gene_4"), sep = ",")
# Combine genes into vector and remove NAs and duplicated genes (genes with multiple peaks)
SlycSSP2_private_genes <- c(SlycSSP2_private_tb$gene_1, SlycSSP2_private_tb$gene_2, SlycSSP2_private_tb$gene_3, SlycSSP2_private_tb$gene_4)
SlycSSP2_private_genes <- unique(SlycSSP2_private_genes[!is.na(SlycSSP2_private_genes)])
SlycSSP2_private_genes <- gsub("\\..*","",SlycSSP2_private_genes)
# Print number of genes in vector
length(SlycSSP2_private_genes) # 230
# Write to .txt file
write(SlycSSP2_private_genes, (file.path(paste0(date, "_", experiment, "/genes/SlycSSP2_genes_private_fc",fc,".txt"))))
# Add description and write to .csv file
SlycSSP2_private_genes_tb <- as.tibble(SlycSSP2_private_genes) %>%
  dplyr::rename(geneid=value) %>% # rename columns
  left_join(desc, by=NULL, suffix = c("geneid", "geneid"), copy = FALSE,  keep = FALSE) %>%
  distinct(geneid, .keep_all = TRUE) %>%
  write_csv(file.path(paste0(date, "_", experiment, "/genes/SlycSSP2_genes_private_fc",fc,".csv")))
length(SlycSSP2_private_genes) # 230
nrow(SlycSSP2_private_genes_tb) # 230




# 4.6 Plot Venn ssp2CR DEGs and SSP2F targets

ssp2CR_deg_genes <- unique(deg_ssp2_lfc$geneid)
length(ssp2CR_deg_genes) #180
length(SlycSSP2_private_genes) #230

# mycolors <- c("#7fc97f", "#beaed4", "#fdc086")
mycolors <- c("#7fc97f", "#beaed4")

# Generate sets
set1 <- ssp2CR_deg_genes
set2 <- SlycSSP2_private_genes

# Chart
venn <- venn.diagram(
  x = list(set1, set2),
  category.names = c("ssp2CR_deg_genes" , "SlycSSP2_private_genes"),
  filename = file.path(paste0(date, "_", experiment, "/Venn/venn_dap_SlycSSP2-deg_ssp2CR.svg")),
  output=FALSE,
  
  # Output features
  imagetype="svg" ,
  height = 4 ,
  width = 4 ,
  # resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = mycolors
)

# Compare private DAP-seq and RNA-seq datasets in Venn-Euler diagrams to determine intersection of genes bound by the bait protein(s) and the genes differentially expressed genes (DEGs) in the mutant(s)
# calculate overlap between vectors
overlap <- calculate.overlap(
  x = list(
    "deg" = ssp2CR_deg_genes,
    "dap" = SlycSSP2_private_genes
  )
);

n.ssp2CR_deg_genes <- length(ssp2CR_deg_genes)
n.SlycSSP2_private_genes <- length(SlycSSP2_private_genes)
str(overlap)

length(overlap$a1) # 180
length(overlap$a2) # 230
length(overlap$a3) # 1 (common to both sets)

# rename the overlap lists
# for a 2-way overlap
names(overlap) <- c("a1", "a2", "a12")

# calculate length of overlaps
a1=length(overlap$a1)
a2=length(overlap$a2)
a12=length(overlap$a12)

# Make Euler object
euler <- euler(c("rnaseq" = a1,
                 "dapseq" = a2,
                 "rnaseq&dapseq" = a12))

# Select colors
myCol <- c("#7fc97f", "#beaed4", "#fdc086") #Option 2: select specific colors, see https://colorbrewer2.org

# plot euler
plot.euler <- plot(euler,
                   quantities = list(type = c("counts")),
                   fills = myCol,
                   alpha=0.5, # change transparency depending on the used color palette
                   edges=F,
                   main=paste0(combi, " (lfc>", lfc, "/ p<", p, "), n(rnaseq)=",n.rnaseq,", n(dapseq)=", n.dapseq)
)
plot.euler

# save plot to file euler
# 1.Open a pdf file
pdf((file.path(paste0(date, "_", experiment, "/DAP-RNAseq_integration/Venn_",combi, "_private_targets_dap_fc",fc,"_rna_lfc", lfc,"_padj", p,".pdf"))), width = 10, height = 8, useDingbats=FALSE)
# 2. Create a plot
plot.euler
# Close the pdf file
dev.off()




# Plot the expression of the overlap (n=1 dehydrogenase genes)

str(overlap)

goi <- c("Solyc11g071290") #Solyc11g071290 (SSP2F169 private)

plot.tpm <- counts %>%
  select(-WT, -ssp, -ssp2, -sspssp2) %>% # remove columns with mean values
  filter(geneid %in% goi) %>% # filter for GOI
  pivot_longer(!geneid, names_to = "sample", values_to = "counts") %>% # transpose
  mutate(genotype = gsub("(.*)_.*_.*..*", "\\1", sample)) %>% # by genotype (string before first underscore in filename)
  mutate(tag = str_remove(geneid, "Solyc")) %>% # shorten geneid to genetag
  ggplot(aes(x = factor(genotype, level = level_order), y = counts)) +
  # geom_point(aes(fill = genotype), color = "white", size = 2,
  #            alpha = 0.8, shape = 21, position = position_jitter(0.1, seed = 666)) +
  geom_point(size = 0.5, alpha = 1, position = position_jitter(0.1, seed = 666)) +
  stat_summary(fun = mean, geom = "point", aes(group = tag, color = genotype), size = 1 ) +
  stat_summary(fun.data = mean_se, geom = "errorbar", aes(group = tag, color = genotype)) +
  facet_wrap(~tag, scales = "free_y", ncol=3) +
  theme_minimal(base_size = 6) +
  # theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), text = element_text(size = 8)) +
  # theme(legend.position = "top", text = element_text(size = 6)) +
  theme(legend.position = "none") +
  labs(x = "", y = "TPM") +
  scale_colour_brewer(palette = "Set1")
plot.tpm

# change file name based on gene list
ggsave((file.path(paste0(date, "_", experiment, "/DAP-RNAseq_integration/tpm_Solyc11g071290_tall.pdf"))), height = 1.4, width = 1.1, bg = "white")

################################################################################



################################################################################
# 5. GO-term analysis

# Load custom org.eg.db library (previously constructed using AnnotationForge)
library(org.SlycopersicumSL40VIBS.eg.db)

# Individual gene lists FC genes

# all targets
SSP2F_targets <- SlycSSP2_genes
length(SSP2F_targets) #1616

SSP2S_targets <- SpimSSP2_genes
length(SSP2S_targets) #5010

# top targets
SSP2F_top_targets_tb <- dat_tb %>%
  subset(logFC.SlycSSP2  >= fc) %>% # filter by fc
  filter(gene != "NA") %>% # filter for genic peaks
  slice_min(PValue.SlycSSP2, n = 1000) %>% # slice top peaks by pvalue
  separate(col = gene, into = c("gene_1", "gene_2", "gene_3", "gene_4"), sep = ",")
SSP2F_top_targets <- c(SSP2F_top_targets_tb$gene_1, SSP2F_top_targets_tb$gene_2, SSP2F_top_targets_tb$gene_3, SSP2F_top_targets_tb$gene_4)
SSP2F_top_targets <- unique(SSP2F_top_targets[!is.na(SSP2F_top_targets)])
SSP2F_top_targets <- gsub("\\..*","",SSP2F_top_targets)  
length(SSP2F_top_targets) #1152 for top 1000 peaks

SSP2S_top_targets_tb <- dat_tb %>%
  subset(logFC.SpimSSP2  >= fc) %>% # filter by fc
  filter(gene != "NA") %>% # filter for genic peaks nrow
  slice_min(PValue.SpimSSP2, n = 1000) %>% # slice top peaks by pvalue
  separate(col = gene, into = c("gene_1", "gene_2", "gene_3", "gene_4"), sep = ",")
SSP2S_top_targets <- c(SSP2S_top_targets_tb$gene_1, SSP2S_top_targets_tb$gene_2, SSP2S_top_targets_tb$gene_3, SSP2S_top_targets_tb$gene_4)
SSP2S_top_targets <- unique(SSP2S_top_targets[!is.na(SSP2S_top_targets)])
SSP2S_top_targets <- gsub("\\..*","",SSP2S_top_targets)  
length(SSP2S_top_targets) #1206 for top 1000 peaks


# private targets
SSP2F_priv_targets <- SlycSSP2_private_genes
length(SSP2F_priv_targets) #230

SSP2S_priv_targets <- SpimSSP2_private_genes
length(SSP2S_priv_targets) #3046

# private top targets
SSP2S_priv_top_targets_tb <- sig_SpimSSP2_vs_SlycSSP2_up_tb %>%
  filter(gene != "NA") %>% # filter for genic peaks
  slice_min(PValue.SpimSSP2_vs_SlycSSP2, n = 1000) %>% # slice top peaks by pvalue
  separate(col = gene, into = c("gene_1", "gene_2", "gene_3", "gene_4"), sep = ",")
SSP2S_priv_top_targets <- c(SSP2S_priv_top_targets_tb$gene_1, SSP2S_priv_top_targets_tb$gene_2, SSP2S_priv_top_targets_tb$gene_3, SSP2S_priv_top_targets_tb$gene_4)
SSP2S_priv_top_targets <- unique(SSP2S_priv_top_targets[!is.na(SSP2S_priv_top_targets)])
SSP2S_priv_top_targets <- gsub("\\..*","",SSP2S_priv_top_targets)  
length(SSP2S_priv_top_targets) #1232 for top 1000 peaks

# direct DEGs
length(direct_deg) #591

# select dataset name and genelist

dataset <- "SSP2S_priv_targets"
genelist <- SSP2S_priv_targets

# dataset <- "SSP2F_priv_targets"
# genelist <- SSP2F_priv_targets # no enrichged GO found

# define "universe" as all genes in genome
genome <- desc$geneid
length(genome) #35614

pvalueCutoff <- 0.05
OrgDb <- "org.SlycopersicumSL40VIBS.eg.db"
category <- "BP" # select GO category (BP, CC, MF, ALL)
ego <- enrichGO(gene          = genelist,
                universe      = genome,
                OrgDb         = OrgDb,
                keyType       = 'GID',
                ont           = category,
                pAdjustMethod = "BH", # select statistical method (see ClusterProfiler vignette)
                pvalueCutoff  = pvalueCutoff,
                # qvalueCutoff  = qvalueCutoff,
                readable      = FALSE)
head(ego)

barplot(
  ego,
  # x = "GeneRatio",
  color = "pvalue", #"pvalue", "p.adjust", "qvalue"
  showCategory = 10,
  size = NULL,
  split = NULL,
  font.size = 10,
  title = paste0(dataset,"\n",OrgDb),
  orderBy = "x",
  label_format = 40) + scale_fill_viridis(direction = -1)

ggsave((file.path(paste0(date, "_", experiment, "/ClusterProfiler/enrichGO_",dataset,"_fc",fc,"_p",pvalueCutoff,"_",category,"_",OrgDb,".pdf"))), width = 5, height = 3.5, useDingbats=FALSE)

# Extract results and save to file
ego.txt <- ego@result
write.csv(as.data.frame(ego.txt), file.path(paste0(date, "_", experiment, "/ClusterProfiler/enrichGO_",dataset,"_fc",fc,"_p",pvalueCutoff,"_",category,"_",OrgDb,".csv")), row.names = FALSE)
