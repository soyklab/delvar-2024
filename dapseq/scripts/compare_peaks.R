#-------------------------------------------------------------------------------------------
#
# PURPOSE: Compare Chip-seq peaks and visualize peak annotations using ChIPseeker  
#
# VERSION: 1
#
# AUTHOR : Giovanna Ambrosini
#          giovanna.ambrosini@epfl.ch 
#
#
# INPUT :  List of tsv-formatted annoated peak files from Differential binding analysis
# OUTPUT:  Venn diagrams representing peak overlaps, peak tables, heatmaps of top peaks for
#          several combinations of overlapping regions, and genomic annotation visualization
#          plots by ChIPseeker (Profile plots, distribution plots around TSS) 
#
#
#---------------------------------------------------------------------------------------------
#
library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library(csaw)
library(eulerr)
library(ggplot2)
library(fastcluster)
library(stringr)                # str_split_fixed()
library(writexl)
library(rapportools)            # is.empty()
library(cowplot)
library(ggrepel)
library(RColorBrewer)
library(clusterProfiler)
library(viridis)
library(ChIPseeker)
library(ComplexHeatmap)

library(extrafont)
#font_import()
loadfonts()
fonts()
citation("csaw")
citation("clusterProfiler")
citation("ChIPseeker")

setwd("/Users/giovanna/tomato_bZIP_sebastien_soyk/")
# Setiting some parameters                                                           ####
chr_to_ignore=unlist(strsplit("chrM,chrY",","))
chr_to_ignore

# Comparing common peaks (vs input)                                                  ####
outputfile="analysis/csaw_DB_comparison/comparison_SSP_vs_SlycSSP2_vs_SpimSSP2_peaks_fdr_0.01.pdf"
outputfile="analysis/csaw_DB_comparison/comparison_SSP_vs_SlycSSP2_vs_SpimSSP2_peaks_fdr_0.05.pdf"
directions=c("up","up", "up")
labels=c("SSP peaks","SlycSSP2 peaks", "SpimSSP2 peaks")
inputs_fdr_0.01 = c("analysis/SSP_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv",
           "analysis/SlycSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv",
           "analysis/SpimSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv")

inputs_fdr_0.05 = c("analysis/SSP_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv",
                    "analysis/SlycSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv",
                    "analysis/SpimSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv")
inputs <- inputs_fdr_0.01
inputs <- inputs_fdr_0.05

# Comparing DB                                                                       ####
##
# SSP vs SlycSSP2                                                                    ####
inputs_db_fdr_0.01 = c("analysis/SSP_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv",
           "analysis/SlycSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv",
           "analysis/SpimSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv",
           "analysis/SSP_vs_SlycSSP2_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv"
           )
# SSP vs SpimSSP2 has 0 db regions for frd=0.01                                      ####
outputfile="analysis/csaw_DB_comparison/comparison_SSP_vs_SlycSSP2_DB_fdr_0.01.pdf"
directions=c("up","up", "up", "all")
labels=c("SSP peaks","SlycSSP2 peaks", "SpimSSP2 peaks", "SSP vs SlycSSP2 (DB)")

inputs <- inputs_db_fdr_0.01

#  SSP vs SlycSSP2 vs SpimSSP2   fdr=0.05                                           ####
inputs_db_fdr_0.05 = c("analysis/SSP_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv",
                       "analysis/SlycSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv",
                       "analysis/SpimSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv",
                       "analysis/SSP_vs_SlycSSP2_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv",
                       "analysis/SSP_vs_SpimSSP2_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv"
)
outputfile="analysis/csaw_DB_comparison/comparison_SSP_vs_SlycSSP2_DB_SSP_vs_SpimSSP2_DB_fdr_0.05.pdf"
directions=c("up","up", "up", "all", "all")
labels=c("SSP peaks","SlycSSP2 peaks", "SpimSSP2 peaks", "SSP vs SlycSSP2 (DB)", "SSP vs SpimSSP2 (DB)")

inputs <- inputs_db_fdr_0.05

# SpimSSP vs SlycSSP2 for frd=0.01                                                   ####
inputs_db_fdr_0.01 = c("analysis/SpimSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv",
                       "analysis/SlycSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv",
                       "analysis/SpimSSP2_vs_SlycSSP2_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv"
)
outputfile="analysis/csaw_DB_comparison/comparison_SpimSSP2_vs_SlycSSP2_DB_fdr_0.01.pdf"
directions=c("up","up", "all")
labels=c("SpimSSP2 peaks","SlycSSP2 peaks", "SpimSSP2 vs SlycSSP2 (DB)")

inputs <- inputs_db_fdr_0.01

# SpimSSP vs SlycSSP2 for frd=0.05                                                   ####
inputs_db_fdr_0.05 = c("analysis/SpimSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv",
                       "analysis/SlycSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv",
                       "analysis/SpimSSP2_vs_SlycSSP2_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv"
)
outputfile="analysis/csaw_DB_comparison/comparison_SpimSSP2_vs_SlycSSP2_DB_fdr_0.05.pdf"
directions=c("up","up", "all")
labels=c("SpimSSP2 peaks","SlycSSP2 peaks", "SpimSSP2 vs SlycSSP2 (DB)")

inputs <- inputs_db_fdr_0.05

#----------------------------------------------------------------------------------------
# Load data                                                                          ####
#----------------------------------------------------------------------------------------
data=lapply(seq_along(inputs),function(i){
    cat("reading",inputs[i],"\n")

    tmp=fread(inputs[i])
    cat(" filtering out chromosomes:",chr_to_ignore,"(",tmp[chr%in%chr_to_ignore,.N],"regions)\n")
    tmp=tmp[!chr%in%chr_to_ignore]
    if(directions[i]=="up")
    {
        ##keep only "up" peaks with positive fold change
        cat(" keeping only regions changing in \"up\" direction\n")
        tmp=tmp[direction=="up"&best.logFC>0]
    }
    if(directions[i]=="down")
    {
        ##keep only "down" peaks with negative fold change
        cat(" keeping only regions changing in \"down\" direction\n")
        tmp=tmp[direction=="down"&best.logFC<0]
    }
    datatmp=GRanges(seqnames=tmp$chr,ranges=IRanges(start=tmp$start,end=tmp$end))
    mcols(datatmp)=tmp[,.(PValue,FDR,best.logFC,overlap)]
    datatmp
})
head(data)
# Get all seq levels                                                                  ####
#?seqlevels
## BASIC USAGE OF THE seqlevels() GETTER AND SETTER
## Operations between 2 or more objects containing genomic ranges (e.g.
## finding overlaps, comparing, or matching) only make sense if the
## operands have the same seqlevels.
## So before performing such operations, it is often necessary to adjust
## the seqlevels in the operands so that they all have the same seqlevels.
## This is typically done with the seqlevels() setter.
## The setter can be used to rename, drop, add and/or reorder seqlevels of an object.
##
all.seqlevels=unique(unlist(lapply(data,seqlevels)))
## Use same seq levels everywhere
for(i in seq_along(data)) {  # 1,2
    seqlevels(data[[i]])=all.seqlevels
}
### Use max end pos as dummy chromome length
##
# Compute max position of every chromosome in both Granges objects (the 2 or more files)
max.seqlen=sapply(all.seqlevels,function(s) {
    max(sapply(data,function(x){max(end(x[seqnames(x)==s]))}))
})
for(i in seq_along(data)) {
    seqlengths(data[[i]])=max.seqlen[names(seqlengths(data[[i]]))]
}
max.seqlen
seqlengths(data[[1]])
seqlengths(data[[2]])
seqlengths(data[[3]])
seqlengths(data[[4]])
seqlengths(data[[5]])

#----------------------------------------------------------------------------------------
# Find overlaps for Venn Diagram                                                     ####
#----------------------------------------------------------------------------------------
## find overlap (in bp) between regions
#
### Makes one data.table from a list of many (rbindlist)
#
# seq_along(inputs)
# [1] 1 2
#
# combn {utils} Generate All Combinations of i Elements, Taken n (1,2,...) at a Time
#
data.overlap=rbindlist(lapply(seq_along(inputs),function(n){
    cat("n",n,"\n")
    rbindlist(lapply(combn(seq_along(inputs),n,simplify=FALSE),function(cmb){
        cat(" cmb",cmb,"\n")
        set1=cmb
        set2=setdiff(seq_along(inputs),cmb)
        region=data[[set1[1]]]
        ##keep regions inside all data[set1]
        for(x in set1)
        {
            region=intersect(region,data[[x]])
        }
        ##keep regions not inside any data[set1]
        for(x in set2)
        {
            region=setdiff(region,data[[x]])
        }
        l=sum(width(reduce(region)))
        data.table(t(setNames(seq_along(inputs)%in%set1,paste0("in",seq_along(inputs)))),overlap.bp=l)
    }))
}))

setnames(data.overlap,paste0("in",seq_along(inputs)),labels)
data.overlap[,overlap.percent:=100*overlap.bp/data.overlap[,sum(overlap.bp)]]

data.overlap


cat("creating",outputfile,"\n")
cairo_pdf(outputfile,10,12, pointsize=10, onefile =TRUE) #cairo_pdf to avoid "-" transformed to "--"

# Venn Diagram                                                                      ####
vd=euler(data.overlap[,.SD,.SDcols=labels],weights=data.overlap[,overlap.bp])
#plot(vd,quantities=list(fontsize=8,col="black",type="percent"),legend = list(side="bottom"))#,labels = list(font = 2,fontsize=8)
eulerr_options(quantities=list(fontsize=24, col="black"), pointsize=24)
#par(mar=c(5,7,2,3))
plot(vd, quantities=paste0(signif(data.overlap[,overlap.percent],2),"%"),
     legend = list(side="bottom"), adjust_labels = TRUE) #,labels = list(font = 2, fontsize=8)

dev.off()
### Done                                                                             ####
cat("Done 1\n")

#----------------------------------------------------------------------------------------
#  Compute Intersections for DB regions  (Violin plots)                              ####
#----------------------------------------------------------------------------------------
### FDR=0.05
inputs_fdr_0.05 = c("analysis/SSP_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv",
                    "analysis/SlycSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv",
                    "analysis/SpimSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv",
                    "analysis/SSP_vs_SlycSSP2_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv",
                    "analysis/SSP_vs_SpimSSP2_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv",
                    "analysis/SpimSSP2_vs_SlycSSP2_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv"
)
directions=c("up","up", "up", "all", "all", "all")
labels=c("SSP peaks","SlycSSP2 peaks", "SpimSSP2 peaks", "SSP vs SlycSSP2 (DB)",  "SSP vs SpimSSP2 (DB)", "SpimSSP2 vs SlycSSP2 (DB)")
### FDR=0.01
inputs_fdr_0.01 = c("analysis/SSP_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv",
                    "analysis/SlycSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv",
                    "analysis/SpimSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv",
                    "analysis/SSP_vs_SlycSSP2_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv",
                    "analysis/SpimSSP2_vs_SlycSSP2_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv"
)
directions=c("up","up", "up", "all", "all")
labels=c("SSP peaks","SlycSSP2 peaks", "SpimSSP2 peaks", "SSP vs SlycSSP2 (DB)", "SpimSSP2 vs SlycSSP2 (DB)")
inputs <- inputs_fdr_0.01

#----------------------------------------------------------------------------------------
# Load data                                                                          ####
#----------------------------------------------------------------------------------------
data=lapply(seq_along(inputs),function(i){
    cat("reading",inputs[i],"\n")

    tmp=fread(inputs[i])
    cat(" filtering out chromosomes:",chr_to_ignore,"(",tmp[chr%in%chr_to_ignore,.N],"regions)\n")
    tmp=tmp[!chr%in%chr_to_ignore]
    if(directions[i]=="up")
    {
        ##keep only "up" peaks with positive fold change
        cat(" keeping only regions changing in \"up\" direction\n")
        tmp=tmp[direction=="up"&best.logFC>0]
    }
    if(directions[i]=="down")
    {
        ##keep only "down" peaks with negative fold change
        cat(" keeping only regions changing in \"down\" direction\n")
        tmp=tmp[direction=="down"&best.logFC<0]
    }
    datatmp=GRanges(seqnames=tmp$chr,ranges=IRanges(start=tmp$start,end=tmp$end))
    mcols(datatmp)=tmp[,.(PValue,FDR,best.logFC,overlap)]
    datatmp
})
head(data)

# Intersect DB data with main peaks                                                  ####
labels
length(data[[4]])
#[1] 42372
lapply(data, length)

### SSP vs SlycSSP (DB)
### Use subsetByOverlaps to keep all columns
tmp <- subsetByOverlaps(data[[1]], data[[4]])
tmp2 <- subsetByOverlaps(data[[2]], data[[4]])
u <- union(tmp, tmp2)
# Add metadata columns
mcols(u)[match(tmp, u), colnames(values(tmp))] <- values(tmp)
mcols(u)[match(tmp2, u), colnames(values(tmp2))] <- values(tmp2)
length(u)
#GRanges object with 33936 ranges and 0 metadata columns
data[[4]] <- u
length(data[[4]])

### SpimSSP2 vs SlycSSP (DB)
tmp <- subsetByOverlaps(data[[2]], data[[5]])
tmp2 <- subsetByOverlaps(data[[3]], data[[5]])
u <- union(tmp, tmp2)
mcols(u)[match(tmp, u), colnames(values(tmp))] <- values(tmp)
mcols(u)[match(tmp2, u), colnames(values(tmp2))] <- values(tmp2)
length(u)
length(data[[5]])
data[[5]] <- u
length(data[[5]])

lapply(data, length)

df <- data.frame()
for(i in seq_along(data)) {
    my_df <- as.data.frame(data[[i]])
    colnames(my_df)=gsub("^seqnames","chr",colnames(my_df))
    my_df$sample <- gsub(" ", "_", labels[i])
    df <- rbind(df, my_df)
}
df
tail(df)
colnames(df)

### Violin plots (LogFC across samples)                                              ####
#
p <- ggplot(df, aes(x=sample, y=best.logFC)) +
    geom_violin(trim=FALSE)
p

df$sample <- gsub("_peaks", "_p", df$sample)
df$sample <- gsub("DB", "", df$sample)
df$sample <- gsub("[()]", "", df$sample)
df$sample <- gsub("SSP_vs_SlycSSP2_", "SSP_vs_SlycSSP2", df$sample)
df$sample <- gsub("SpimSSP2_vs_SlycSSP2_", "SpimSSP2_vs_SlycSSP2", df$sample)

# violin plot with mean points
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)
# violin plot with median points
p + stat_summary(fun.y=median, geom="point", size=2, color="red")

# Add mean and standard deviation
p + stat_summary(fun.data="mean_sdl", mult=1,
                 geom="crossbar", width=0.2 )

p + stat_summary(fun.data=mean_sdl, mult=1,
                 geom="pointrange", color="red")

data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
}

p + stat_summary(fun.data=data_summary) + geom_boxplot(width=0.1) + theme_minimal()

p <- ggplot(df, aes(x=sample, y=best.logFC, fill=sample)) +
    geom_violin(trim=FALSE) +
    stat_summary(fun.data=data_summary) + geom_boxplot(width=0.1) + theme_minimal() +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "darkred","darkgreen")) +
    theme(legend.position="top") + # Remove legend
    theme(
        axis.title.x = element_blank(),
        #axis.text.x = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 32),
        #panel.grid.major = element_blank(),
        panel.border = element_blank(),
        #legend.justification = c(1, 0),
        #legend.position = c(0.6, 0.7),
        #legend.title = element_text(size = 22),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        )

p
ggsave(filename = "analysis/csaw_DB_comparison/ViolinPlot_comp_logFC_fdr_0.01.png", plot = p, device = "png", width = 14, height = 10, dpi = 300)
ggsave(filename = "analysis/csaw_DB_comparison/ViolinPlot_comp_logFC_fdr_0.01.pdf", plot = p, device = "pdf", width = 14, height = 10, dpi = 300, useDingbats =F )

### Violin plots (Pval across samples)                                               ####
#
p <- ggplot(df, aes(x=sample, y=PValue, fill=sample)) +
    geom_violin(trim=FALSE) +
    stat_summary(fun.data=data_summary) + geom_boxplot(width=0.1) + theme_minimal() +
    scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "darkred","darkgreen")) +
    theme(legend.position="top") + # Remove legend
    theme(
        axis.title.x = element_blank(),
        #axis.text.x = element_text(size = 24),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24),
        axis.title.y = element_text(size = 32),
        #panel.grid.major = element_blank(),
        panel.border = element_blank(),
        #legend.justification = c(1, 0),
        #legend.position = c(0.6, 0.7),
        #legend.title = element_text(size = 22),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
    )

p
ggsave(filename = "analysis/csaw_DB_comparison/ViolinPlot_comp_Pval_fdr_0.01.png", plot = p, device = "png", width = 14, height = 10, dpi = 300)
ggsave(filename = "analysis/csaw_DB_comparison/ViolinPlot_comp_Pval_fdr_0.01.pdf", plot = p, device = "pdf", width = 14, height = 10, dpi = 300, useDingbats =F )

### Done                                                                              ####
cat("Done 2\n")

### Histograms                                                                        ####
x <- df[df$sample == "SSP_peaks", "best.logFC"]
hist(x, cex.axis=2, main=df[df$sample == "SSP_peaks", "sample"][1])
abline(v = mean(x),                       # Add line for mean
      col = "red",
      lwd = 3)

text(x = mean(x) * 1.7,                   # Add text for mean
     y = mean(x) * 1.7,
     paste("Mean =", round(mean(x), 3)),
     col = "red",
     cex = 2)

x <- df[df$sample == "SlycSSP2_peaks", "best.logFC"]
hist(x, cex.axis=2, main=df[df$sample == "SlycSSP2_peaks", "sample"][1])
abline(v = mean(x),                       # Add line for mean
       col = "red",
       lwd = 3)

text(x = mean(x) * 1.7,                   # Add text for mean
     y = mean(x) * 1.7,
     paste("Mean =", round(mean(x), 3)),
     col = "red",
     cex = 2)

x <- df[df$sample == "SpimSSP2_peaks", "best.logFC"]
hist(x, cex.axis=2, main=df[df$sample == "SpimSSP2_peaks", "sample"][1])
abline(v = mean(x),                       # Add line for mean
       col = "red",
       lwd = 3)

text(x = mean(x) * 1.7,                   # Add text for mean
     y = mean(x) * 1.7,
     paste("Mean =", round(mean(x), 3)),
     col = "red",
     cex = 2)

x <- df[df$sample == "SSP_p", "PValue"]
hist(x, cex.axis=2, main=df[df$sample == "SSP_p", "sample"][1])
abline(v = mean(x),                       # Add line for mean
       col = "red",
       lwd = 3)

text(x = mean(x) * 1.7,                   # Add text for mean
     y = mean(x) * 1.7,
     paste("Mean =", round(mean(x), 3)),
     col = "red",
     cex = 2)

x <- df[df$sample == "SlycSSP2_p", "PValue"]
hist(x, cex.axis=2, main=df[df$sample == "SlycSSP2_p", "sample"][1])
abline(v = mean(x),                       # Add line for mean
       col = "red",
       lwd = 3)

text(x = mean(x) * 1.7,                   # Add text for mean
     y = mean(x) * 1.7,
     paste("Mean =", round(mean(x), 5)),
     col = "red",
     cex = 2)

x <- df[df$sample == "SpimSSP2_p", "PValue"]
hist(x, cex.axis=2, main=df[df$sample == "SpimSSP2_p", "sample"][1])
abline(v = mean(x),                       # Add line for mean
       col = "red",
       lwd = 3)

text(x = mean(x) * 1.7,                   # Add text for mean
     y = mean(x) * 1.7,
     paste("Mean =", round(mean(x), 5)),
     col = "red",
     cex = 2)

### Done                                                                              ####
cat("Done 3\n")

#----------------------------------------------------------------------------------------
#  Prepare Summary table for all peak lists                                          ####
#----------------------------------------------------------------------------------------
### FDR=0.05
inputs_fdr_0.05 = c("analysis/SSP_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv",
                    "analysis/SlycSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv",
                    "analysis/SpimSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv",
                    "analysis/SSP_vs_SlycSSP2_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv",
                    "analysis/SSP_vs_SpimSSP2_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv",
                    "analysis/SpimSSP2_vs_SlycSSP2_filt_fc_4_minq.20_tol_100_fdr_0.05_db_anno_results.tsv"
)
directions=c("up","up", "up", "all", "all", "all")
labels=c("SSP_peaks","SlycSSP2_peaks", "SpimSSP2_peaks", "SSP_vs_SlycSSP2",  "SSP_vs_SpimSSP2", "SpimSSP2_vs_SlycSSP2")
### FDR=0.01
inputs_fdr_0.01 = c("analysis/SSP_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv",
                    "analysis/SlycSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv",
                    "analysis/SpimSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv",
                    "analysis/SSP_vs_SlycSSP2_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv",
                    "analysis/SpimSSP2_vs_SlycSSP2_filt_fc_4_minq.20_tol_100_fdr_0.01_db_anno_results.tsv"
)
directions=c("up","up", "up", "all", "all")
labels=c("SSP_peaks","SlycSSP2_peaks", "SpimSSP2_peaks", "SSP_vs_SlycSSP2", "SpimSSP2_vs_SlycSSP2")
inputs <- inputs_fdr_0.01

###
#----------------------------------------------------------------------------------------
# Load data                                                                          ####
#----------------------------------------------------------------------------------------
data=lapply(seq_along(inputs),function(i){
    cat("reading",inputs[i],"\n")

    tmp=fread(inputs[i])
    cat(" filtering out chromosomes:",chr_to_ignore,"(",tmp[chr%in%chr_to_ignore,.N],"regions)\n")
    tmp=tmp[!chr%in%chr_to_ignore]
    if(directions[i]=="up")
    {
        ##keep only "up" peaks with positive fold change
        cat(" keeping only regions changing in \"up\" direction\n")
        tmp=tmp[direction=="up"&best.logFC>0]
    }
    if(directions[i]=="down")
    {
        ##keep only "down" peaks with negative fold change
        cat(" keeping only regions changing in \"down\" direction\n")
        tmp=tmp[direction=="down"&best.logFC<0]
    }
    datatmp=GRanges(seqnames=tmp$chr,ranges=IRanges(start=tmp$start,end=tmp$end))
    mcols(datatmp)=tmp[,.(PValue,FDR,best.logFC,overlap)]
    colnames(mcols(datatmp)) <- c("PValue", "FDR", "logFC", "Annotation")
    datatmp
})
head(data)

### Change Annotation Column                                                         ####
for (j in 1:length(data)) {
    an <- vector()
    an <- unname(sapply(mcols(data[[j]])$Annotation, function(x){
        l <- strsplit(unlist(strsplit(x, ",")), ":")
        gene <- ""
        strand <- ""
        type <- ""
        tmp <- ""
        if (length(l) != 0) {
            for (i in seq(length(l)) ){ gene<-paste(l[[i]][1], gene, sep = ",")}
            for (i in seq(length(l)) ){ strand<-paste(l[[i]][2], strand, sep = ",")}
            for (i in seq(length(l)) ){ type<-paste(l[[i]][3], type, sep = ",")}
            tmp <- paste(gsub('.{1}$', '', gene), gsub('.{1}$', '', strand), gsub('.{1}$', '', type))
        }
        tmp
    }))
    mcols(data[[j]])$Annotation <- an
}

# Intersect DB data with main peaks                                                  ####
labels
#[1] "SSP_peaks"            "SlycSSP2_peaks"       "SpimSSP2_peaks"       "SSP_vs_SlycSSP2"
#[5] "SpimSSP2_vs_SlycSSP2"
length(data[[4]])
#[1] 42372
lapply(data, length)

### SSP vs SlycSSP (DB) - extract valid peaks
### Use subsetByOverlaps to keep all columns
tmp <- subsetByOverlaps(data[[1]], data[[4]])
tmp2 <- subsetByOverlaps(data[[2]], data[[4]])
u <- union(tmp, tmp2)
#GRanges object with 33936 ranges and 0 metadata columns
# Add metadata columns
mcols(u)[match(tmp, u), colnames(values(tmp))] <- values(tmp)
mcols(u)[match(tmp2, u), colnames(values(tmp2))] <- values(tmp2)
length(u)
SSP_vs_SlycSSP2 <- u
length(data[[4]])
length(SSP_vs_SlycSSP2)

### SpimSSP2 vs SlycSSP (DB)
tmp <- subsetByOverlaps(data[[2]], data[[5]])
tmp2 <- subsetByOverlaps(data[[3]], data[[5]])
u <- union(tmp, tmp2)
mcols(u)[match(tmp, u), colnames(values(tmp))] <- values(tmp)
mcols(u)[match(tmp2, u), colnames(values(tmp2))] <- values(tmp2)
length(u)
#GRanges object with 33936 ranges and 0 metadata columns
length(data[[5]])
SpimSSP2_vs_SlycSSP2 <- u
length(SpimSSP2_vs_SlycSSP2)

### Common peaks between SSP and SlycSSP2
###
length(data[[1]])
length(data[[2]])
common_SSP_SlycSSP2 <-  subsetByOverlaps(data[[1]], data[[2]])
length(common_SSP_SlycSSP2)
SSP_specs <- data[[1]][-queryHits(findOverlaps(data[[1]], common_SSP_SlycSSP2, type="any")),]
length(SSP_specs)
head(SSP_specs)
SlycSSP2_specs <- data[[2]][-queryHits(findOverlaps(data[[2]], common_SSP_SlycSSP2, type="any")),]
head(SlycSSP2_specs)
length(SlycSSP2_specs)
labels

### Common peaks between SSP and SpimSSP2
###
length(data[[1]])
length(data[[3]])
common_SSP_SpimSSP2 <-  subsetByOverlaps(data[[1]], data[[3]])
length(common_SSP_SpimSSP2)
length(data[[1]][-queryHits(findOverlaps(data[[1]], common_SSP_SpimSSP2, type="any")),] )
length(SSP_specs)
head(SSP_specs)
SpimSSP2_specs <- data[[3]][-queryHits(findOverlaps(data[[3]], common_SSP_SpimSSP2, type="any")),]
head(SpimSSP2_specs)
length(SpimSSP2_specs)
labels

lapply(data, length)

length(SSP_vs_SlycSSP2)
length(SpimSSP2_vs_SlycSSP2)
length(common_SSP_SlycSSP2)
length(common_SSP_SpimSSP2)
length(SSP_specs)

labels
SSP_p_df <- as.data.frame(data[[1]])
SlycSSP2_p_df <- as.data.frame(data[[2]])
SpimSSP2_p_df <- as.data.frame(data[[3]])
SSP_vs_SlycSSP2_df <- as.data.frame(data[[4]])
SpimSSP2_vs_SlycSSP2_df <- as.data.frame(data[[5]])

common_SSP_SlycSSP2_df <- as.data.frame(common_SSP_SlycSSP2)
common_SSP_SpimSSP2_df <- as.data.frame(common_SSP_SpimSSP2)
SSP_spec_df <- as.data.frame(SSP_specs)
SlycSSP2_spec_df <-  as.data.frame(SlycSSP2_specs)
SSP_vs_SlycSSP2_valid_df <- as.data.frame(SSP_vs_SlycSSP2)
SpimSSP2_vs_SlycSSP2_valid_df <- as.data.frame(SpimSSP2_vs_SlycSSP2)

dim(SSP_p_df)
dim(SlycSSP2_p_df)
dim(SpimSSP2_p_df)
dim(SSP_vs_SlycSSP2_df)
dim(SpimSSP2_vs_SlycSSP2_df)


dim(common_SSP_SlycSSP2_df)
dim(common_SSP_SpimSSP2_df)
dim(SSP_spec_df)
dim(SlycSSP2_spec_df)
dim(SSP_vs_SlycSSP2_valid_df)
dim(SpimSSP2_vs_SlycSSP2_valid_df)

SSP_p_df$id <- paste(sep = "_", SSP_p_df$seqnames, SSP_p_df$start, SSP_p_df$end)
SlycSSP2_p_df$id <- paste(sep = "_", SlycSSP2_p_df$seqnames, SlycSSP2_p_df$start, SlycSSP2_p_df$end)
SpimSSP2_p_df$id <- paste(sep = "_", SpimSSP2_p_df$seqnames, SpimSSP2_p_df$start, SpimSSP2_p_df$end)
SSP_vs_SlycSSP2_df$id <- paste(sep = "_", SSP_vs_SlycSSP2_df$seqnames, SSP_vs_SlycSSP2_df$start, SSP_vs_SlycSSP2_df$end)
SpimSSP2_vs_SlycSSP2_df$id <- paste(sep = "_", SpimSSP2_vs_SlycSSP2_df$seqnames, SpimSSP2_vs_SlycSSP2_df$start, SpimSSP2_vs_SlycSSP2_df$end)

dim(SSP_p_df)
dim(SlycSSP2_p_df)
dim(SpimSSP2_p_df)
dim(SSP_vs_SlycSSP2_df)
dim(SpimSSP2_vs_SlycSSP2_df)

head(SSP_p_df)
head(SlycSSP2_p_df)
head(SpimSSP2_p_df)
head(SSP_vs_SlycSSP2_df)
head(SpimSSP2_vs_SlycSSP2_df)

### Merge data frames                                                                 ####
colnames(SSP_p_df)
head(SSP_p_df[, c(6,7,8,10)])
SSP_p_df_t <- SSP_p_df[, c(6,7,8,10)]
SlycSSP2_p_df_t <- SlycSSP2_p_df[, c(6,7,8,10)]
SpimSSP2_p_df_t <- SpimSSP2_p_df[, c(6,7,8,10)]
SSP_vs_SlycSSP2_df_t <- SSP_vs_SlycSSP2_df[, c(6,7,8,10)]
SpimSSP2_vs_SlycSSP2_df_t <- SpimSSP2_vs_SlycSSP2_df[, c(6,7,8,10)]

head(SSP_p_df_t)
head(SlycSSP2_p_df_t)
head(SpimSSP2_p_df_t)
head(SSP_vs_SlycSSP2_df_t)
head(SpimSSP2_vs_SlycSSP2_df_t)

dim(SSP_p_df_t)
dim(SlycSSP2_p_df_t)
dim(SpimSSP2_p_df_t)
dim(SSP_vs_SlycSSP2_df_t)
dim(SpimSSP2_vs_SlycSSP2_df_t)

dim(merge(SSP_p_df_t, SlycSSP2_p_df_t, all=TRUE, by="id"))
head(merge(SSP_p_df_t, SlycSSP2_p_df_t, all=TRUE, by="id"))

m1 <- merge(SSP_p_df_t, SlycSSP2_p_df_t, all=TRUE, by="id")
dim(m1)
head(m1)
colnames(m1)
colnames(m1) <- c("id", "PValue.SSP", "FDR.SPP", "logFC.SSP", "PValue.SlycSSP2", "FDR.SlycSSP2", "logFC.SlycSSP2")

m2 <- merge(m1, SpimSSP2_p_df_t, all=TRUE, by="id")
dim(m2)
head(m2)
colnames(m2) <- c("id", "PValue.SSP", "FDR.SPP", "logFC.SSP", "PValue.SlycSSP2", "FDR.SlycSSP2", "logFC.SlycSSP2", "PValue.SpimSSP2", "FDR.SpimSSP2", "logFC.SpimSSP2")
colnames(m2)

m3 <- merge(m2, SSP_vs_SlycSSP2_df_t, all=TRUE, by="id")
dim(m3)
head(m3)
colnames(m3) <- c("id", "PValue.SSP", "FDR.SPP", "logFC.SSP", "PValue.SlycSSP2", "FDR.SlycSSP2", "logFC.SlycSSP2",
                  "PValue.SpimSSP2", "FDR.SpimSSP2", "logFC.SpimSSP2", "PValue.SSP_vs_SlycSSP2",
                  "FDR.SSP_vs_SlycSSP2", "logFC.SSP_vs_SlycSSP2" )
colnames(m3)

m4 <- merge(m3, SpimSSP2_vs_SlycSSP2_df_t, all=TRUE, by="id")
dim(m4)
head(m4)
colnames(m4) <- c("id", "PValue.SSP", "FDR.SPP", "logFC.SSP", "PValue.SlycSSP2", "FDR.SlycSSP2", "logFC.SlycSSP2",
                  "PValue.SpimSSP2", "FDR.SpimSSP2", "logFC.SpimSSP2", "PValue.SSP_vs_SlycSSP2",
                  "FDR.SSP_vs_SlycSSP2", "logFC.SSP_vs_SlycSSP2", "PValue.SpimSSP2_vs_SlycSSP2",
                  "FDR.SpimSSP2_vs_SlycSSP2", "logFC.SpimSSP2_vs_SlycSSP2")

dim(Reduce(function(x,y) merge(x = x, y = y, by = "id", all=TRUE),
       list(SSP_p_df_t, SlycSSP2_p_df_t, SpimSSP2_p_df_t)))
head(Reduce(function(x,y) merge(x = x, y = y, by = "id", all=TRUE),
       list(SSP_p_df_t, SlycSSP2_p_df_t, SpimSSP2_p_df_t)))

rownames(m4) <- m4$id

sum(!is.na(m4$logFC.SSP_vs_SlycSSP2))

V_flag_SSP_vs_SlycSSP2 <- rep(1, dim(m4)[1])
for(i in seq(1:dim(m4)[1])) {
    if(is.na(m4$logFC.SSP_vs_SlycSSP2[i]) || (is.na(m4$logFC.SSP[i]) && is.na(m4$logFC.SlycSSP2[i]))) {
        V_flag_SSP_vs_SlycSSP2[i] <- 0
    }
}
sum(V_flag_SSP_vs_SlycSSP2 == 1)

V_flag_SpimSSP2_vs_SlycSSP2 <- rep(1, dim(m4)[1])
for(i in seq(1:dim(m4)[1])) {
    if(is.na(m4$logFC.SpimSSP2_vs_SlycSSP2[i]) || (is.na(m4$logFC.SpimSSP2[i]) && is.na(m4$logFC.SlycSSP2[i]))) {
        V_flag_SpimSSP2_vs_SlycSSP2[i] <- 0
    }
}
sum(V_flag_SpimSSP2_vs_SlycSSP2 == 1)

length(V_flag_SSP_vs_SlycSSP2)
length(V_flag_SpimSSP2_vs_SlycSSP2)

head(m4)
m4$T_flag_SSP_vs_SlycSSP2 <- V_flag_SSP_vs_SlycSSP2
m4$T_flag_SpimSSP2_vs_SlycSSP2 <- V_flag_SpimSSP2_vs_SlycSSP2
head(m4)

m4[c('chr', 'start', 'end')] <- str_split_fixed(m4$id, '_', 3)

head(m4)
colnames(m4)
colnames(m4[,c(19,20,21,2:13,17,14:16,18)])
m4 <- m4[,c(19,20,21,2:13,17,14:16,18)]
head(m4)
dim(m4)
#[1] 49170    20
m4[,2] <- as.integer(m4[,2])
m4[,3] <- as.integer(m4[,3])

#  Order m4 by chromosome and start position                                         ####
m4 <- m4[order( m4[,1], m4[,2] ),]

# Add Annotation                                                                     ####
#
gffFile <- "/Users/giovanna/tomato_bZIP_sebastien_soyk/data/genome/S100_v2.1.0/SollycSweet-100_genes_v2.1.1.new.gff3"

## Generate annotation object
txdb <- makeTxDbFromGFF(file=gffFile,
                        dataSource="Alonge et al., 2021",
                        organism="Solanum lycopersicum", dbxrefTag="Alias")

library(org.SlycopersicumNatalia.eg.db)
org.SlycopersicumNatalia.eg.db

m4_gr <- makeGRangesFromDataFrame(m4, keep.extra.columns = T, ignore.strand = T)
m4_gr
anno <- detailRanges(m4_gr, txdb=txdb,
                     orgdb=org.SlycopersicumNatalia.eg.db, key.field="GID", promoter=c(3000, 2000))
head(anno$overlap)

### Change annotation vector
an <- vector()
an <- unname(sapply(anno$overlap, function(x){
    l <- strsplit(unlist(strsplit(x, ",")), ":")
    gene <- ""
    strand <- ""
    type <- ""
    tmp <- ""
    if (length(l) != 0) {
        for (i in seq(length(l)) ){ gene<-paste(l[[i]][1], gene, sep = ",")}
        for (i in seq(length(l)) ){ strand<-paste(l[[i]][2], strand, sep = ",")}
        for (i in seq(length(l)) ){ type<-paste(l[[i]][3], type, sep = ",")}
        tmp <- paste(gsub('.{1}$', '', gene), gsub('.{1}$', '', strand), gsub('.{1}$', '', type))
    }
    tmp
}))
head(an)
length(an)
#[1] 49170
m4[c('gene', 'strand', 'type')] <- str_split_fixed(an, ' ', 3)

head(m4)
write.table(m4,sep="\t",quote=FALSE,row.names=FALSE,file="analysis/peak_table_fdr_0.01.tsv")

#### Save m4 as xlsx                                   ####
write_xlsx(m4,path="analysis/peak_table_fdr_0.01.xlsx")

### Done                                                                              ####
cat("Done 4\n")

# Heatmap of top peaks based on csaw normalized counts                                ####
# Import Data with Filtered and normalized counts                                     ####
###
### FDR=0.01
r_inputs_fdr_0.01 = c("analysis/SSP_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.01_filtered_data.rds",
                    "analysis/SlycSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.01_filtered_data.rds",
                    "analysis/SpimSSP2_vs_Input_filt_fc_4_minq.20_tol_100_fdr_0.01_filtered_data.rds",
                    "analysis/SSP_vs_SlycSSP2_filt_fc_4_minq.20_tol_100_fdr_0.01_filtered_data.rds",
                    "analysis/SpimSSP2_vs_SlycSSP2_filt_fc_4_minq.20_tol_100_fdr_0.01_filtered_data.rds"
)
directions=c("up","up", "up", "all", "all")
labels=c("SSP_peaks","SlycSSP2_peaks", "SpimSSP2_peaks", "SSP_vs_SlycSSP2", "SpimSSP2_vs_SlycSSP2")
inputs <- r_inputs_fdr_0.01

metadatafile <- "/Users/giovanna/tomato_bZIP_sebastien_soyk/data/metadata.txt"
getwd()

### Set bam files path                       ####
inputdir <- "/Users/giovanna/tomato_bZIP_sebastien_soyk/data/bams"
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

# Load filtered count data                                                            ####
### Read filtered data for calculating normalized counts                              ####
filtered.data <- readRDS('analysis/filtered_data.rds')
head(filtered.data)
mcols(filtered.data)
rowData(filtered.data)
#DataFrame with 391063 rows and 0 columns (the entire csaw-scanned regions/windows) - win.data -
# bin size = window size = 10 bp
assay(filtered.data)

metadata
#sample condition replicate
#/Users/giovanna/tomato_bZIP_sebastien_soyk/data/bams/Input_1.bam       Input_1     input         1
#/Users/giovanna/tomato_bZIP_sebastien_soyk/data/bams/Input_2.bam       Input_2     input         2
#/Users/giovanna/tomato_bZIP_sebastien_soyk/data/bams/SSP_1.bam           SSP_1       SSP         1
#/Users/giovanna/tomato_bZIP_sebastien_soyk/data/bams/SSP_2.bam           SSP_2       SSP         2
#/Users/giovanna/tomato_bZIP_sebastien_soyk/data/bams/SlycSSP2_1.bam SlycSSP2_1  SlycSSP2         1
#/Users/giovanna/tomato_bZIP_sebastien_soyk/data/bams/SlycSSP2_2.bam SlycSSP2_2  SlycSSP2         2
#/Users/giovanna/tomato_bZIP_sebastien_soyk/data/bams/SpimSSP2_1.bam SpimSSP2_1  SpimSSP2         1
#/Users/giovanna/tomato_bZIP_sebastien_soyk/data/bams/SpimSSP2_2.bam SpimSSP2_2  SpimSSP2         2

data.normalized.all.mat <- calculateCPM(filtered.data,use.norm.factors=TRUE,log=FALSE)
# create list
data.normalized.all=list()
out.data=rowRanges(filtered.data)
for(i in 1:ncol(filtered.data)) {
    mcols(out.data)=data.frame(score=data.normalized.all.mat[,i])
    data.normalized.all[[i]]=out.data
}
data.normalized.all

### Order SSP peaks by FDR                                                    ####
o <-order(SSP_p_df$FDR, decreasing=FALSE)
length(o)
#o <- order(differential.binding.regions$PValue)
head(SSP_p_df[o,])
dim(SSP_p_df)
#[1] 38239    9

### 500 top SSP peaks                                                         ####
#
SSP_p_df[o[500],]
#      chr  start    end       PValue          FDR    logFC               gene type strand
#34945 chr8 65562301 65563110 1.554068e-08 1.808055e-06 6.835306 Solyc08g081620.4.1   PE      -
top_SSP_p <- SSP_p_df[o[1:500],]
dim(top_SSP_p)
#[1] 500  9
head(top_SSP_p)

### 500 top SlycSSP2 peaks                                                    ####
#
o <-order(SlycSSP2_p_df$FDR, decreasing=FALSE)
length(o)
head(SlycSSP2_p_df[o,])
dim(SlycSSP2_p_df)
#[1] 5475   9
SlycSSP2_p_df[o[500],]
top_SlycSSP2_p <- SlycSSP2_p_df[o[1:500],]
dim(top_SlycSSP2_p)
head(top_SlycSSP2_p)

### 500 top SlycSSP2 peaks                                                    ####
#
o <-order(SpimSSP2_p_df$FDR, decreasing=FALSE)
length(o)
head(SpimSSP2_p_df[o,])
dim(SpimSSP2_p_df)
SpimSSP2_p_df[o[500],]
top_SpimSSP2_p <- SpimSSP2_p_df[o[1:500],]
dim(top_SpimSSP2_p)
head(top_SpimSSP2_p)

top_SSP_p$id <- paste(sep = "_", top_SSP_p$chr, top_SSP_p$start, top_SSP_p$end)
top_SlycSSP2_p$id <- paste(sep = "_", top_SlycSSP2_p$chr, top_SlycSSP2_p$start, top_SlycSSP2_p$end)
top_SpimSSP2_p$id <- paste(sep = "_", top_SpimSSP2_p$chr, top_SpimSSP2_p$start, top_SpimSSP2_p$end)

top_SSP_p <- data.frame(id=top_SSP_p[, c(10)])
top_SlycSSP2_p <- data.frame(id=top_SlycSSP2_p[, c(10)])
top_SpimSSP2_p <- data.frame(id=top_SpimSSP2_p[, c(10)])

dim(merge(top_SSP_p, top_SlycSSP2_p, all=TRUE, by="id"))
m1 <- merge(top_SSP_p, top_SlycSSP2_p, all=TRUE, by="id")
dim(m1)
m2 <- merge(m1, top_SpimSSP2_p, all=TRUE, by="id")
dim(m2)
colnames(m2)
head(m2)
data.normalized.all

f <- str_split_fixed(m2$id[1], '_', 3)
cur.region <- GRanges(f[1], IRanges(as.integer(f[2]), as.integer(f[3])))
SSP_1 <- sum(subsetByOverlaps(data.normalized.all[[3]], cur.region)$score)
SSP_2 <- sum(subsetByOverlaps(data.normalized.all[[4]], cur.region)$score)
SlycSSP2_1 <-sum(subsetByOverlaps(data.normalized.all[[5]], cur.region)$score)
SlycSSP2_2 <-sum(subsetByOverlaps(data.normalized.all[[6]], cur.region)$score)
SpimSSP2_1 <- sum(subsetByOverlaps(data.normalized.all[[7]], cur.region)$score)
SpimSSP2_2 <- sum(subsetByOverlaps(data.normalized.all[[8]], cur.region)$score)

###
SSP_1 <- vector()
SSP_2 <- vector()
SlycSSP2_1 <- vector()
SlycSSP2_2 <- vector()
SpimSSP2_1 <- vector()
SpimSSP2_2 <- vector()
for(i in seq(length(m2$id))) {
    f <- str_split_fixed(m2$id[i], '_', 3)
    cur.region <- GRanges(f[1], IRanges(as.integer(f[2]), as.integer(f[3])))
    SSP_1[i] <- sum(subsetByOverlaps(data.normalized.all[[3]], cur.region)$score)
    SSP_2[i] <- sum(subsetByOverlaps(data.normalized.all[[4]], cur.region)$score)
    SlycSSP2_1[i] <-sum(subsetByOverlaps(data.normalized.all[[5]], cur.region)$score)
    SlycSSP2_2[i] <-sum(subsetByOverlaps(data.normalized.all[[6]], cur.region)$score)
    SpimSSP2_1[i] <- sum(subsetByOverlaps(data.normalized.all[[7]], cur.region)$score)
    SpimSSP2_2[i] <- sum(subsetByOverlaps(data.normalized.all[[8]], cur.region)$score)
}
length(SSP_1)
length(SSP_2)
length(SlycSSP2_1)
length(SlycSSP2_2)
length(SpimSSP2_1)
length(SpimSSP2_2)

head(cbind(m2, SSP_1 = SSP_1, SSP_2 = SSP_2, SlycSSP2_1 = SlycSSP2_1, SlycSSP2_2 = SlycSSP2_2,
      SpimSSP2_1 = SpimSSP2_1, SpimSSP2_2 = SpimSSP2_2))

m2 <- cbind(m2, SSP_1 = SSP_1, SSP_2 = SSP_2, SlycSSP2_1 = SlycSSP2_1, SlycSSP2_2 = SlycSSP2_2,
            SpimSSP2_1 = SpimSSP2_1, SpimSSP2_2 = SpimSSP2_2)
head(m2)
norm_counts <- m2[, -c(1)]

head(norm_counts)
dim(norm_counts)

top_peak_table <- m2
head(top_peak_table)
top_peak_table[c('chr', 'start', 'end')] <- str_split_fixed(top_peak_table$id, '_', 3)
head(top_peak_table)
top_peak_table <- top_peak_table[, c(1,8:10,2:7)]
dim(top_peak_table)
write.table(top_peak_table,sep="\t",quote=FALSE,row.names=FALSE,file="analysis/Top-500peaks_norm_counts.tsv")
write.table(top_peak_table,sep="\t",quote=FALSE,row.names=FALSE,file="results/Top-500peaks_norm_counts.tsv")

### Generate heatmap of 500-top significant SSP/SSP2 peaks                        ####
#
### Define some graphical parameters
ht_opt$legend_title_gp = gpar(fontsize = 20)
ht_opt$legend_labels_gp = gpar(fontsize = 16)
ht_opt$simple_anno_size = unit(5, "mm")
ht_opt(HEATMAP_LEGEND_PADDING = unit(3, "mm"), ADD = TRUE)
ht_opt(ANNOTATION_LEGEND_PADDING = unit(3, "mm"), ADD = TRUE)

ha_column1 = HeatmapAnnotation(gp = gpar(col="black", fontsize=20),
                               col=list(Samples=c(SSP_1='red', SSP_2='darkred', SlycSSP2_1='green', SlycSSP2_2='darkgreen',
                                                  SpimSSP2_1='orange', SpimSSP2_2='orange3')),
                               Samples=c('SSP_1','SSP_2','SlycSSP2_1', 'SlycSSP2_2', 'SpimSSP2_1', 'SpimSSP2_2'),
                               annotation_name_gp= gpar(fontsize = 20)
                               )
my_alpha <- c(rep(0, dim(norm_counts)[1]))
# Heatmaps                                                                            ####
hmap<-Heatmap(t(scale(t(norm_counts))), name = "Z-Score",
              clustering_distance_rows = "pearson",  clustering_distance_columns = "pearson",
              top_annotation = ha_column1,
              #show_row_names = TRUE,
              show_row_names = FALSE,
              row_names_gp = gpar(col = "blue",
                                  alpha=my_alpha,
                                  fontsize = 10, fontface="bold"),
              cluster_columns = TRUE, cluster_rows = TRUE,
              show_column_names = FALSE)

#  Draw heatmap
draw(hmap, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)
decorate_heatmap_body("Z-Score", {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))})

png("results/visualization/Top-500peaks_6_samples_heatmap.png", width=1000, height=1000, res = 150)
draw(hmap, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)
decorate_heatmap_body("Z-Score",{grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))})
dev.off()

pdf("results/visualization/Top-500peaks_6_samples_heatmap.pdf", width=10, height=10, useDingbats =F)
draw(hmap, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)
decorate_heatmap_body("Z-Score",{grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))})
dev.off()

#### Heatmap for top 500 most significant annotated peaks from all samples            ####
#
### Order SSP peaks by FDR and filter the top 500 peaks                               ####
colnames(SSP_p_df)
top_SSP_p  <- SSP_p_df %>% arrange(FDR, decreasing=FALSE) %>% filter(!is.empty(gene)) %>% slice(1:500)
dim(top_SSP_p)
#[1] 500  9
head(top_SSP_p, 20)

### 500 top SlycSSP2 annotated peaks                                                             ####
#
top_SlycSSP2_p  <- SlycSSP2_p_df %>% arrange(FDR, decreasing=FALSE) %>% filter(!is.empty(gene)) %>% slice(1:500)
dim(top_SlycSSP2_p)
#[1] 500  9
head(top_SlycSSP2_p)

### 500 top SlycSSP2 peaks                                                           ####
#
top_SpimSSP2_p  <- SpimSSP2_p_df %>% arrange(FDR, decreasing=FALSE) %>% filter(!is.empty(gene)) %>% slice(1:500)
dim(top_SpimSSP2_p)
#[1] 500  10
head(top_SpimSSP2_p)

top_SSP_p$id <- paste(sep = "_", top_SSP_p$chr, top_SSP_p$start, top_SSP_p$end)
top_SlycSSP2_p$id <- paste(sep = "_", top_SlycSSP2_p$chr, top_SlycSSP2_p$start, top_SlycSSP2_p$end)
top_SpimSSP2_p$id <- paste(sep = "_", top_SpimSSP2_p$chr, top_SpimSSP2_p$start, top_SpimSSP2_p$end)

top_SSP_p.sav <- top_SSP_p
top_SlycSSP2_p.sav <- top_SlycSSP2_p
top_SpimSSP2_p.sav <- top_SpimSSP2_p

colnames(top_SSP_p)
colnames(top_SSP_p.sav)
top_SSP_p <- data.frame(id=top_SSP_p[, c(10)])
top_SlycSSP2_p <- data.frame(id=top_SlycSSP2_p[, c(10)])
top_SpimSSP2_p <- data.frame(id=top_SpimSSP2_p[, c(10)])

dim(merge(top_SSP_p, top_SlycSSP2_p, all=TRUE, by="id"))
m1 <- merge(top_SSP_p, top_SlycSSP2_p, all=TRUE, by="id")
dim(m1)
m2 <- merge(m1, top_SpimSSP2_p, all=TRUE, by="id")
dim(m2)
colnames(m2)
head(m2)
data.normalized.all

###
SSP_1 <- vector()
SSP_2 <- vector()
SlycSSP2_1 <- vector()
SlycSSP2_2 <- vector()
SpimSSP2_1 <- vector()
SpimSSP2_2 <- vector()
for(i in seq(length(m2$id))) {
    f <- str_split_fixed(m2$id[i], '_', 3)
    cur.region <- GRanges(f[1], IRanges(as.integer(f[2]), as.integer(f[3])))
    SSP_1[i] <- sum(subsetByOverlaps(data.normalized.all[[3]], cur.region)$score)
    SSP_2[i] <- sum(subsetByOverlaps(data.normalized.all[[4]], cur.region)$score)
    SlycSSP2_1[i] <-sum(subsetByOverlaps(data.normalized.all[[5]], cur.region)$score)
    SlycSSP2_2[i] <-sum(subsetByOverlaps(data.normalized.all[[6]], cur.region)$score)
    SpimSSP2_1[i] <- sum(subsetByOverlaps(data.normalized.all[[7]], cur.region)$score)
    SpimSSP2_2[i] <- sum(subsetByOverlaps(data.normalized.all[[8]], cur.region)$score)
}
length(SSP_1)
length(SSP_2)
length(SlycSSP2_1)
length(SlycSSP2_2)
length(SpimSSP2_1)
length(SpimSSP2_2)
#[1] 844

head(cbind(m2, SSP_1 = SSP_1, SSP_2 = SSP_2, SlycSSP2_1 = SlycSSP2_1, SlycSSP2_2 = SlycSSP2_2,
           SpimSSP2_1 = SpimSSP2_1, SpimSSP2_2 = SpimSSP2_2))

m2 <- cbind(m2, SSP_1 = SSP_1, SSP_2 = SSP_2, SlycSSP2_1 = SlycSSP2_1, SlycSSP2_2 = SlycSSP2_2,
            SpimSSP2_1 = SpimSSP2_1, SpimSSP2_2 = SpimSSP2_2)
head(m2)
norm_counts <- m2[, -c(1)]

head(norm_counts)
dim(norm_counts)
#[1] 844   6

top_peak_table <- m2
head(top_peak_table)

### Add annotation
colnames(top_SSP_p.sav)
top_SSP_p.sav <- top_SSP_p.sav[,c(10,7:9)]
top_SlycSSP2_p.sav <- top_SlycSSP2_p.sav[,c(10,7:9)]
top_SpimSSP2_p.sav <- top_SpimSSP2_p.sav[,c(10,7:9)]

head(top_SSP_p.sav)
head(top_peak_table)
dim(top_SSP_p.sav)
dim(top_peak_table)
head(merge(top_SSP_p.sav, top_SlycSSP2_p.sav,  all=TRUE, by="id"), 50)
m <- merge(top_SSP_p.sav, top_SlycSSP2_p.sav,  all=TRUE, by="id")
dim(m)
#
# copy cols 5,6,7 into col 2,3,4 when 2,3,4=<NA>
m[is.na(m[, 2]), c(2:4)] <- m[is.na(m[, 2]), c(5:7)]
head(m,50)
m <- m[, c(1:4)]
colnames(m)[2:4] <- c("gene", "type", "strand")
head(m,50)
#
# Merge with top_SpimSSP2_p.sav
m1 <- merge(m, top_SpimSSP2_p.sav,  all=TRUE, by="id")
dim(m1)
head(m1,50)
m1[is.na(m1[, 2]), c(2:4)] <- m1[is.na(m1[, 2]), c(5:7)]
m1 <- m1[, c(1:4)]
colnames(m1)[2:4] <- c("gene", "type", "strand")
head(m1,50)
dim(m1)

# Merge Annotation (m1) with top_peak_table
head(top_peak_table)
dim(top_peak_table)
top_peak_table <- merge(top_peak_table, m1, by="id")
dim(top_peak_table)
head(top_peak_table)

top_peak_table[c('chr', 'start', 'end')] <- str_split_fixed(top_peak_table$id, '_', 3)
head(top_peak_table)
top_peak_table <- top_peak_table[, c(1,11:13,2:10)]
dim(top_peak_table)
#[1] 844  13
write.table(top_peak_table,sep="\t",quote=FALSE,row.names=FALSE,file="analysis/Top-500peaks_anno_norm_counts.tsv")
write.table(top_peak_table,sep="\t",quote=FALSE,row.names=FALSE,file="results/Top-500peaks_anno_norm_counts.tsv")
#### Save as xlsx
write_xlsx(top_peak_table,path="results/Top-500peaks_anno_norm_counts.xlsx")
write_xlsx(top_peak_table,path="analysis/Top-500peaks_anno_norm_counts.xlsx")

### Generate heatmap of 500-top significant SSP/SSP2 peaks                           ####
#
### Define some graphical parameters
ht_opt$legend_title_gp = gpar(fontsize = 20)
ht_opt$legend_labels_gp = gpar(fontsize = 16)
ht_opt$simple_anno_size = unit(5, "mm")
ht_opt(HEATMAP_LEGEND_PADDING = unit(3, "mm"), ADD = TRUE)
ht_opt(ANNOTATION_LEGEND_PADDING = unit(3, "mm"), ADD = TRUE)

ha_column1 = HeatmapAnnotation(gp = gpar(col="black", fontsize=20),
                               col=list(Samples=c(SSP_1='red', SSP_2='darkred', SlycSSP2_1='green', SlycSSP2_2='darkgreen',
                                                  SpimSSP2_1='orange', SpimSSP2_2='orange3')),
                               Samples=c('SSP_1','SSP_2','SlycSSP2_1', 'SlycSSP2_2', 'SpimSSP2_1', 'SpimSSP2_2'),
                               annotation_name_gp= gpar(fontsize = 20)
)
my_alpha <- c(rep(0, dim(norm_counts)[1]))
# Heatmap
hmap<-Heatmap(t(scale(t(norm_counts))), name = "Z-Score",
              clustering_distance_rows = "pearson",  clustering_distance_columns = "pearson",
              top_annotation = ha_column1,
              #show_row_names = TRUE,
              show_row_names = FALSE,
              row_names_gp = gpar(col = "blue",
                                  alpha=my_alpha,
                                  fontsize = 10, fontface="bold"),
              cluster_columns = TRUE, cluster_rows = TRUE,
              show_column_names = FALSE)

#  Draw heatmap
draw(hmap, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)
decorate_heatmap_body("Z-Score", {grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))})

png("results/visualization/Top-500peaks_anno_6_samples_heatmap.png", width=1000, height=1000, res = 150)
draw(hmap, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)
decorate_heatmap_body("Z-Score",{grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))})
dev.off()

pdf("results/visualization/Top-500peaks_anno_6_samples_heatmap.pdf", width=10, height=10, useDingbats =F)
draw(hmap, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)
decorate_heatmap_body("Z-Score",{grid.rect(gp = gpar(fill = "transparent", col = "black", lwd = 2))})
dev.off()
### Done                                                                          ####
cat("Done 5\n")

#----------------------------------------------------------------------------------------
#  Visualization using ChIPseeker                                                    ####
#----------------------------------------------------------------------------------------
#

### Density of read count frequency to see where binding is relative to the TSS or each sample  ####
experiment = "bZIP-DAPseq" ### change name according to data
date <- format(Sys.time(), "%Y%m%d")#

## Generate annotation object
gffFile <- "/Users/giovanna/tomato_bZIP_sebastien_soyk/data/genome/S100_v2.1.0/SollycSweet-100_genes_v2.1.1.new.gff3"
txdb <- makeTxDbFromGFF(file=gffFile,
                        dataSource="Alonge et al., 2021",
                        organism="Solanum lycopersicum", dbxrefTag="Alias")

library(org.SlycopersicumNatalia.eg.db)
org.SlycopersicumNatalia.eg.db

# Load peak files                                                                                ####
#
getwd()
list.files("results/visualization/beds", pattern= "Input.*fdr_0.01.*.narrowPeak", full.names=T)
list.files("results/visualization/beds", pattern= "Input.*fdr_0.01.*.results.bed", full.names=T)
samplefiles <- as.list(list.files("results/visualization/beds", pattern= "Input.*fdr_0.01.*.narrowPeak", full.names=T))
samplefiles <- as.list(list.files("results/visualization/beds", pattern= "Input.*fdr_0.01.*.results.bed", full.names=T))
samplefiles
names(samplefiles) <- c("SlycSSP2", "SpimSSP2", "SSP")
samplefiles

SlycSSP2 <- readPeakFile(samplefiles$SlycSSP2)
SpimSSP2 <- readPeakFile(samplefiles$SpimSSP2)
SSP <- readPeakFile(samplefiles$SSP)
SSP
covplot(SSP, weightCol="V5")
covplot(SlycSSP2, weightCol="V5")
covplot(SpimSSP2, weightCol="V5")

length(SlycSSP2)
length(SSP)
length(SpimSSP2)

# Prepare the promoter regions
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
length(promoter)

### Calculate the tag matrix
# getTagMatrix() function can produce the matrix for visualization. peak stands for the peak file.
# window stands for a collection of regions that users want to look into.
# weightCol refers to column in peak file. This column acts as a weight value. Def=NULL
tagMatrixList <- lapply(samplefiles, getTagMatrix, windows=promoter, weightCol="V5")
tagMatrixList <- lapply(samplefiles, getTagMatrix, windows=promoter)

dim(tagMatrixList$SlycSSP2)
#[1] 3393 6001
dim(tagMatrixList$SSP)
#[1] 16022  6001
dim(tagMatrixList$SpimSSP2)
#[1] 9663 6001

###  Profile plots
###  Average Profile of ChIP peaks binding to TSS region in all samples
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95, resample=500, facet="row")

plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf = 0.95, resample = 500, facet = "row",
            xlab="Genomic Region (5'->3')", ylab = "Peak Count Frequency", ) +
    theme(axis.text=element_text(size=22),
           axis.title=element_text(size=24,face="bold")) +
    scale_color_manual(values = c(SSP="darkred",  SlycSSP2="darkgreen", SpimSSP2="orange3"))
ggsave(filename = "results/visualization/AvgProf_peaks.pdf", device = "pdf", width = 10, height = 10, dpi = 150, useDingbats=FALSE)

tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL) +
    theme(axis.text=element_text(size=22),
          axis.title=element_text(size=24,face="bold"))

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrix <- getTagMatrix(SSP, windows=promoter)
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 500,
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

tagMatrix <- getTagMatrix(SlycSSP2, windows=promoter)
plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 500,
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")

# Peak heatmaps
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)
ggsave((file.path(paste0("./ChIPseeker_",date,"/tagHeatmap_peaks_q5.pdf"))), width = 3, height = 3, useDingbats=FALSE)

### Profile of ChIP peaks binding to body regions

plotPeakProf2(peak = SpimSSP2, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 800,
              TxDb = txdb, weightCol = "V5",ignore_strand = F) +
    theme(axis.text=element_text(size=22),
          axis.title=element_text(size=24,face="bold")) +
    scale_color_manual(values = c(SSP="darkred",  SlycSSP2="darkgreen", SpimSSP2="orange3"))
ggsave((file.path(paste0("results/visualization/PeakProf2_conf_peaks_SSP.pdf"))), width = 4, height = 3, useDingbats=FALSE)
ggsave((file.path(paste0("results/visualization/PeakProf2_conf_peaks_SlycSSP2.pdf"))), width = 4, height = 3, useDingbats=FALSE)
ggsave((file.path(paste0("results/visualization/PeakProf2_conf_peaks_SpimSSP2.pdf"))), width = 4, height = 3, useDingbats=FALSE)

### For all samples
#
plotPeakProf2(files, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body",
              TxDb = txdb, facet = "row", nbin = 800) +
    theme(axis.text=element_text(size=22),
          axis.title=element_text(size=24,face="bold")) +
    scale_color_manual(values = c(SSP="darkred",  SlycSSP2="darkgreen", SpimSSP2="orange3"))
ggsave((file.path(paste0("results/visualization/PeakProf2_conf_peaks_2.pdf"))), width = 4, height = 3, useDingbats=FALSE)

# Get annotations
peakAnnoList <- lapply(samplefiles, annotatePeak, TxDb=txdb,
                       tssRegion=c(-5000, 5000), verbose=FALSE)

peakAnnoList


#### Visualize annotation distribution                                               ####

plotAnnoBar(peakAnnoList, title="Feature distribution of transcription factor-binding loci") +
    theme(axis.text=element_text(size=22),
          axis.title=element_text(size=24,face="bold"),
          legend.text=element_text(size=22),
          legend.title=element_text(size=24,face="bold"))

ggsave((file.path(paste0("results/visualization/AnnoBar_peaks_all_samples.pdf"))), width = 10, height = 7, useDingbats=FALSE)

upsetplot(peakAnno, vennpie=TRUE)

# Visualize distribution of TF-binding loci relative to TSS                          ####

plotDistToTSS(peakAnnoList, title="Distribution of TF binding loci relative to TSS")  +
    theme(axis.text=element_text(size=22),
          axis.title=element_text(size=24,face="bold"),
          legend.text=element_text(size=22),
          legend.title=element_text(size=24,face="bold"))

ggsave((file.path(paste0("results/visualization/DistToTSS_peaks_all_samples.pdf"))), width = 11, height = 7, useDingbats=FALSE)

SSP
SlycSSP2
SpimSSP2

#http://bioconductor.org/packages/devel/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html
### Profile of ChIP peaks binding to different regions                                       ####
#
# getBioRegion can prepare the different regions for ChIP peaks to bind.
# getTagMatrix can accept type, by, upstream and downstream parameters to get tagmatrix according
# to different needs.
# plotPeakProf and plotPeakProf2 supports the plotting of profiles of peaks binding to different regions.
# Users can also create heatmap or average profile of ChIP peaks binding to these regions.
#
# In order to plot body regions, a new method binning,was introduced to getTagMatrix
# binning scaled the regions having different lengths to the equal length by deviding the regions into
# the same amounts of boxes.
# Because the amount of boxes is equal, the regions can be thought of scaling to equal length.
# binning method can speed up the getTagMatrix by changing the precision from bp to box(several bps).

## The results of binning method and normal method are nearly the same.
tagMatrix_binning <- getTagMatrix(peak = SSP, TxDb = txdb,
                                  upstream = 3000, downstream = 3000,
                                  type = "start_site", by = "gene",
                                  weightCol = "V5", nbin = 800)

plotAvgProf.binning(tagMatrix_binning, conf = 0.95, resample = 500,
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
# could not find function "plotAvgProf.binning"

### Here uses `plotPeakProf2` to do all things in one step.
## Gene body regions having lengths smaller than nbin will be filtered
## A message will be given to warning users about that.
## >> 9 peaks(0.872093%), having lengths smaller than 800bp, are filtered...

## The ignore_strand is FALSE in default. We put here to emphasize that.
## We will not show it again in the below example
##
SSP
plotPeakProf2(peak = SSP, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 800,
              TxDb = txdb, weightCol = "V5",ignore_strand = F)

genebody <- getBioRegion(TxDb = txdb,
                         by = "gene",
                         type = "body")

matrix_no_flankextension <- getTagMatrix(SSP, windows = genebody, nbin = 800)

plotPeakProf(matrix_no_flankextension,conf = 0.95)

matrix_actual_extension <- getTagMatrix(SSP, windows = genebody, nbin = 800,
                                        upstream = 1000,downstream = 1000)
plotPeakProf(matrix_actual_extension,conf = 0.95)
### Done                                                                             ####
cat("Done 6\n")

#
#----------------------------------------------------------------------------------------
# Merge results from FIMO analysis                                                   ####
#----------------------------------------------------------------------------------------
#
getwd()
#m4$Annotation <- anno$overlap
peak_tab <- read.table(file="analysis/peak_table_fdr_0.01.tsv",sep="\t",stringsAsFactors=FALSE,header=TRUE)
head(peak_tab)

fimo_tab <- read.table(file="results/fimo/results/summary_table_final.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)
head(fimo_tab)

peak_tab$id <- paste(sep = "_", peak_tab$chr, peak_tab$start, peak_tab$end)
fimo_tab$id <- paste(sep = "_", fimo_tab$chr, fimo_tab$start, fimo_tab$end)

### Merge
#
dim(peak_tab)
#[1] 49170    24
dim(fimo_tab)
#[1] 428   9
head(fimo_tab)

head(fimo_tab[, -c(1:3)])
fimo_tab <- fimo_tab[, -c(1:3)]
dim(fimo_tab)
#[1] 428   6
fimo_tab$id %in% peak_tab$id
# TRUE...

dim(merge(peak_tab, fimo_tab, all=TRUE, by="id"))
# [1] 49170    29

m <- merge(peak_tab, fimo_tab, all=TRUE, by="id")
rownames(m)

head(m)
head(peak_tab)
### It doesn't keep the order of the original data frame (peak_tab)

### Keep Original Row Order when Merging Data Using join Functions of dplyr Package
#
#mo <- data_join <- inner_join(peak_tab, fimo_tab, keep = TRUE)
# joining, by = "id"
mo <- data_join <- full_join(peak_tab, fimo_tab, by="id")
dim(mo)
# [1] 49170    29
#
head(mo)

write.table(mo,sep="\t",quote=FALSE,row.names=FALSE,file="analysis/peak_table_fdr_0.01_FIMO.tsv")
write.table(mo,sep="\t",quote=FALSE,row.names=FALSE,file="results/peak_table_fdr_0.01_FIMO.tsv")

#### Save m4 as xlsx                                   ####
write_xlsx(mo,path="analysis/peak_table_fdr_0.01_FIMO.xlsx")
write_xlsx(mo,path="results/peak_table_fdr_0.01_FIMO.xlsx")

### Done                                                                          ####
cat("Done 7\n")