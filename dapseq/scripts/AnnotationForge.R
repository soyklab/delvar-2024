# AnnotationForge
# build custom gene ontology databases
# https://bioconductor.org/packages/release/bioc/html/AnnotationForge.html

setwd("./")
home <- getwd()
# BiocManager::install("AnnotationForge")

# Load libraries
library(AnnotationForge)
library(tidyverse)

# Load GO (SL4.0)
# vib plaza:https://bioinformatics.psb.ugent.be/plaza/versions/plaza_v5_dicots/download/download
# vib plaza: https://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_05/GO/go.sly.csv.gz

godat <- read_delim("go.sly.csv", skip = 8, delim="\t") # downloaded from https://ftp.psb.ugent.be/pub/plaza/plaza_public_dicots_05/GO/go.sly.csv.gz
head(godat)
nrow(godat)

# Load gene descriptions
descdat <- read.delim(file="ITAG4.0.descriptions.shortid", header=FALSE, sep='\t')
head(descdat)

#################
# Make Sym df
Sym <- descdat[,c(1,2,7)]
head(Sym)
colnames(Sym) <- c("GID", "SYMBOL", "GENENAME")
Sym$SYMBOL <- gsub("\\..*","",Sym$SYMBOL)
length(unique(Sym$GID))
length(unique(Sym$SYMBOL))
Sym$GID <- Sym$SYMBOL
head(Sym)
tail(Sym)
nrow(Sym)
length(unique(Sym$GID))

#################
# Make Chr df
head(descdat)
Chr <- descdat[,c(1,3)]
head(Chr)
colnames(Chr) <- c("GID", "CHROMOSOME")
Chr$CHROMOSOME <- gsub("SL4.0","",Chr$CHROMOSOME)

length(unique(Chr$GID))
Chr$GID <- gsub("\\..*","",Chr$GID)
length(unique(Chr$GID))

head(Chr)
tail(Chr)

#################
# Make GO df
head(godat)
sum(duplicated(godat))

GO <- godat[,c(1,3,4)]
head(GO)
sum(duplicated(GO)) # number of duplicated rows
colnames(GO) <- c("GID","GO", "EVIDENCE")
nrow(GO)
GO <- GO %>% drop_na(GO) # remove NAs
nrow(GO)

length(unique(GO$GID))
GO$GID <- gsub("\\..*","",GO$GID)
length(unique(GO$GID))

head(GO)
tail(GO)

#################
# check files dfs
nrow(Sym)
nrow(Chr)
nrow(GO)

str(Sym)
str(Chr)
str(GO)

head(Sym)
head(Chr)
head(GO)


# check overlaps between Sym and GO dataframes
SymGO_overlap <- Sym %>%
  filter(GID %in% GO$GID)
nrow(SymGO_overlap)  # number of genes WITH annotation
nrow(Sym) - nrow(SymGO_overlap) # number of genes WITHOUT annotation

Sym_noGO <- Sym %>%
  filter(!GID %in% GO$GID)
nrow(Sym_noGO)
head(Sym_noGO) # genes WITHOUT annotation
tail(Sym_noGO) # genes WITHOUT annotation


# check overlaps between Chr and GO dataframes
ChrGO_overlap <- Chr %>%
  filter(GID %in% GO$GID)
nrow(ChrGO_overlap)
nrow(Chr)
nrow(ChrGO_overlap) # number of genes WITH annotation
nrow(Chr) - nrow(ChrGO_overlap) # number of genes WITHOUT annotation
 
Chr_noGO <- Chr %>%
  filter(!GID %in% GO$GID)
nrow(Chr_noGO)
head(Chr_noGO) # genes WITHOUT annotation
tail(Chr_noGO) # genes WITHOUT annotation

# check overlaps between GO and Sym dataframes
GOSym_overlap <- GO %>%
  filter(GID %in% Sym$GID)
nrow(GOSym_overlap) 
nrow(GO)

noGOSym_overlap <- GO %>%
  filter(!GID %in% Sym$GID)
nrow(noGOSym_overlap) # shoud be 0
nrow(GO)

noGOChr_overlap <- GO %>%
  filter(!GID %in% Chr$GID)
nrow(noGOChr_overlap) # should be 0
nrow(GO)

sum(duplicated(GO))
sum(duplicated(Sym))
sum(duplicated(Chr))

GO <- distinct(GO) # remove duplicated rows

# double-check overlaps between Sym and GO dataframes
SymGO_overlap_distinct <- Sym %>%
  filter(GID %in% GO$GID)
nrow(SymGO_overlap)  # number of genes WITH annotation before removing duplicate rows
nrow(SymGO_overlap_distinct)  # number of genes WITH annotation after removing duplicate rows

# Stats:
nrow(Sym) # 34075 (number of all genes)
nrow(SymGO_overlap) # 25750 (number of annotated genes)
nrow(GO) # 1458692 (number of GO-terms; multiple terms per gene possible)



## Call function to generate package (change package name)
makeOrgPackage(gene_info=Sym, chromosome=Chr, go=GO,
# makeOrgPackage(gene_info=Sym_dist, chromosome=Chr_dist, go=GO_dist,
               version="0.1",
               maintainer="ssoyk <sebastian.soyk@unil.ch>",
               author="ssoyk <sebastian.soyk@unil.ch>",
               outputDir = ".",
               tax_id="4081",
               genus="Solanum",
               species="lycopersicumSL40VIBS", # DO NOT USE UNDERSCORE IN NAME
               goTable="go")

## then you can call install.packages based on the return value
install.packages("./org.SlycopersicumSL40VIBS.eg.db", repos = NULL, type = "source") # = "g(genus)species" DO NOT USE UNDERSCORE IN NAME

