###################################
# This R script is used to plot the distribution of deleterious variants



###################################
# determine experiment name and date
experiment = "del_var" ### change name according to data
date <- format(Sys.time(), "%Y%m%d")

###################################
# Load libraries
suppressMessages(library(tidyverse))
suppressMessages(library("ggpubr"))
suppressMessages(library(ggplot2))
suppressMessages(library(ggthemes))
suppressMessages(library(reshape2))
suppressMessages(library(stringr))
suppressMessages(library(plyr))
suppressMessages(library(dplyr))
suppressMessages(library(viridis))

###################################
# Set directories
setwd("./")
home <- getwd()
filedir <- ("/FILE-PATH-TO_SIFT4G-OUTPUT-DIRECTORY/") # change file path to directory that contains output files from sift4g analysis

# create dated parent directory for file exports
date <- format(Sys.time(), "%Y%m%d")
dir.create(file.path(home,paste0(date,"_",experiment)), showWarnings = FALSE)

###################################
# Load raw data and process

# Load accession data
dat_acc <- read.table(file.path(filedir, "accessions"), header = FALSE, sep = '\t')
colnames(dat_acc) <- c("SAMPLE", "SPECIES")
head(dat_acc)

dat_acc %>% 
  group_by(SPECIES) %>%
  dplyr::summarise(COUNT = n())

# Load del var with geneid data
dat_del <- read.table(file.path(filedir, "82_acc_Spim0.1_filtered_SIFTannotations.xls"), header = TRUE, sep = '\t')
head(dat_del)
nrow(dat_del)

# Clean gene names
dat_del$TRANSCRIPT_ID <- gsub("mRNA.","",dat_del$TRANSCRIPT_ID)
dat_del$GENE_ID <- gsub("gene.","",dat_del$GENE_ID)
dat_del$GENE_ID <- gsub("\\..*","",dat_del$GENE_ID)
head(dat_del)
tail(dat_del)

# Add variant IDs (CHROM:POS)
dat_del_tb <- dat_del %>%
  unite("VARIANT_ID", CHROM:POS, remove = FALSE) %>%
  unite("RESIDUE_ID", c(REF_AMINO,AMINO_POS,ALT_AMINO), remove = FALSE, sep=c("")) %>%
  unite("NONSYN_ID", c(GENE_ID,RESIDUE_ID), remove = FALSE, sep=c("_")) %>%
  as_tibble()
head(dat_del_tb)



#######################
# load all genotype counts for all nonsynonymous variants (LARGE FILE)
dat_var <- read.table(file.path(filedir, "82_acc_Spim0.1_filtered_SIFTpredictions.list"), header = FALSE, sep = '\t')
colnames(dat_var) <- c("CHROM", "POS", "REF", "ALT", "SAMPLE", "GT", "SIFTINFO")
head(dat_var)

# Add variant IDs (CHROM:POS)
dat_var_tb <- dat_var %>%
  unite("VARIANT_ID", CHROM:POS, remove = FALSE) %>%
  as_tibble()

# Add species info to variant list
dat_var <- dat_var_tb %>%
  left_join(dat_acc, by=NULL, suffix = c("SAMPLE", "SPECIES"), copy = FALSE,  keep = FALSE)

# Filter for DELETERIOUS variants
dat_del <- dat_var %>%
  filter(grepl("DELETERIOUS", SIFTINFO))
nrow(dat_del)



############################################
# PLOTTING GENOME WIDE
############################################

# Global picture of CDS variants

# subset data for CDS
dat_cds <- subset(dat_del_tb, REGION =="CDS")
nrow(dat_cds)

# plot VARIANT_TYPE in cds
p1 <- dat_cds %>%
  group_by(VARIANT_TYPE) %>%
  dplyr::summarise(count = n()) %>%
  ggplot(aes(x = reorder(VARIANT_TYPE,(count)), y = count, fill = 1)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = count, hjust = +0.75)) +
  coord_flip() +
  theme_light() +
  theme(axis.title.y = element_blank(), legend.position = "none") +
  scale_fill_viridis()

# plot SIFT_PREDICTION in cds
p2 <- dat_cds %>%
  group_by(SIFT_PREDICTION) %>%
  dplyr::summarise(count = n()) %>%
  ggplot(aes(x = reorder(SIFT_PREDICTION,(count)), y = count, fill = 1)) +
  geom_bar(stat = 'identity') +
  geom_text(aes(label = count, hjust = +0.75)) +
  coord_flip() +
  theme_light() +
  theme(axis.title.y = element_blank(), legend.position = "none") +
  scale_fill_viridis()

figure <- ggarrange(p1, p2,
                    labels = c("a", "b"),
                    ncol = 1, nrow = 2)

figure

# write plot to pdf
ggsave(paste0(date, "_", experiment, "/barplot_CDS_variants.pdf"), width = 5, height = 4, useDingbats=FALSE)



############################################
# Filter for CDS variants with effect on protein sequence (NONSYNONYMOUS)

# subset data for NONSYNONYMOUS
dat_nonsyn <- subset(dat_del_tb, VARIANT_TYPE =="NONSYNONYMOUS")
nrow(dat_nonsyn)

######
### Boxplots DELETERIOUS variants

# plot only alternative DELETERIOUS variants by genotype

# filter for ALTERNATIVE DELETERIOUS variants
dat_del_alt <- dat_del %>%
  filter(GT %in% c("0/1", "1/1"))

# plot (upright)
p3 <- dat_del_alt %>%
  group_by(GT , SAMPLE, SPECIES) %>% 
  dplyr::summarise(count = n()) %>%
  ggplot(aes(x=GT, y=count, color=SPECIES)) +  
  geom_jitter(width=0.25, height = 0.0, size=0.25, color="grey") +
  geom_boxplot(alpha=0, size = 0.5) +
  facet_wrap(~factor(SPECIES, levels=c('SP', 'SLC', 'SLL'))) +
  theme_light(base_size = 7) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0,NA) +
  scale_color_manual(values=c('#d95f02','#7570b3','#1b9e77'),guide=FALSE)

p3

# write plot to pdf
ggsave(paste0(date, "_", experiment, "/boxplot_ALTGT_counts_per_accession_upright.pdf"), width = 2, height = 2, useDingbats=FALSE)

# write list to csv
altcount_list <- dat_del_alt %>%
  group_by(GT , SAMPLE, SPECIES) %>% 
  dplyr::summarise(count = n())
write.csv(altcount_list, file = (paste0(date, "_", experiment, "/boxplot_ALTGT_counts_per_accession.csv")))



######
### Barplots DELETERIOUS variants

# color-code by species

altcount_list %>%
  ggplot(aes(x= reorder(SAMPLE, count), y=count, fill=SPECIES) )+
  geom_bar(position="stack", stat="identity") +
  theme_light(base_size = 6) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top") +
  scale_fill_manual(values=c('#d95f02','#7570b3','#1b9e77'))

ggsave(paste0(date, "_", experiment, "/barplot_ALTGT_counts_per_accession_SPECIES.pdf"), width = 5, height = 5, useDingbats=FALSE)


# color-code by genotype (homo/het)
altcount_list %>%
  ggplot(aes(fill=GT, y=count, x= reorder(SAMPLE, count))) + 
  geom_bar(position="stack", stat="identity") +
  theme_light(base_size = 6) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top") +
  scale_fill_manual(values=c('#9FAFC3','#3F5F88'))

ggsave(paste0(date, "_", experiment, "/barplot_ALTGT_counts_per_accession_GT.pdf"), width = 5, height = 5, useDingbats=FALSE)





############################################
# PLOTTING GENES OF INTEREST
############################################

#########
# Load goi lists
goi_bZIPA <- read.table(file.path("gene_list/bZIP_GroupA"), header = F)
nrow(goi_bZIPA)
colnames(goi_bZIPA) <- c("GENE_ID")

goi_CETS <- read.table(file.path("gene_list/CETS"), header = F)
nrow(goi_CETS)
colnames(goi_CETS) <- c("GENE_ID")


###
# Filter NONSYNONYMOUS variants for GOIs
dat_nonsyn_CETS <- dat_nonsyn %>%
  filter(GENE_ID %in% goi_CETS$GENE_ID)
nrow(dat_nonsyn_CETS)

dat_nonsyn_bZIPA <- dat_nonsyn %>%
  filter(GENE_ID %in% goi_bZIPA$GENE_ID)
nrow(dat_nonsyn_bZIPA)

###
# Calculate number of accessions per alt allele
head(dat_var_tb)

bZIPA_var.counts <- dat_var_tb %>%
  filter(grepl(paste(goi_bZIPA$geneid, collapse="|"), SIFTINFO)) %>%
  filter(GT %in% c("1/1", "0/1")) %>%
  group_by(VARIANT_ID) %>% 
  dplyr::summarise(count = n())
bZIPA_dat.var.counts <- dat_nonsyn_bZIPA %>%
  left_join(bZIPA_var.counts, by=NULL, suffix = c("VARIANT_ID", "VARIANT_ID"), copy = FALSE,  keep = FALSE)
write.csv(bZIPA_dat.var.counts, file = (paste0(date, "_", experiment, "/table_var_counts_bZIPA.csv")))


CETS_var.counts <- dat_var_tb %>%
  filter(grepl(paste(goi_CETS$geneid, collapse="|"), SIFTINFO)) %>%
  filter(GT %in% c("1/1", "0/1")) %>%
  group_by(VARIANT_ID) %>% 
  dplyr::summarise(count = n())
CETS_dat.var.counts <- dat_nonsyn_CETS %>%
  left_join(CETS_var.counts, by=NULL, suffix = c("VARIANT_ID", "VARIANT_ID"), copy = FALSE,  keep = FALSE)
write.csv(CETS_dat.var.counts, file = (paste0(date, "_", experiment, "/table_var_counts_CETS.csv")))


# plot SIFT_SCORE
p1 <-  ggplot(CETS_dat.var.counts) +
  geom_point(aes(reorder(GENE_ID, count, mean), y=SIFT_SCORE, size = count, color = -SIFT_SCORE)) +
  theme_light() +
  scale_y_reverse() +
  coord_flip() +
  theme(axis.title.y = element_blank(), axis.text.x = element_text(angle = 0)) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1) +
  theme(legend.position = "bottom") +
  scale_colour_viridis() +
  ggtitle("CETS")

p2 <- ggplot(bZIPA_dat.var.counts) +
  geom_point(aes(reorder(GENE_ID, count, mean), y=SIFT_SCORE, size = count, color = -SIFT_SCORE)) +
  theme_light() +
  scale_y_reverse() +
  coord_flip() +
  theme(axis.title.y = element_blank(), axis.text.x = element_text(angle = 0)) +
  geom_hline(yintercept=0.05, linetype="dashed", color = "red", size=1) +
  theme(legend.position = "bottom") +
  scale_colour_viridis() +
  ggtitle("bZIP-A")

figure <- ggarrange(p1, p2,
                    labels = c("a", "b"),
                    ncol = 1, nrow = 3)
figure

# write plot to pdf
ggsave(paste0(date, "_", experiment, "/dotplot_siftscore_CETS_bZIPA.pdf"), width = 5, height = 11, useDingbats=FALSE)
p1
ggsave(paste0(date, "_", experiment, "/dotplot_siftscore_CETS.pdf"), width = 3, height = 2.5, useDingbats=FALSE)
p2
ggsave(paste0(date, "_", experiment, "/dotplot_siftscore_bZIPA.pdf"), width = 3, height = 3, useDingbats=FALSE)




###
# Calculate number of accessions categories per alt allele
head(dat_var_tb)
head(dat_nonsyn_CETS)


CETS_var <- dat_var_tb %>%
  filter(grepl(paste(goi_CETS$GENE_ID, collapse="|"), SIFTINFO)) %>%
  filter(grepl("NONSYNONYMOUS", SIFTINFO)) %>%
  select(!SIFTINFO)

bZIPA_var <- dat_var_tb %>%
  filter(grepl(paste(goi_bZIPA$GENE_ID, collapse="|"), SIFTINFO)) %>%
  filter(grepl("NONSYNONYMOUS", SIFTINFO)) %>%
  select(!SIFTINFO)


# Add species category info to variant list
head(CETS_var)
head(bZIPA_var)
head(dat_acc)

var_CETS_cat <- CETS_var %>%
  left_join(dat_acc, by=NULL, suffix = c("SAMPLE", "SPECIES"), copy = FALSE,  keep = FALSE) %>%
  as_tibble()
head(var_CETS_cat)

var_bZIPA_cat <- bZIPA_var %>%
  left_join(dat_acc, by=NULL, suffix = c("SAMPLE", "SPECIES"), copy = FALSE,  keep = FALSE) %>%
  as_tibble()
head(var_bZIPA_cat)



# Add gene info to variant / category list
head(var_CETS_cat)
head(var_bZIPA_cat)
head(dat_del_tb)

dat_CETS <- var_CETS_cat %>%
  left_join(dat_del_tb, by=NULL, suffix = c("NONSYN_ID", "GENE_ID"), copy = FALSE,  keep = FALSE)
head(dat_CETS)

dat_bZIPA <- var_bZIPA_cat %>%
  left_join(dat_del_tb, by=NULL, suffix = c("NONSYN_ID", "GENE_ID"), copy = FALSE,  keep = FALSE)
head(dat_bZIPA)


################################
# Plots
################################
# CETS

# all synonymous
CETS_counts <- dat_CETS %>%
  group_by(NONSYN_ID, SPECIES, GT) %>% 
  dplyr::summarise(COUNT = n())

CETS_counts %>%
  filter(GT %in% c("1/1", "0/1")) %>%
  ggplot(aes(fill=SPECIES, y=COUNT, x = NONSYN_ID)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~factor(SPECIES, levels=c('SP','SLC','SLL')), ncol=3) +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top") +
  theme_light(base_size = 6) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c('#d95f02','#7570b3','#1b9e77'))+
  ggtitle("CETS synonymous variants")

ggsave(paste0(date, "_", experiment, "/barplot_CETS_counts_per_species_nonsyn.pdf"), width = 2.5, height = 3, useDingbats=FALSE)

# write list to csv
write.csv(CETS_counts, file = (paste0(date, "_", experiment, "/barplot_CETS_counts_per_species_nonsyn.csv")))



# deleterious only
CETS_counts_del <- dat_CETS %>%
  filter(!grepl("TOLERATED", SIFT_PREDICTION)) %>%
  group_by(NONSYN_ID, SPECIES, GT) %>% 
  dplyr::summarise(COUNT = n())

CETS_counts_del %>%
  filter(GT %in% c("1/1", "0/1")) %>%
  ggplot(aes(fill=SPECIES, y=COUNT, x = NONSYN_ID)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~factor(SPECIES, levels=c('SP','SLC','SLL')), ncol=3) +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top") +
  theme_light(base_size = 6) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c('#d95f02','#7570b3','#1b9e77'))+
  ggtitle("CETS deleterious variants")

ggsave(paste0(date, "_", experiment, "/barplot_CETS_counts_per_species_del.pdf"), width = 2.5, height = 1, useDingbats=FALSE)

# write list to csv
write.csv(CETS_counts_del, file = (paste0(date, "_", experiment, "/barplot_CETS_counts_per_species_del.csv")))


################################
# bZIPs

# all synonymous
bZIPA_counts <- dat_bZIPA %>%
  group_by(NONSYN_ID, SPECIES, GT) %>% 
  dplyr::summarise(COUNT = n())

bZIPA_counts %>%
  filter(GT %in% c("1/1", "0/1")) %>%
  ggplot(aes(fill=SPECIES, y=COUNT, x = NONSYN_ID)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~factor(SPECIES, levels=c('SP','SLC','SLL')), ncol=3) +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top") +
  theme_light(base_size = 6) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c('#d95f02','#7570b3','#1b9e77'))+
  ggtitle("bZIP synonymous variants")

ggsave(paste0(date, "_", experiment, "/barplot_bZIPA_counts_per_species_nonsyn.pdf"), width = 2.5, height = 3, useDingbats=FALSE)

# write list to csv
write.csv(bZIPA_counts, file = (paste0(date, "_", experiment, "/barplot_bZIPA_counts_per_species_nonsyn.csv")))


# all deleterious
bZIPA_counts_del <- dat_bZIPA %>%
  filter(!grepl("TOLERATED", SIFT_PREDICTION)) %>%
  group_by(NONSYN_ID, SPECIES, GT) %>% 
  dplyr::summarise(COUNT = n())

bZIPA_counts_del %>%
  filter(GT %in% c("1/1", "0/1")) %>%
  ggplot(aes(fill=SPECIES, y=COUNT, x = NONSYN_ID)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~factor(SPECIES, levels=c('SP','SLC','SLL')), ncol=3) +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), legend.position="top") +
  theme_light(base_size = 6) +
  theme(legend.position = "none") +
  scale_fill_manual(values=c('#d95f02','#7570b3','#1b9e77'))+
  ggtitle("bZIP deleterious variants")

ggsave(paste0(date, "_", experiment, "/barplot_bZIPA_counts_per_species_del.pdf"), width = 2.5, height = 1.2, useDingbats=FALSE)

# write list to csv
write.csv(bZIPA_counts_del, file = (paste0(date, "_", experiment, "/barplot_bZIPA_counts_per_species_del.csv")))

