########################################
# ONLY WORKS FOR CIS INTERACTIONS!
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: supplemental table with significant z-scrore interactions (from 3Dflow output) and SNP rsIDs
#            SNP rsID for SNP to be examined
#            window size, numeric, bp, (to look for interactions within this window around the SNP)
#            gene list (from winiona)
#            bin size used for Hi-C data (numeric, bp)
########################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
SNP <- args[2]
window_size <- as.numeric(as.character(args[3]))
genes <- args[4]
bin_size <- as.numeric(as.character(args[5]))


##########
library(dplyr)
library(tidyr)
library(stringr)
#library(multcomp) # for anova and tukey test
.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2")
#install_github("vqv/ggbiplot")
##remotes::install_github("R-CoderDotCom/ridgeline@main")
#library(ridgeline)
##########

#########################################################################
#########################################################################
print("#read in files")
####interaction data
##options(scipen = 999)
#dat_file <- "CM_cis50Kb_SNPs_zscoreSig_interactions_and_SNP_rsID.txt"
##SNP <- "rs11191568"
#SNP <- "rs4631439"
#window_size <- 1000000
#genes <- "nearby_rs4631439genes.txt"
##germlayer_file <- "germlayer_info.txt"
#bin_size <- 50000

#read in interaction df
dat <- read.table(dat_file, header = TRUE)
print("summary of ALL sig zscores per cell type")
summary(dat)
#only interactions with SNPs
dat <- dat %>% filter(!is.na(SNP_rsID_A)|!is.na(SNP_rsID_B))
print("#split ID col")
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
dat <- dat %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
print("#remove A and B from chrom names")
dat$chrA <- gsub("A", "", dat$chrA)
dat$chrB <- gsub("B", "", dat$chrB)
head(dat)
#filter for inters with SNP of interest
SNP_inters <- dat %>% filter(grepl(SNP, SNP_rsID_A) | grepl(SNP,SNP_rsID_B))
SNP_inters$st1 <- as.numeric(as.character(SNP_inters$st1))
SNP_inters$end1 <- as.numeric(as.character(SNP_inters$end1))
SNP_inters$st2 <- as.numeric(as.character(SNP_inters$st2))
SNP_inters$end2 <- as.numeric(as.character(SNP_inters$end2))
summary(SNP_inters)
#get interaction bin with SNP of interest
SNP_inter_bin <- SNP_inters %>% filter(grepl(SNP, SNP_rsID_A)) %>% select(chrA,st1,end1) %>% unique()
SNP_inter_bin
wl <- as.numeric(as.character(SNP_inter_bin$st1[1])) - window_size
wr <- as.numeric(as.character(SNP_inter_bin$end1[1])) + window_size
#get inters within Xbp window of SNP, both regions have to be within window
SNP_inters_window <- SNP_inters %>% filter((st1 >= wl && st1 <= wr) && (st2 >= wl && st2 <= wr))
SNP_inters_window$st1 <- as.numeric(as.character(SNP_inters_window$st1))
SNP_inters_window$st2 <- as.numeric(as.character(SNP_inters_window$st2))
SNP_inters_window$end1 <- as.numeric(as.character(SNP_inters_window$end1))
SNP_inters_window$end2 <- as.numeric(as.character(SNP_inters_window$end2))

#read in gene list file from winona
genes_df <- read.table(genes, sep = "\t", header = T)
#bin start and end that genes are in, NOTE: they often overlap multiple bins!
genes_df$bin_st <- plyr::round_any(genes_df$start, bin_size, f = floor)
genes_df$bin_end <- plyr::round_any(genes_df$end, bin_size, f = ceiling)
#splitting up genes that cover multiple bin (repeating them for join with interactions later)
genes_df <- genes_df %>% mutate(dif = bin_end - bin_st)
mgenes_df <- data.frame(seqnamess=character(),
                        start=numeric(),
                        end=numeric(),
                        width=numeric(),
                        strand=character(),
                        gene_id=character(),
                        ensembl=character(),
                        bin_st=numeric(),
                        bin_end=numeric(),
                        diff=numeric(),
                        stringsAsFactors=FALSE)
for (r in 1:nrow(genes_df)){
  #r=2
  td <- genes_df[r,]
  #if diff is the same as the bin size
  if (td$dif[1] == bin_size){
    mgenes_df <- rbind(mgenes_df, td)
  } else {
    #bin starts that cover gene
    sseq <- seq(td$bin_st[1],td$bin_st[1]+td$dif[1]-bin_size, by=bin_size)
    #bin ends that cover gene
    eseq <- seq(td$bin_end[1]-td$dif[1]+bin_size,td$bin_end[1], by=bin_size)
    #add each bin that covers the gene to the df
    for (i in 1:length(sseq)){
      trow <- as.data.frame(t(c(td$seqnames[1],td$start[1],td$end[1],td$width[1],td$strand[1],td$gene_id[1],td$ensembl[1],sseq[i],eseq[i],td$dif[1])))
      trow$V2 <- as.numeric(as.character(trow$V2))
      trow$V3 <- as.numeric(as.character(trow$V3))
      trow$V4 <- as.numeric(as.character(trow$V4))
      trow$V8 <- as.numeric(as.character(trow$V8))
      trow$V9 <- as.numeric(as.character(trow$V9))
      trow$V10 <- as.numeric(as.character(trow$V10))
      colnames(trow) <- colnames(genes_df)
      mgenes_df <- rbind(mgenes_df, trow)
    }#for
  }#if else
}#for each nearby gene
mgenes_df <- mgenes_df %>% select(-dif)
head(mgenes_df)

#list all gene IDs in one column
genes_id_start <- mgenes_df %>% 
  select(bin_st,gene_id) %>%
  group_by(bin_st) %>% 
  dplyr::mutate(gene_id_st = paste0(gene_id, collapse = ",")) %>%
  select(-gene_id) %>%
  unique()
genes_id_end <- mgenes_df %>% 
  select(bin_end,gene_id) %>%
  group_by(bin_end) %>% 
  dplyr::mutate(gene_id_end = paste0(gene_id, collapse = ",")) %>%
  select(-gene_id) %>%
  unique()
ensembl_start <- mgenes_df %>% 
  select(bin_st,ensembl) %>%
  group_by(bin_st) %>% 
  dplyr::mutate(ensembl_st = paste0(ensembl, collapse = ",")) %>%
  select(-ensembl) %>%
  unique()
ensembl_end <- mgenes_df %>% 
  select(bin_end,ensembl) %>%
  group_by(bin_end) %>% 
  dplyr::mutate(ensembl_end = paste0(ensembl, collapse = ",")) %>%
  select(-ensembl) %>%
  unique()
#merge above dfs
genes_comb_st <- full_join(genes_id_start, ensembl_start, by = "bin_st")
genes_comb_end <- full_join(genes_id_end, ensembl_end, by = "bin_end")
#add gene names to interactions df
#add genes in anchor interaction (chrA)
colnames(genes_comb_st) <- c("st1","gene_id_st1","ensembl_st1")
SNP_inters_window <- full_join(SNP_inters_window, genes_comb_st, by = "st1")
colnames(genes_comb_end) <- c("end1","gene_id_end1","ensembl_end1")
SNP_inters_window <- full_join(SNP_inters_window, genes_comb_end, by = "end1")
#add genes in anchor interaction (chrA)
colnames(genes_comb_st) <- c("st2","gene_id_st2","ensembl_st2")
SNP_inters_window <- full_join(SNP_inters_window, genes_comb_st, by = "st2")
colnames(genes_comb_end) <- c("end2","gene_id_end2","ensembl_end2")
SNP_inters_window <- full_join(SNP_inters_window, genes_comb_end, by = "end2")
#combin redundant cols
SNP_inters_window <- SNP_inters_window %>% mutate(gene_idA = paste0(gene_id_st1,",",gene_id_end1)) %>% select(-gene_id_st1,-gene_id_end1)
SNP_inters_window <- SNP_inters_window %>% mutate(gene_idB = paste0(gene_id_st2,",",gene_id_end2)) %>% select(-gene_id_st2,-gene_id_end2)
SNP_inters_window <- SNP_inters_window %>% mutate(ensemblA = paste0(ensembl_st1,",",ensembl_end1)) %>% select(-ensembl_st1,-ensembl_end1)
SNP_inters_window <- SNP_inters_window %>% mutate(ensemblB = paste0(ensembl_st2,",",ensembl_end2)) %>% select(-ensembl_st2,-ensembl_end2)
head(SNP_inters_window)
tail(SNP_inters_window)

#rank/order interactions based on zscore (highest to lowest)
#write to file per cell type
cells <- colnames(SNP_inters_window)
cells <- cells[8:(ncol(SNP_inters_window)-6)]
cells
for (c in cells){
  #c = "H9hESC_day00_Zhang"
  nc <- ncol(SNP_inters_window)
  #select only cell of interest col
  tmpd <- SNP_inters_window[,c(1,which(colnames(SNP_inters_window) %in% c), nc-5, nc-4, nc-3,nc-2,nc-1,nc)]
  #make cell col numeric
  tmpd[,2] <- as.numeric(as.character(tmpd[,2]))
  summary(tmpd)
  tmpd <- tmpd[order(tmpd[,2], decreasing = TRUE),]
  #write to file
  ws <- window_size /1000000
  write.table(tmpd, file = as.character(paste0(SNP,"_",ws,"Mb_window_",c,"_ranked_interactions_SNP_rsID_genes.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = T)
}#for
print("DONE")
#
