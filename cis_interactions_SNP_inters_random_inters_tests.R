########################################
#validating known trans-chromosomal interactions
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: reg SNPs to be tested filtered interactions file (tab separated, first 6 columns are a bed-like format, last few columns are other info)
#            reg SNPs (all) filtered interactions file (tab separated, first 6 columns are a bed-like format, last few columns are other info)
#            non-reg SNPs filtered interactions file (tab separated, first 6 columns are a bed-like format, last few columns are other info)
#            z-scores, 3Dflow output, ALL sig interactions (tab separated,)
#            cell order for graphs (file with name of cells as it matches the data frame, one cell per line)
#            out file prefix
#            number of permutations to perform (numeric)
#            bin size (bp)
#            original file with SNP info (tsv) ONLY needed when the test inters is inters with at least one SNP
########################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
testSNPs_intersfile <- args[1]
regSNPs_intersfile <- args[2]
nonRegSNPs_intersfile <- args[3]
all_intersfile <- args[4]
cellsfile <- args[5]
outfile <- args[6]
permnum <- as.numeric(as.character(args[7]))
bin_size <- as.numeric(as.character(args[8]))
SNP_file <- args[9]
print("# SNP_file")
SNP_file

##########
library(tidyr)
library(purrr)
library(forcats)
library(dplyr)
library(ggplot2)
#library(nortest, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.3.6.1") #for normality test with large sample size
#library(hexbin)
library(ggforce)#for ridgeline
library(ggridges)#for ridgeline
.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.3.6.1")
##########

#########################################################################
#set graph theme
theme_set(theme_bw() + theme(strip.background =element_rect(fill="#e7e5e2")) +
            #change size of facet header text
            theme(strip.text = element_text(size =10.49)) +
            theme(plot.title = element_text(hjust = 0.5, size = 18),
                  panel.background = element_rect(fill = "white", colour = NA),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.spacing = unit(0.25, "lines"),
                  axis.text=element_text(size=18),
                  axis.title = element_text(size = 18),
                  #plot margins
                  #plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
                  #for second legend on y-axis
                  axis.text.y.right = element_text(size=18),
                 # legend.title = element_blank(),
                  legend.text = element_text(size = 18),
                  #change the colour of facet label background
                  strip.background = element_rect(fill = "#E6E1EA"),
                  #remove space between facest
                  panel.spacing.x=unit(0, "lines"),
                  #                  legend.key = element_blank(),
                  legend.background=element_blank(),
                  #legend background
                  legend.key = element_rect(fill = NA),
                  #                  legend.position="none")
                  legend.position="top")
)
#########################################################################


print("#read in files")
##interaction data
#testSNPs_intersfile <- "VSMC_cis50Kb_SNPs_in_both_inter_regions_overlapping_intearctions_for_circos.txt"
#regSNPs_intersfile <- "VSMC_cis50Kb_SNPs_in_both_inter_regions_overlapping_intearctions_for_circos.txt"
#nonRegSNPs_intersfile <- "nonRegSNPs_all_cells_cis50Kb_SNPs_MERGED_overlapping_intearctions_for_circos.txt"
#all_intersfile <- "test_3Dflow_zscores_cis_noCrossValid2.txt"
#bin_size <- 50000
#outfile <- "test_cis_SNP_blood_pressure"
#cellsfile <- "cell_subset2.txt"
#permnum <- 50
#SNP_file <- "VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt"
#library(harrypotter)
#library(factoextra)
#library(hexbin)
#library(ggVennDiagram)
#library(ggupset)

cells_sub <- read.table(cellsfile)
cells_sub

#read test SNP inters files (to be tested)
testSNP_inters <- read.table(testSNPs_intersfile, header = F)
colnames(testSNP_inters) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
#filter by desired cells
testSNP_inters <- testSNP_inters %>% filter(cell %in% cells_sub$V1)
testSNP_inters <- testSNP_inters %>% mutate(ID = paste0("A",chrA,".",st1,".",end1,".B",chrB,".",st2,".",end2))
head(testSNP_inters)
#read reg inters files (all)
regSNP_inters <- read.table(regSNPs_intersfile, header = F)
colnames(regSNP_inters) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
#filter by desired cells
regSNP_inters <- regSNP_inters %>% filter(cell %in% cells_sub$V1)
regSNP_inters <- regSNP_inters %>% mutate(ID = paste0("A",chrA,".",st1,".",end1,".B",chrB,".",st2,".",end2))
head(regSNP_inters)
#read non-reg inters files
nonRegSNP_inters <- read.table(nonRegSNPs_intersfile, header = F)
colnames(nonRegSNP_inters) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
#filter by desired cells
nonRegSNP_inters <- nonRegSNP_inters %>% filter(cell %in% cells_sub$V1)
nonRegSNP_inters <- nonRegSNP_inters %>% mutate(ID = paste0("A",chrA,".",st1,".",end1,".B",chrB,".",st2,".",end2))
head(nonRegSNP_inters)


#num inters per bin
print("#counting each interaction twice (once for each chrom in interaction)")
anchD <- testSNP_inters
anchD$AllChr <- anchD$chrA
anchD$AllSt <- anchD$st1
anchD$AllEnd <- anchD$end1
tarD <- testSNP_inters
tarD$AllChr <- tarD$chrB
tarD$AllSt <- tarD$st2
tarD$AllEnd <- tarD$end2
testSNP_inters_bin <- rbind(anchD,tarD) %>% select(cell,AllChr,AllSt, zscore)

#separating pos and neg zscores
testSNP_inters_bin_pos <- testSNP_inters_bin %>%
  group_by(AllChr,AllSt,cell) %>%
  filter(zscore >0) %>%
  summarise(mPosZscore = mean(zscore), .groups = "keep")
head(testSNP_inters_bin_pos)
testSNP_inters_bin_neg <- testSNP_inters_bin %>%
  group_by(AllChr,AllSt,cell) %>%
  filter(zscore <0) %>%
  summarise(mNegZscore = mean(zscore), .groups = "keep")
head(testSNP_inters_bin_neg)
testSNP_inters_bin_all <- full_join(testSNP_inters_bin_pos, testSNP_inters_bin_neg, by = c("AllChr","AllSt","cell"))
head(testSNP_inters_bin_all)

testSNP_inters_count <- rbind(anchD,tarD) %>% select(cell,AllChr,AllSt) %>%
  group_by(AllChr,AllSt,cell) %>%
  summarise(numInters = n(), .groups = "keep")
head(testSNP_inters_count)

#read in all inters 3Dflow output file
ALL_inters <- read.table(all_intersfile, header = T, sep = "\t")
#split ID col
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
ALL_inters$ID <- sub("B", "\\.B", as.character(ALL_inters$ID))
ALL_inters <- ALL_inters %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
#remove A and B from chrom names
ALL_inters$chrA <- gsub("A", "", ALL_inters$chrA)
ALL_inters$chrB <- gsub("B", "", ALL_inters$chrB)
#wide to long
ALL_inters <- ALL_inters %>% gather(key = "cell", value = "zscore", 8:ncol(ALL_inters))
head(ALL_inters)

#random interactions pool/universe
#filtered so it does not contain reg inters (anchor and target)
inters_univ <- ALL_inters %>% filter(!ID %in% unique(regSNP_inters$ID)) %>%
  #filtered so it does not contain non-reg inters (anchor and target)
  filter(!ID %in% unique(nonRegSNP_inters$ID)) %>%
  #filtered so it does not contain non-reg inters (anchor and target)
  filter(!ID %in% unique(testSNP_inters$ID)) %>%
  #mutate to create cols for anchor and target IDs
  mutate(idA = paste0(chrA, ".", st1,".",end1))%>%
  mutate(idB = paste0(chrB, ".", st2,".",end2))
head(inters_univ)
bin_inters_univ <- unique(c(inters_univ$idA, inters_univ$idB))

print("#list of bin involved in interactions")
testSNP_inters <- testSNP_inters %>% mutate(idA = paste0(chrA, ".", st1,".",end1)) %>%
  mutate(idB = paste0(chrB, ".", st2,".",end2))
list_SNP_bins <- c()
if (!is.na(SNP_file)){ #at least one inter region has a SNP
  #read in SNP info
  SNP <- read.table(SNP_file, header = TRUE,sep = "\t")
  #accounting for merged SNP files
  SNP_df <- SNP
  if ("logFC_comp" %in% colnames(SNP)){
    SNP_df <- SNP %>% dplyr::select(chr_b38, start_b38,end_b38, logFC_comp)
  } else {
    SNP_df$logFC_comp.x <- as.numeric(as.character(SNP_df$logFC_comp.x))
    SNP_df$logFC_comp.y <- as.numeric(as.character(SNP_df$logFC_comp.y))
    SNP_df$logFC_comp <- mean(c(abs(SNP_df$logFC_comp.x),abs(SNP_df$logFC_comp.y)))
  }#if else
  #bins that SNPs are in
  SNP_df$bin <- plyr::round_any(SNP_df$start_b38, as.numeric(as.character(bin_size)), f = floor)
  colnames(SNP_df) <- c("AllChr","SNPstart","SNPend","logFC","AllSt")
  SNP_ID_df <- SNP_df
  SNP_ID_df$ID <- format(SNP_ID_df$ID, scientific = FALSE)
  SNP_ID_df$ID <- paste0(SNP_ID_df$AllChr,".",format(SNP_ID_df$AllSt, scientific = FALSE),".",format(SNP_ID_df$AllSt + as.numeric(as.character(bin_size)),scientific = FALSE))
  SNP_ID_df$ID <- gsub(" ","", SNP_ID_df$ID)
  #decrease SNP list to just unique IDs, find mean of |logFC|
  SNP_uniq <- SNP_ID_df %>%
    group_by(ID) %>%
    dplyr::summarise(mlogFC = mean(abs(logFC)))
  list_SNP_bins <- unique(SNP_uniq$ID)
} else { #dealing with input that is both inter regions have SNPs
  list_SNP_bins <- unique(c(list_SNP_bins,testSNP_inters$idA, testSNP_inters$idB))
}
list_SNP_bins

#functions for permutations
#mean of positive zscores
permdat_poszscore <- function(x,y,z) {
  rbins <- sample(x,y, replace = F)
  rinters <- z %>% filter(cell == c) %>% filter(idA %in% rbins | idB %in% rbins) %>% na.omit() %>%
    filter(zscore > 0)
  mean(rinters$zscore)
}
permdat_negzscore <- function(x,y,z) {
  rbins <- sample(x,y, replace = F)
  rinters <- z %>% filter(cell == c) %>% filter(idA %in% rbins | idB %in% rbins) %>% na.omit() %>%
    filter(zscore < 0)
  mean(rinters$zscore)
}
permdat_numInters <- function(x,y,z) {
  rbins <- sample(x,y, replace = F)
  rinters <- z %>% filter(cell == c) %>% filter(idA %in% rbins | idB %in% rbins) %>% na.omit()
  #calculate values per bin for random inters
  anchD <- rinters
  anchD$AllChr <- anchD$chrA
  anchD$AllSt <- anchD$st1
  anchD$AllEnd <- anchD$end1
  tarD <- rinters
  tarD$AllChr <- tarD$chrB
  tarD$AllSt <- tarD$st2
  tarD$AllEnd <- tarD$end2
  rinters_bin <- rbind(anchD,tarD) %>% select(cell,AllChr,AllSt, zscore)
  rinters_count <- rbind(anchD,tarD) %>% select(cell,AllChr,AllSt) %>%
    group_by(AllChr,AllSt,cell) %>%
    summarise(numInters = n(), .groups = "keep")
  mean(rinters_count$numInters)
}


#test between random inters and SNP inters
#df to save pvals
perm_pvals <- data.frame(cell=character(),
                          pval=integer(),
                          testType=character(),
                          stringsAsFactors=FALSE)
for (c in cells_sub$V1){
  #c = "OmniC_pooled_M_0d"
  print("############")
  print(c)
  print("############")
  ##num inters in SNPs df
  #numInters <- testSNP_inters %>% filter(cell == c) %>% nrow()
  print("#num bins with SNPs in test reg inters")
  numBins <- sum(unique(c(testSNP_inters$idA, testSNP_inters$idB)) %in% list_SNP_bins)
  print(numBins)
  #SNPs info
  tSNP_inters_count <- testSNP_inters_count %>% filter(cell == c)
  tSNP_inters_bin_all <- testSNP_inters_bin_all %>% filter(cell == c)
  
  set.seed(369)
  #get same num of inters from univ without replacement
  hisdat_posZ <- as.data.frame(replicate(permnum,permdat_poszscore(bin_inters_univ, numBins,inters_univ)))
  colnames(hisdat_posZ) <- c("posZ")
  p <- (ggplot(hisdat_posZ,aes(x=posZ))
        + geom_histogram()
        + geom_vline(aes(xintercept = mean(tSNP_inters_bin_all$mPosZscore, na.rm = T)), colour = "blue")
        + labs(y = "Frequency",
               x = "Mean positive z-score per interaction",
               title = c)
        + scale_x_continuous(expand = c(0, 0))
        + scale_y_continuous(expand = c(0, 0))
  )
  pdf(paste0(outfile,"_",c,"_histogram_positive_zscore.pdf"), width = 14, height = 4)
  print(p)
  dev.off()
  # p-value for positive zscore
  pospval <- ((hisdat_posZ %>% filter(posZ >= mean(tSNP_inters_bin_all$mPosZscore, na.rm = T)) %>% nrow())+1) / (permnum+1)
  #add to table
  tr <- c(c,pospval,"posZscore")
  perm_pvals[nrow(perm_pvals) + 1, ] <- tr
  hisdat_negZ <- as.data.frame(replicate(permnum,permdat_negzscore(bin_inters_univ, numBins,inters_univ)))
  colnames(hisdat_negZ) <- c("negZ")
  p <- (ggplot(hisdat_negZ,aes(x=negZ))
        + geom_histogram()
        + geom_vline(aes(xintercept = mean(tSNP_inters_bin_all$mNegZscore, na.rm = T)), colour = "blue")
        + labs(y = "Frequency",
               x = "Mean negative z-score per interaction",
               title = c)
        + scale_x_continuous(expand = c(0, 0))
        + scale_y_continuous(expand = c(0, 0))
  )
  pdf(paste0(outfile,"_",c,"_histogram_negative_zscore.pdf"), width = 14, height = 4)
  print(p)
  dev.off()
  # p-value for negative zscore
  negpval <- ((hisdat_negZ %>% filter(negZ <= mean(tSNP_inters_bin_all$mNegZscore, na.rm = T)) %>% nrow())+1) / (permnum+1)
  #add to table
  tr <- c(c,negpval,"negZscore")
  perm_pvals[nrow(perm_pvals) + 1, ] <- tr
  hisdat_num <- as.data.frame(replicate(permnum,permdat_numInters(bin_inters_univ,numBins,inters_univ)))
  colnames(hisdat_num) <- c("numInters")
  p <- (ggplot(hisdat_num,aes(x=numInters))
        + geom_histogram()
        + geom_vline(aes(xintercept = mean(tSNP_inters_count$numInters)), colour = "blue")
        + labs(y = "Frequency",
               x = "Mean # of interactions per bin",
               title = c)
        + scale_x_continuous(expand = c(0, 0))
        + scale_y_continuous(expand = c(0, 0))
  )
  pdf(paste0(outfile,"_",c,"_histogram_num_inters.pdf"), width = 14, height = 4)
  print(p)
  dev.off()
  # p-value for number of inters per bin
  numpval <- ((hisdat_num %>% filter(numInters >= mean(tSNP_inters_count$numInters, na.rm = T)) %>% nrow())+1) / (permnum+1)
  #add to table
  tr <- c(c,numpval,"numInters")
  perm_pvals[nrow(perm_pvals) + 1, ] <- tr
}
perm_pvals

#write permutation pvalues to df
write.table(perm_pvals, file = as.character(paste0(outfile,"_permutation_pvalues.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = T)




##test between random inters and SNP inters
#for (c in cells_sub$V1){
#  c = "OmniC_pooled_M_0d"
#  print("############")
#  print(c)
#  print("############")
#  #num inters in SNPs df
#  numInters <- testSNP_inters %>% filter(cell == c) %>% nrow()
#  #SNPs info
#  tSNP_inters_count <- testSNP_inters_count %>% filter(cell == c)
#  tSNP_inters_bin_all <- testSNP_inters_bin_all %>% filter(cell == c)
#  
#  set.seed(369)
#  #get same num of inters from univ without replacement
#  rinters <- sample_n(inters_univ %>% filter(cell == c) %>% na.omit(), numInters, replace = F)
#  rinters
#  print("#write random inters to df")
#  write.table(rinters %>% select(-ID), file = as.character(paste0(outfile,"_",c,"_random_intearctions_for_circos.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#  #calculate values per bin for random inters
#  anchD <- rinters
#  anchD$AllChr <- anchD$chrA
#  anchD$AllSt <- anchD$st1
#  anchD$AllEnd <- anchD$end1
#  tarD <- rinters
#  tarD$AllChr <- tarD$chrB
#  tarD$AllSt <- tarD$st2
#  tarD$AllEnd <- tarD$end2
#  rinters_bin <- rbind(anchD,tarD) %>% select(cell,AllChr,AllSt, zscore)
#  #separating pos and neg zscores
#  rinters_bin_pos <- rinters_bin %>%
#    group_by(AllChr,AllSt,cell) %>%
#    filter(zscore >0) %>%
#    summarise(mPosZscore = mean(zscore), .groups = "keep")
#  rinters_bin_neg <- rinters_bin %>%
#    group_by(AllChr,AllSt,cell) %>%
#    filter(zscore <0) %>%
#    summarise(mNegZscore = mean(zscore), .groups = "keep")
#  rinters_bin_all <- full_join(rinters_bin_pos, rinters_bin_neg, by = c("AllChr","AllSt","cell"))
#  rinters_count <- rbind(anchD,tarD) %>% select(cell,AllChr,AllSt) %>%
#    group_by(AllChr,AllSt,cell) %>%
#    summarise(numInters = n(), .groups = "keep")
#  
#  print("#Mann-Whitney Test between num inters of SNP inters and random inters")
#  print(wilcox.test(tSNP_inters_count$numInters, rinters_count$numInters))
#  print("summary: SNP inters")
#  print(summary(tSNP_inters_count))
#  print("summary: random inters")
#  print(summary(rinters_count))
#  
#  print("summary: SNP inters zscore")
#  print(summary(tSNP_inters_bin_all))
#  print("summary: random inters zscore")
#  print(summary(rinters_bin_all))
#  print("#Mann-Whitney Test between POSITIVE z-score of SNP inters and random inters")
#  print(wilcox.test(tSNP_inters_bin_all$mPosZscore, rinters_bin_all$mPosZscore))
#  print("#Mann-Whitney Test between NEGATIVE z-score of SNP inters and random inters")
#  print(wilcox.test(tSNP_inters_bin_all$mNegZscore, rinters_bin_all$mNegZscore))
#}
print("DONE")
