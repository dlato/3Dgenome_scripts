########################################
#validating known trans-chromosomal interactions
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: reg SNPs filtered interactions file (tab separated, first 6 columns are a bed-like format, last few columns are other info)
#            z-scores, 3Dflow output, ALL sig interactions (tab separated,)
#            cell order for graphs (file with name of cells as it matches the data frame, one cell per line)
########################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
regSNPs_intersfile <- args[1]
all_intersfile <- args[2]
cellsfile <- args[3]

##########
library(tidyr)
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
#regSNPs_intersfile <- "CardioMyo_cis50Kb_filtered_inters_SNP_CFDP1_overlapping_intearctions_for_circos.txt"
#all_intersfile <- "test_3Dflow_zscores_cis_noCrossValid2.txt"
#bin_size <- 50000
#outprefix <- "test_cis_SNP_blood_pressure"
#cellsfile <- "cell_subset_CM.txt"
#library(harrypotter)
#library(factoextra)
#library(hexbin)
#library(ggVennDiagram)
#library(ggupset)

cells_sub <- read.table(cellsfile)
cells_sub

#read inters files
SNP_inters <- read.table(regSNPs_intersfile, header = F)
colnames(SNP_inters) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
#filter by desired cells
SNP_inters <- SNP_inters %>% filter(cell %in% cells_sub$V1)
SNP_inters <- SNP_inters %>% mutate(ID = paste0("A",chrA,".",st1,".",end1,".B",chrB,".",st2,".",end2))
head(SNP_inters)
#num inters per bin
print("#counting each interaction twice (once for each chrom in interaction)")
anchD <- SNP_inters
anchD$AllChr <- anchD$chrA
anchD$AllSt <- anchD$st1
anchD$AllEnd <- anchD$end1
tarD <- SNP_inters
tarD$AllChr <- tarD$chrB
tarD$AllSt <- tarD$st2
tarD$AllEnd <- tarD$end2
SNP_inters_bin <- rbind(anchD,tarD) %>% select(cell,AllChr,AllSt, zscore)
print("TEST SUMMARY")
summary(SNP_inters_bin)
SNP_inters_bin %>% filter(cell == "OmniC_pooled_M_0d") %>% summary()

#separating pos and neg zscores
SNP_inters_bin_pos <- SNP_inters_bin %>%
  group_by(AllChr,AllSt,cell) %>%
  filter(zscore >0) %>%
  summarise(mPosZscore = mean(zscore), .groups = "keep")
head(SNP_inters_bin_pos)
SNP_inters_bin_neg <- SNP_inters_bin %>%
  group_by(AllChr,AllSt,cell) %>%
  filter(zscore <0) %>%
  summarise(mNegZscore = mean(zscore), .groups = "keep")
head(SNP_inters_bin_neg)
SNP_inters_bin_all <- merge(SNP_inters_bin_pos, SNP_inters_bin_neg, by = c("AllChr","AllSt","cell"))
head(SNP_inters_bin_all)

SNP_inters_count <- rbind(anchD,tarD) %>% select(cell,AllChr,AllSt) %>%
  group_by(AllChr,AllSt,cell) %>%
  summarise(numInters = n(), .groups = "keep")
head(SNP_inters_count)

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
inters_univ <- ALL_inters %>% filter(!ID %in% unique(SNP_inters$ID)) 
head(inters_univ)


print("TEST OmniC dat")
SNP_inters_count %>% filter(cell == "OmniC_pooled_M_0d") %>% na.omit()
SNP_inters_bin_all %>% filter(cell == "OmniC_pooled_M_0d") %>% na.omit()

#test between random inters and SNP inters
for (c in cells_sub$V1){
  #c = "Cardiac_mesoderm_cell_day05_Zhang"
  print("############")
  print(c)
  print("############")
  #num inters in SNPs df
  numInters <- SNP_inters %>% filter(cell == c) %>% nrow()
  #SNPs info
  tSNP_inters_count <- SNP_inters_count %>% filter(cell == c)
  tSNP_inters_bin_all <- SNP_inters_bin_all %>% filter(cell == c)
  
  set.seed(369)
  #get same num of inters from univ
  rinters <- sample_n(inters_univ %>% filter(cell == c) %>% na.omit(), numInters)
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
  #separating pos and neg zscores
  rinters_bin_pos <- rinters_bin %>%
    group_by(AllChr,AllSt,cell) %>%
    filter(zscore >0) %>%
    summarise(mPosZscore = mean(zscore), .groups = "keep")
  rinters_bin_neg <- rinters_bin %>%
    group_by(AllChr,AllSt,cell) %>%
    filter(zscore <0) %>%
    summarise(mNegZscore = mean(zscore), .groups = "keep")
  rinters_bin_all <- merge(rinters_bin_pos, rinters_bin_neg, by = c("AllChr","AllSt","cell"))
  rinters_count <- rbind(anchD,tarD) %>% select(cell,AllChr,AllSt) %>%
    group_by(AllChr,AllSt,cell) %>%
    summarise(numInters = n(), .groups = "keep")
  
  print("#Mann-Whitney Test between num inters of SNP inters and random inters")
  print(wilcox.test(tSNP_inters_count$numInters, rinters_count$numInters))
  print("summary: SNP inters")
  print(summary(tSNP_inters_count))
  print("summary: random inters")
  print(summary(rinters_count))
  
  print("summary: SNP inters zscore")
  print(summary(tSNP_inters_bin_all))
  print("summary: random inters zscore")
  print(summary(rinters_bin_all))
  print("#Mann-Whitney Test between POSITIVE z-score of SNP inters and random inters")
  print(wilcox.test(tSNP_inters_bin_all$mPosZscore, rinters_bin_all$mPosZscore))
  print("#Mann-Whitney Test between NEGATIVE z-score of SNP inters and random inters")
  print(wilcox.test(tSNP_inters_bin_all$mNegZscore, rinters_bin_all$mNegZscore))
}
print("DONE")
