########################################
# create supplemental data frame with inversions info and SNP rsIDs per bin (blood pressure)
# will create two new columns: one for SNPs in each region of the interaction (i.e. one list for chrA region and one list for chrB region)
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: significant z-scrore interactions (from 3Dflow output)
#            cells to look for SNPs in (text file with each cell name on a different row)
#            output file prefix
#            SNP file (tsv)
#            bin size (numeric, bp)
########################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
cells_file <- args[2]
outfile <- args[3]
SNP_file <- args[4]
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
#dat_file <- "test_1vsAll_dat.txt"
#cells_file <- "cell_subset.txt"
##germlayer_file <- "germlayer_info.txt"
#bin_size <- 1000000
#outfile <- "test_cell_list"
#SNP_file <- "VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt"
#library(factoextra)#for PCA
#library(harrypotter) #for colours

#read in cells to filter
cells_sub <- read.table(cells_file)
cells_sub
#read in SNP info
SNP <- read.table(SNP_file, header = TRUE,sep = "\t")
#accounting for merged SNP files
head(SNP)
print("#bins that SNPs are in")
SNP$bin <- plyr::round_any(SNP$start_b38, bin_size, f = floor)
SNP_df <- SNP %>% select(snp_info,chr_b38,bin) %>%
          mutate(b_end = bin + bin_size) %>%
          mutate(ID = paste0(chr_b38,".",bin,".",b_end))# %>%
          #select(-chr_b38,-bin,-b_end)
#list all SNP IDs in one column
SNP_df <- SNP_df %>% 
  select(ID,snp_info) %>%
  group_by(ID) %>% 
  dplyr::mutate(SNP_rsIDs = paste0(snp_info, collapse = ",")) %>%
  select(-snp_info) %>%
  unique()
head(SNP_df)
#read in interaction df
dat <- read.table(dat_file, header = TRUE)
print("summary of ALL sig zscores per cell type")
summary(dat)
print("#split ID col")
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
dat$ID <- sub("B", "\\.B", as.character(dat$ID))
dat <- dat %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
print("#remove A and B from chrom names")
dat$chrA <- gsub("A", "", dat$chrA)
dat$chrB <- gsub("B", "", dat$chrB)
dat <- dat %>% mutate(IDa = paste0(chrA,".",st1,".",end1)) %>%
       mutate(IDb = paste0(chrB,".",st2,".",end2)) %>%
       select(-chrA, -st1, -end1,-chrB, -st2, -end2)
head(dat)
#combine SNPs to interaction df
colnames(SNP_df) <- c("IDa","SNP_rsID_A")
datC <- full_join(dat, SNP_df, by="IDa")
colnames(SNP_df) <- c("IDb","SNP_rsID_B")
datC <- full_join(datC, SNP_df, by="IDb") %>% select(-IDa, -IDb)
head(datC)
#select only the cells we are interested in
datt <- datC[,c(1,which(colnames(datC) %in% cells_sub$V1),ncol(datC)-1,ncol(datC))]
head(datt)
print("#write to df")
write.table(datt, file = as.character(paste0(outfile,"_interactions_and_SNP_rsID.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = T)

##########
print("DONE")
#
