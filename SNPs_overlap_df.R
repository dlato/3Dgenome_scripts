########################################
# find overlapping SNPs between two files (blood pressure project)
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: SNP file 1 (tsv)
#            SNP file 2 (tsv)
#            output file name (character, no path)
########################################

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
dat_file1 <- args[1]
dat_file2 <- args[2]
outfile <- args[3]

##########
library(tidyr)
library(dplyr)
#.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2")
#library(circlize,lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2/") # for circos
#library(regioneR,lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2/")#for permutation
#library(factoextra, lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3/")#for PCA
#library(harrypotter, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3") #for colours
#install_github("vqv/ggbiplot")
##remotes::install_github("R-CoderDotCom/ridgeline@main")
#library(ridgeline)
##########

#########################################################################
print("#read in files")
####interaction data
#dat_file1 <- "CM_test_SNPs.txt"
#dat_file2 <- "VSMC_test_SNPs.txt"
#outfile <- "CM_VSMC_merged.txt"


dat1 <- read.table(dat_file1, header = TRUE, sep = "\t")
dat2 <- read.table(dat_file2, header = TRUE, sep = "\t")
head(dat1)
head(dat2)
nrow(dat1)

mdf <- merge(dat1, dat2, by = c("snp_info", "pos", "chr_b38", "start_b38", "end_b38"))
nrow(mdf)
mdf
#write merged df to file
write.table(mdf, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



print("DONE")
