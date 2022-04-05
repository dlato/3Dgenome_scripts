########################################
# combin normalized cooler output files (for blood pressure project)
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: file with full path and name of raw cooler inter freq files (text file with each file/path on a different row)
#            output file name (character)
#######################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
cells_file <- args[1]
outfile <- args[2]

##########
library(tidyr)
library(dplyr)
library(stringr)

#.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2")
#library(circlize,lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2/") # for circos
#library(regioneR,lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2/")#for permutation
#library(factoextra, lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3/")#for PCA
#library(harrypotter, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3") #for colours
##install_github("vqv/ggbiplot")
###remotes::install_github("R-CoderDotCom/ridgeline@main")
##library(ridgeline)
###########

#########################################################################


print("#read in files")
#####interaction data
###tissue_file <- "tissue_system_info.txt"
#dat_file <- "test_trans_raw.txt"
#cells_file <- "cells_path.txt"
##germlayer_file <- "germlayer_info.txt"
#bin_size <- 1000000
#anno_file <- "hg38_p13_v32_annotation.txt"
#outfile <- "test_cell_list"
#SNP_file <- "VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt"


#read in file with path to data
cellp <- read.table(cells_file,sep ="\t")
head(cellp)
fdf <- data.frame()

#loop through each cell
for (c in 1:nrow(cellp)) {
  #get cell 
#  c = 1
  cell <- cellp[c,]
  cell <- data.frame(do.call("rbind", strsplit(as.character(cell), "/", fixed = TRUE)))
  cn <- cell[,ncol(cell)-1]
  dat <- read.table(cellp[c,], header = FALSE, sep = "\t")
  #dat <- read.table(dat_file, header = FALSE, sep = "\t")
  dat$V1 <- sub("anchor_","A", as.character(dat$V1))
  dat$V1 <- sub("_target_","B", as.character(dat$V1))
  dat$V1 <- sub("_","", as.character(dat$V1))
  dat$V1 <- gsub("\\/","\\.", as.character(dat$V1))
  head(dat)
  colnames(dat) <- c("ID",cn)
  dat[,2] <- as.numeric(format(dat[,2],scientific = F))
  head(dat)
  #if the final df (fdf) is empty
  if (c == 1){
    fdf <- rbind(fdf,dat)
  } else {
    fdf <- merge(fdf, dat, by ="ID", all=TRUE)
  }#ifelse
}#for
head(fdf)
summary(fdf)

#write data to table
write.table(fdf, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

print("DONE")
