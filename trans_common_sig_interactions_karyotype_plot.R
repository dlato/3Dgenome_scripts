########################################
# trans common interactions karyotype plot
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: bed-like file with location of common interaction bin (chrom, start, end) and the mean zscore for that bin (tsv)
########################################

options(echo=F)
options(scipen = 1000000)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]

##########
.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.4.1.2")
library(tidyr)
library(dplyr)
library(karyoploteR)#for karyotype plot
#library(circlize,lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2/") # for circos
#library(regioneR,lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2/")#for permutation
#library(factoextra, lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3/")#for PCA
#library(harrypotter, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3") #for colours
#install_github("vqv/ggbiplot")
##remotes::install_github("R-CoderDotCom/ridgeline@main")
#library(ridgeline)
##########

#########################################################################
#########################################################################


print("#read in files")
#interaction data
#dat_file <- "merged.bed"

dat <- read.table(dat_file, header = FALSE)
colnames(dat) <- c("chrom","start","end","mzscore")
dat
summary(dat)
type(dat)
Gdat <- toGRanges(dat)
Gdat
###############
# karyotype of common interactions
##############
pdf("trans_common_interactions_mzscore_karyotype.pdf", width = 14, height = 8)
kp <- plotKaryotype(genome = "hg38")
kpAddBaseNumbers(kp)
kpHeatmap(kp, data=Gdat,y=Gdat$mzscore, colors = c("#FFAA00","#5F0B32"),r0=0.05, r1=0.6)
dev.off()
#kpPlotRegions(kp, data=dat, col="#F26419")

print("DONE")
