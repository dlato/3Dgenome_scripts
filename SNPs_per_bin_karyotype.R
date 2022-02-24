########################################
# all genes present in common trans-chromosomal interactions
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: CM SNP file (tsv)
#            VSMC SNP file (tsv)
#            output file prefix
#            bin size (numeric, bp)
########################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
CM_SNP_file <- args[1]
VSMC_SNP_file <- args[2]
outfile <- args[3]
bin_size <- args[4]

##########
.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.4.1.2")
library(tidyr)
library(dplyr)
library(karyoploteR)#for karyotype plot
library(GenomicRanges)
#install_github("vqv/ggbiplot")
##remotes::install_github("R-CoderDotCom/ridgeline@main")
#library(ridgeline)
##########

#########################################################################
print("#read in files")
####interaction data
#options(scipen = 999)
#bin_size <- 1000000
#CM_SNP_file <- "CM_test_SNPs.txt"
#outfile <- "test_num_snps"
#VSMC_SNP_file <- "VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt"
#library(tidyr)#for PCA
#library(dplyr) #for colours
#library(karyoploteR)#for PCA
#library(GenomicRanges)#for plotting

#read in SNP info
CM_SNP <- read.table(CM_SNP_file, header = TRUE,sep = "\t")
print("#bins that SNPs are in")
CM_SNP$bin <- plyr::round_any(CM_SNP$start_b38, as.numeric(as.character(bin_size)), f = floor)
CM_num_df <- CM_SNP %>% dplyr::select(chr_b38, bin) %>% group_by(chr_b38,bin) %>% summarise(CM_num = n(), .groups = "keep") 
colnames(CM_num_df) <- c("chrom", "start","CM_num")
CM_num_df$end <- CM_num_df$start + bin_size
CM_num_df <- CM_num_df %>% dplyr::select(chrom, start, end, CM_num)
VSMC_SNP <- read.table(VSMC_SNP_file, header = TRUE,sep = "\t")
print("#bins that SNPs are in")
VSMC_SNP$bin <- plyr::round_any(VSMC_SNP$start_b38, as.numeric(as.character(bin_size)), f = floor)
VSMC_num_df <- VSMC_SNP %>% dplyr::select(chr_b38, bin) %>% group_by(chr_b38,bin) %>% summarise(VSMC_num = n(), .groups = "keep") 
colnames(VSMC_num_df) <- c("chrom", "start","VSMC_num")
VSMC_num_df$end <- VSMC_num_df$start + bin_size
VSMC_num_df <- VSMC_num_df %>% dplyr::select(chrom, start, end, VSMC_num)
#merge dfs
dat <- merge(VSMC_num_df, CM_num_df, by = c("chrom", "start", "end"))
write.table(dat, file = as.character(paste0(outprefix,"_num_SNPs_per_",bin_size,"bp_bin.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = true)
#convert to GRanges
gdat <- toGRanges(as.data.frame(dat))
head(gdat)

###############
# karyotype of number of SNPs per bin
##############
pdf(paste0(outprefix,"_num_SNPs_per_",bin_size,"bp_bin_karyotype.pdf"), width = 14, height = 8)
kp <- plotKaryotype(genome = "hg38", plot.type = 1)
kpAddMainTitle(kp, paste0("Number of variants per ",bin_size,"bp bin"), cex=0.8)
kpAddBaseNumbers(kp)
at <- autotrack(current.track = 1, total.tracks = 2, margin = 0.2)
kpAddLabels(kp, labels = "VSMC", r0=at$r0, r1=at$r1, cex = 0.7)
kpHeatmap(kp, data=gdat,y=gdat$VSMC_num, colors = c("#C3BACA","#800080"),r0=0, r1=0.305)
at <- autotrack(current.track = 2, total.tracks = 2, margin = 0.2)
kpAddLabels(kp, labels = "CM", r0=at$r0, r1=at$r1, cex = 0.7)
kpHeatmap(kp, data=gdat,y=gdat$CM_num, colors = c("#FFD199","#FF8C00"),r0=0.355, r1=0.66)
dev.off()

for (c in unique(dat$chrom)){
  #c = "chr1"
  pdf(paste0(outprefix,"_num_SNPs_per_",bin_size,"bp_bin_",c,"_karyotype.pdf"), width = 14, height = 4)
  kp <- plotKaryotype(genome = "hg38", plot.type = 1, chromosomes = c(c))
  kpAddMainTitle(kp, paste0("Number of variants per ",bin_size,"bp bin"), cex=0.8)
  kpAddBaseNumbers(kp)
  at <- autotrack(current.track = 1, total.tracks = 2, margin = 0.2)
  kpAddLabels(kp, labels = "VSMC", r0=at$r0, r1=at$r1, cex = 0.7)
  kpHeatmap(kp, data=gdat,y=gdat$VSMC_num, colors = c("#C3BACA","#800080"),r0=0, r1=0.45)
  at <- autotrack(current.track = 2, total.tracks = 2, margin = 0.2)
  kpAddLabels(kp, labels = "CM", r0=at$r0, r1=at$r1, cex = 0.7)
  kpHeatmap(kp, data=gdat,y=gdat$CM_num, colors = c("#FFD199","#FF8C00"),r0=0.5, r1=1)
  dev.off()
}
print("DONE")
#
