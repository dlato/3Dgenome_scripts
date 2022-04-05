########################################
# all genes present in common trans-chromosomal interactions (blood pressure)
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: interactions filtered by SNPs (bed-like tab separated file)
#            range to filter interactions (bed file)
#            bin size (numeric, bp)
########################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
bed_file <- args[2]
bin_size <- args[3]
outfile <- args[4]

##########
library(tidyr)
library(dplyr)
#library(stringr)
#library(GenomicRanges)
#library(ggplot2)
#library(ggforce)#for ridgeline
#library(ggridges)#for ridgeline
#library(ggbiplot)#for PCA
#library(devtools)#for PCA
##library(multcomp) # for anova and tukey test
#.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2")
#library(circlize,lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2/") # for circos
#library(regioneR,lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2/")#for permutation
#library(factoextra, lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3/")#for PCA
#library(harrypotter, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3") #for colours
##install_github("vqv/ggbiplot")
###remotes::install_github("R-CoderDotCom/ridgeline@main")
##library(ridgeline)
##########

#########################################################################
##set graph theme
#theme_set(theme_bw() + theme(strip.background =element_rect(fill="#e7e5e2")) +
#            #change size of facet header text
#            theme(strip.text = element_text(size =10.49)) +
#            theme(plot.title = element_text(hjust = 0.5, size = 18),
#                  panel.background = element_rect(fill = "white", colour = NA),
#                  panel.grid.major = element_blank(),
#                  panel.grid.minor = element_blank(),
#                  panel.spacing = unit(0.25, "lines"),
#                  axis.text=element_text(size=18),
#                  axis.title = element_text(size = 18),
#                  #plot margins
#                  #plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
#                  #for second legend on y-axis
#                  axis.text.y.right = element_text(size=18),
#                  #legend.title = element_blank(),
#                  legend.title = element_text(size = 16),
#                  legend.text = element_text(size = 18),
#                  #change the colour of facet label background
#                  strip.background = element_rect(fill = "#E6E1EA"),
#                  #remove space between facet
#                  panel.spacing.x=unit(0, "lines"),
#                  #                  legend.key = element_blank(),
#                  legend.background=element_blank(),
#                  #legend background
#                  legend.key = element_rect(fill = NA),
#                  #                  legend.position="none")
#                  legend.position="top")
#)
#########################################################################


print("#read in files")
####interaction data
##options(scipen = 999)
#Atype <- "1_vs_All"
##tissue_file <- "tissue_system_info.txt"
#dat_file <- "cis_arch_plot_test_dat.txt"
#bed_file <- "chr10q24.32.bed"
##germlayer_file <- "germlayer_info.txt"
#bin_size <- 1000000
#anno_file <- "hg38_p13_v32_annotation.txt"
#outfile <- "test_cell_list"
#SNP_file <- "VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt"
#library(factoextra)#for PCA
#library(harrypotter) #for colours

dat <- read.table(dat_file, header = FALSE)
colnames(dat) <- c("chr1","st1","end1","chr2","st2","end2","cell","zscore")
dat$cell <- as.factor(dat$cell)
print("summary of ALL sig zscores per cell type")
summary(dat)

#read in bed file with filtering region
bed <- read.table(bed_file, header = FALSE)
colnames(bed) <- c("chr","start","end")
#calculate bins that range falls in
bed$bin_s <- plyr::round_any(bed$start, as.numeric(as.character(bin_size)), f = floor)
bed$bin_e <- plyr::round_any(bed$end, as.numeric(as.character(bin_size)), f = floor)
bed

#filter inters for ones within bed file region
finters <- dat %>% filter(chr1 %in% bed$chr | chr2 %in% bed$chr) %>%
           filter((st1 >= bed$bin_s & st1 <= bed$bin_e) | (st2 >= bed$bin_s & st2 <= bed$bin_e))
head(finters)
summary(finters)
#write to table
write.table(finters, file = as.character(paste0(outfile,"_inters.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

print("DONE")
#
