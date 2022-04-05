########################################
# all genes present in common trans-chromosomal interactions (blood pressure)
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: SNP list of examples info (filtered from large SNP list)
#            interactions filtered for all SNPs list 
#            bin size (bp) 
#            outfile prefix (character) 
########################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
SNPex_file <- args[1]
inter_file <- args[2]
bin_size <- as.numeric(as.character(args[3]))
outfile <- args[4]

##########
library(tidyr)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(ggplot2)
library(ggforce)#for ridgeline
library(ggridges)#for ridgeline
library(ggbiplot)#for PCA
library(devtools)#for PCA
#library(multcomp) # for anova and tukey test
.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2")
library(circlize,lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2/") # for circos
library(regioneR,lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2/")#for permutation
library(factoextra, lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3/")#for PCA
library(harrypotter, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3") #for colours
#install_github("vqv/ggbiplot")
##remotes::install_github("R-CoderDotCom/ridgeline@main")
#library(ridgeline)
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
                  #legend.title = element_blank(),
                  legend.title = element_text(size = 16),
                  legend.text = element_text(size = 18),
                  #change the colour of facet label background
                  strip.background = element_rect(fill = "#E6E1EA"),
                  #remove space between facet
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
####interaction data
##options(scipen = 999)
#SNPex_file <- "VSMC_SNP_examples_info.txt"
#inter_file <- "VSMC_trans1Mb_SNPs_MERGED_overlapping_intearctions_for_circos.txt"
##germlayer_file <- "germlayer_info.txt"
#bin_size <- 1000000
#library(factoextra)#for PCA
#library(harrypotter) #for colours


#read in SNP info
SNP <- read.table(SNPex_file, header = TRUE,sep = "\t")
head(SNP)
#accounting for merged SNP files
SNP_df <- SNP
if ("logFC_comp" %in% colnames(SNP)){
  SNP_df <- SNP %>% select(chr_b38, start_b38,end_b38, logFC_comp)
} else {
  SNP_df$logFC_comp <- mean(abs(SNP_df$logFC_comp.x),abs(SNP_df$logFC_comp.y))
}#if else
print("#new col for bins that SNPs are in")
SNP_df$bin <- plyr::round_any(SNP_df$start_b38, as.numeric(as.character(bin_size)), f = floor)
colnames(SNP_df) <- c("AllChr","SNPstart","SNPend","logFC","AllSt")
head(SNP_df)
# absolute value of mean logFC per bin
lfc <- SNP_df
lfc$logFC <- abs(lfc$logFC)
lfc <- lfc %>%
  group_by(AllChr,AllSt) %>%
  dplyr::summarise(mlogFC = mean(logFC))
lfc$AllEnd <- lfc$AllSt + bin_size
head(lfc)
lfc$ID <- paste0(lfc$AllChr,".",format(lfc$AllSt, scientific = FALSE),".",format(lfc$AllSt + bin_size,scientific = FALSE))
lfc$ID <- gsub(" ","", lfc$ID)
head(lfc)
print("done formatting SNP example data")
dat <- read.table(inter_file, header = FALSE)
colnames(dat) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
dat$ID <- paste0(dat$chrA,".",format(dat$st1, scientific = FALSE),".",format(dat$end1, scientific = FALSE),".",dat$chrB,".",format(dat$st2, scientific = FALSE),".",format(dat$end2, scientific = FALSE))
dat$ID <- gsub(" ","", dat$ID)
head(dat)
print("done formatting filtered inters data")
selection <- as.logical( # 7
  rowSums( # 6
    matrix( # 5
      unlist( # 4
        lapply(lfc$ID, function(x) { #3
          as.numeric( # 2
            str_detect(dat$ID,x) # 1
          )
        })
      ), 
      ncol = length(lfc$ID), byrow = F)
  )
)

odat <- unique(dat[selection,] %>% select(-ID)) %>% na.omit()
head(odat)
print("#write to df")
write.table(odat, file = as.character(paste0(outfile,"_overlapping_intearctions_for_circos.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

print("# getting largest region of SNPs for zoomed in archplot")
beddf <- as.data.frame(t(c(unique(lfc$AllChr),min(lfc$AllSt),max(lfc$AllEnd))))
print("#write to df")
write.table(beddf, file = as.character(paste0(outfile,"_largest_SNP_region_for_zoomed_plot.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
print("DONE")
#
