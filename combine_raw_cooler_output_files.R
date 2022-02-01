########################################
# all genes present in common trans-chromosomal interactions
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: file with full path and name of raw cooler inter freq files (text file with each file/path on a different row)
#            output file name (character)
#######################################

options(echo=F)
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
  colnames(dat) <- c("ID",cn)
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
write.table(fdf, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

print("DONE")
