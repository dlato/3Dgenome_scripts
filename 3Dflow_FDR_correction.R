########################################
# applying FDR correction on 3Dflow p-values (following Sanyal et al. 2012: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3555147/ )
######
# Developer: Daniella F. Lato
#            email:  DaniellaLato@gmail.com
#            github: https://github.com/dlato
######
# arguments: 3Dflow p-value output data (tsv) ** make sure it is ALL pvalues! (not just sig)
########################################

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
pdat_file <- args[1]

##########
library(dplyr)
library(tidyr)
#library(GenomicRanges)
#install_github("vqv/ggbiplot")
##remotes::install_github("R-CoderDotCom/ridgeline@main")
#library(ridgeline)
##########

#########################################################################
#########################################################################


print("#read in files")
###interaction data
##zdat_file <- "test_1vsAll_dat.txt"
#pdat_file <- "test_1vsAll_pvalues.txt"
##pval <- 0.01
##outfile <- "test_outfile_filtering"

pdat <- read.table(pdat_file, header = TRUE)
head(pdat)
#FDR correction per cell
for (c in seq(2,ncol(pdat))){
pdat[,c] <- p.adjust(as.vector(pdat[,c]),"BH")
}#for
head(pdat)

#write to files
write.table(pdat, file = as.character(paste0(gsub(".txt","",pdat_file),"_FDR_correction.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
##################
print("DONE")
