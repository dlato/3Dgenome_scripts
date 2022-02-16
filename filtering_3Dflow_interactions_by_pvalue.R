########################################
# filtering 3Dflow zscores based on a new p-value cutoff (so 3Dflow does not need to be run again)
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: 3Dflow z-score output data (tsv) ** make sure it is ALL zscores! (not just sig)
#            3Dflow p-value output data (tsv) ** make sure it is ALL pvalues! (not just sig)
#            new p-value cutoff 
#            outfile prefix (character) 
########################################

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
zdat_file <- args[1]
pdat_file <- args[2]
pval <- as.numeric(as.character(args[3]))
outfile <- args[4]

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
##interaction data
#zdat_file <- "test_1vsAll_dat.txt"
#pdat_file <- "test_1vsAll_pvalues.txt"
#pval <- 0.01
#outfile <- "test_outfile_filtering"

zdat <- read.table(zdat_file, header = TRUE)
pdat <- read.table(pdat_file, header = TRUE)
#combine zscore and pvalue into one df
zdatL <- gather(zdat, key = "cell", value = "zscore", 2:length(colnames(zdat)))
pdatL <- gather(pdat, key = "cell", value = "pvalue", 2:length(colnames(pdat)))
dat <- merge(zdatL, pdatL, by=c("ID","cell"))
head(dat)

#filter for input pvalue cutoff
p_df <- dat %>% filter(pvalue <= pval)
head(p_df)
summary(p_df)
#change back to wide format
pv_df <- p_df %>% select(ID,cell,pvalue) %>% spread(cell,pvalue)
z_df <- p_df %>% select(ID,cell,zscore) %>% spread(cell,zscore)
head(z_df)
summary(z_df)

#write to files
write.table(pv_df, file = as.character(paste0(outfile,"_zscore_pvalue_cutoff_",pval,".txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(z_df, file = as.character(paste0(outfile,"_pvalues_pvalue_cutoff_",pval,".txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

##################
print("DONE")
