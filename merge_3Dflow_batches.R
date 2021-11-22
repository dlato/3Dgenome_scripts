########################################
#merging zscore or pvalue df output from 3Dflow
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: full path and file name for output merged dataframe
#            all dataframes to merge with full path, each must be separated by space. Can merge unlimited dfs
########################################

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

##########
library(tidyverse)
##########

#########################################################################


print("#read in files")
#args <- list("testFileName.txt", "test_df1.txt","test_df2.txt")
dfargs <- args[-1]
#list of dfs
l_dfs <- lapply(dfargs, read.table,header =TRUE)
mdat <- reduce(l_dfs, full_join, by = "ID")
write.table(mdat, file = args[[1]], sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
