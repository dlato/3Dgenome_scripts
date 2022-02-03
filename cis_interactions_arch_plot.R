########################################
#plotting cis interactions in arch plot format
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: bed-like file for interaction data (tsv, first three columns are for the anchor region (chr, start, end) second three columns are for the target region (chr,start,end), last two columns are cell and zscore information)
#            germlayer df (tsv)
#            all interactions file (tsv)
#            tissue/system df (tsv)
#            bin size (bp)
#            full path and name of output file
########################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]

##########
library(Sushi)

library(tidyr)
library(dplyr)
.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2")
library(circlize,lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2/") # for circos
#library(karyoploteR)#for karyotype plot
#library(BRGenomics)#for karyotype plot
#install_github("vqv/ggbiplot")
##remotes::install_github("R-CoderDotCom/ridgeline@main")
#library(ridgeline)
##########

#########################################################################
print("#read in files")
###interaction data
library(circlize) # for circos
options(scipen = 999)
dat_file <- "cis_arch_plot_test_dat.txt"



dat <- read.table(dat_file, header = FALSE)
head(dat)
colnames(dat) <- c("chr1","st1","end1","chr2","st2","end2","cell","zscore")
dat$cell <- as.factor(dat$cell)
summary(dat)
subdat <- dat %>% filter(cell %in% c("Liver","Lung","IMR90"))
subdat$cellnum <- subdat$cell
#change cells to numbers for plotting
subdat <- subdat %>%
  mutate(cellnum = factor(cellnum, labels = c("1","2", "3")))
subdat$cellnum <- as.numeric(as.character(subdat$cellnum))
subdat

#plot arch plot for cis interactions
chrom = "chr10"
#chromstart = min(c(min(dat$st1),min(dat$st2))) - 10000
#chromend = max(c(max(dat$st1),max(dat$st2))) + 10000
chromstart = 0
chromend = 133000000
pbpe = plotBedpe(subdat,chrom,chromstart,chromend,
                   heights = subdat$zscore,plottype="loops",
                   colorby=subdat$cellnum,
                   colorbycol=SushiColors(length(unique(subdat$cellnum)))
)
#x-axis
labelgenome(chrom, chromstart,chromend,n=3,scale="Mb")
#y-axis
axis(side=2,las=2,tcl=.2)
#y-axis title
mtext("Z-score",side=2,line=3,cex=1,font=2)
#legend
legend("topright",inset =0.01,legend=c("K562","HeLa","GM12878"),
       col=SushiColors(3)(3),pch=19,bty='n',text.font=2)


print("DONE")
