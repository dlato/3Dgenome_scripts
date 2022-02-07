########################################
#plotting trans interactions in circos plot
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: bed-like file for interaction data (tsv, first three columns are for the anchor region (chr, start, end) second three columns are for the target region (chr,start,end))
#            outfile prefix (character)
########################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
outprefix <- args[2]

##########
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
#library(circlize) # for circos
#options(scipen = 999)
#dat_file <- "circos_in_data.txt"

dat <- read.table(dat_file, header = FALSE)
colnames(dat) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
summary(dat)

#separate the df into two dfs
dfA <- dat %>% select(chrA,st1,end1)
dfB <- dat %>% select(chrB,st2,end2)

#chrom info
#re-order chroms based on chrom len
chrs_len_ord <- c("chr1","chr2",
                  "chr3","chr4",
                  "chr5","chr6",
                  "chr7","chrX",
                  "chr8","chr9",
                  "chr11","chr10",
                  "chr12","chr13",
                  "chr14","chr15",
                  "chr16","chr17",
                  "chr18","chr20",
                  "chr19","chrY",
                  "chr22","chr21")
rev_chrs_len_ord <- rev(chrs_len_ord)
#chrom info: centromere (midpoint calculated from UCSC, aprox), chrom class
chrInf <- data.frame( chrom = chrs_len_ord,
                      centromere = c(123252373.5,93787431.5,
                                     90856062,50074452.5,
                                     48585285.5,60557102.5,
                                     60058972.5,61016889,
                                     45249872,43893383.5,
                                     53454152,39800499.5,
                                     35764400,17692000.5,
                                     17117352,19037747.5,
                                     36878628.5,25067566.5,
                                     18464134,28099979.5,
                                     26161912,10470308,
                                     15520235.5,11917946),
                      chrClass = c("Metacentric","Metacentric",
                                   "Metacentric","Submetacentric",
                                   "Submetacentric","Submetacentric",
                                   "Submetacentric","Submetacentric",
                                   "Submetacentric","Submetacentric",
                                   "Submetacentric","Submetacentric",
                                   "Submetacentric","Acrocentric",
                                   "Acrocentric","Acrocentric",
                                   "Metacentric","Submetacentric",
                                   "Submetacentric","Metacentric",
                                   "Metacentric","Acrocentric",
                                   "Acrocentric","Acrocentric"),
                      size = c(248956422,242193529,
                               198295559,190214555,
                               181538259,170805979,
                               159345973,156040895,145138636,
                               138394717,135086622,
                               133797422,133275309,
                               114364328,107043718,
                               101991189,90338345,
                               83257441,80373285,
                               64444167,58617616,57227415,
                               50818468,46709983)
                      
)
#re-order chrInf df based on order Phil wants
p_chr_ord <- c("chr1","chr2",
               "chr3","chr4",
               "chr5","chr6",
               "chr7",
               "chr8","chr9",
               "chr11","chr10",
               "chr12","chr13",
               "chr14","chr15",
               "chr16","chr17",
               "chr18","chr20",
               "chr19","chr22",
               "chr21","chrX","chrY")
chrInf$chrom <- factor(chrInf$chrom, levels=p_chr_ord)

chrInf

print("#circos plot of common interactions")
dfA$chrA <- gsub("chr","",dfA$chrA)
dfB$chrB <- gsub("chr","",dfB$chrB)
colnames(dfA) <- c("chr","start","end")
colnames(dfB) <- c("chr","start","end")
head(dfA)
head(dfB)
unique(dfA$chr)
unique(dfB$chr)
nrow(dfA)
nrow(dfB)
pdf(paste0(outprefix,"_circos_trans_inters.pdf"), width = 14, height = 8)
circos.clear()
col_text <- "grey40"
circos.par("track.height"=0.8,gap.degree=5,cell.padding=c(0,0,0,0))
circos.initialize(factors=gsub("chr","",chrInf$chrom),
                  xlim=matrix(c(rep(0,length(chrInf$chrom)),chrInf$size),ncol=2))

print("# genomes")
circos.track(ylim=c(0,1),panel.fun=function(x,y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim),mean(ylim),chr,cex=0.5,col=col_text,
              facing="bending.inside",niceFacing=TRUE)
},bg.col="grey90",bg.border=F,track.height=0.06)
print("# genomes x axis")
brk <- seq(0,250, 50)*10^6
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.axis(h="top",major.at=brk,labels=round(brk/10^6,1),labels.cex=0.4,
              col=col_text,labels.col=col_text,lwd=0.7,labels.facing="clockwise")
},bg.border=F)
print("# add interactions to plot")
#rcols <- scales::alpha(ifelse(sign(nuc1$st1-nuc1$end1)!=sign(nuc2$st2-nuc2$end2),"#f46d43","#66c2a5"),alpha=0.4)
#rcols <- scales::alpha(ifelse(sign(nuc1$st1-nuc1$end1)!=sign(nuc2$st2-nuc2$end2),"black","red"))
#circos.genomicLink(nuc1,nuc2,col=rcols,border=NA)
circos.genomicLink(dfA,dfB)
dev.off()
#################
print("DONE")
