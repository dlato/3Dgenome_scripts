########################################
# circos plot of interactions
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: 3Dflow z-score output data (tsv) ** make sure it is ALL zscores! (not just sig)
#            3Dflow p-value output data (tsv) ** make sure it is ALL pvalues! (not just sig)
#            germlayer df (tsv)
#            all interactions df (tsv)
########################################

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)

#Atype <- args[4]
# test is on some common interactions
# following this example:
# http://www.royfrancis.com/beautiful-circos-plots-in-r/
library(circlize) # for circos
library(tidyr)
library(dplyr)

#read in common interactions test df
dat_file <- "test_pairwise_dat.txt"
dat <- read.table(dat_file, header = TRUE)
dat$ID <- as.character(dat$ID)
#remove NAs and therefore only look at common inters
df <- na.omit(dat)
head(df)
#get just the ID col
df <- df %>% select(ID)
#split ID col
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
df$ID <- sub("B", "\\.B", as.character(df$ID))
dat2 <- df %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)

#remove A and B from chrom names
dat2$chrA <- gsub("A", "", dat2$chrA)
dat2$chrB <- gsub("B", "", dat2$chrB)
dat2 <- dat2 %>% select(-ID)
dat2$chrA <- gsub("chr","",dat2$chrA)
dat2$chrB <- gsub("chr","",dat2$chrB)
nuc1 <- dat2 %>% select(chrA, st1, end1)
nuc2 <- dat2 %>% select(chrB, st2, end2)
nuc1$st1 <- as.numeric(as.character(nuc1$st1))
nuc1$end1 <- as.numeric(as.character(nuc1$end1))
nuc2$st2 <- as.numeric(as.character(nuc2$st2))
nuc2$end2 <- as.numeric(as.character(nuc2$end2))
summary(nuc1)
#re-order chroms based on chrom len (except x and y which are at the end)
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
#chrom info: centromere (midpoint calculated from UCSC, aprox), chrom class
chrInf <- data.frame( chrom = p_chr_ord,
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
                               159345973,145138636,
                               138394717,135086622,
                               133797422,133275309,
                               114364328,107043718,
                               101991189,90338345,
                               83257441,80373285,
                               64444167,58617616,
                               50818468,46709983,
                               156040895,57227415)
                      
)



pdf("circos_test.pdf", width = 14, height = 8)
circos.clear()
col_text <- "grey40"
circos.par("track.height"=0.8,gap.degree=5,cell.padding=c(0,0,0,0))
circos.initialize(factors=gsub("chr","",chrInf$chrom),
                  xlim=matrix(c(rep(0,length(chrInf$chrom)),chrInf$size),ncol=2))

# genomes
circos.track(ylim=c(0,1),panel.fun=function(x,y) {
  chr=CELL_META$sector.index
  xlim=CELL_META$xlim
  ylim=CELL_META$ylim
  circos.text(mean(xlim),mean(ylim),chr,cex=0.5,col=col_text,
              facing="bending.inside",niceFacing=TRUE)
},bg.col="grey90",bg.border=F,track.height=0.06)
# genomes x axis
brk <- seq(0,250, 50)*10^6
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.axis(h="top",major.at=brk,labels=round(brk/10^6,1),labels.cex=0.4,
              col=col_text,labels.col=col_text,lwd=0.7,labels.facing="clockwise")
},bg.border=F)
# add interactions to plot
#rcols <- scales::alpha(ifelse(sign(nuc1$st1-nuc1$end1)!=sign(nuc2$st2-nuc2$end2),"#f46d43","#66c2a5"),alpha=0.4)
#rcols <- scales::alpha(ifelse(sign(nuc1$st1-nuc1$end1)!=sign(nuc2$st2-nuc2$end2),"black","red"))
#circos.genomicLink(nuc1,nuc2,col=rcols,border=NA)
circos.genomicLink(nuc1,nuc2)
dev.off()

