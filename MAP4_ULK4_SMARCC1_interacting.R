########################################
#plotting cis interactions in arch plot format
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: bed-like file for interaction data (tsv, first three columns are for the anchor region (chr, start, end) second three columns are for the target region (chr,start,end), last two columns are cell and zscore information)
#            chromosome to plot (character, i.e. chr1)
#            outfile prefix (character)
#            window for zoomed in archplots (numeric, bp) (i.e. value of 3000000 would create a 3Mb window around the SNP example)
#            bed file with SNP bin information (one large region for whole SNP example) (made from previous R script that filters interactions for specific SNP example)
########################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
chrom <- args[2]
outprefix <- args[3]
window_size <- as.numeric(as.character(args[4]))
SNP_bed <- args[5]

##########
.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.4.1.2")
library(tidyr)
library(forcats)
library(dplyr)
library(Sushi)
library(ggplot2)
library(harrypotter)
#library(karyoploteR)#for karyotype plot
#library(BRGenomics)#for karyotype plot
#install_github("vqv/ggbiplot")
##remotes::install_github("R-CoderDotCom/ridgeline@main")
#library(ridgeline)
##########

#########################################################################
print("#read in files")
CM_dat_file <- "CardioMyo_cis50Kb_SNPs_MERGED_overlapping_intearctions_for_circos.txt"
VSMC_dat_file <- "VSMC_cis50Kb_SNPs_MERGED_overlapping_intearctions_for_circos.txt"
outprefix <- "test_MAP_ULK_SMARC"

#CM SNPs interactions
dat <- read.table(CM_dat_file)
colnames(dat) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
#dat$chrA <- as.factor(dat$chrA)
#dat$chrB <- as.factor(dat$chrB)
head(dat)

MAP_ULK <- dat %>% filter(chrA == "chr3") %>%
  #MAP4 filter
  filter(((st1 >= 47850000 & st1 <= 48100000) | (st2 >= 47850000 & st2 <= 48100000)) | ((end1 >= 47850000 & end1 <= 48100000) | (end2 >= 47850000 & end2 <= 48100000))) %>%
  #ULK4 filter
  filter(((st1 >= 41700000 & st1 <= 42000000) | (st2 >= 41700000 & st2 <= 42000000)) | ((end1 >= 41700000 & end1 <= 42000000) | (end2 >= 41700000 & end2 <= 42000000)))
MAP_ULK %>% dplyr::select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")

MAP_SMARCC <- dat %>% filter(chrA == "chr3") %>%
  #MAP4 filter
  filter(((st1 >= 47850000 & st1 <= 48100000) | (st2 >= 47850000 & st2 <= 48100000)) | ((end1 >= 47850000 & end1 <= 48100000) | (end2 >= 47850000 & end2 <= 48100000))) %>%
  #SMARCC1 filter
  filter(((st1 >= 47700000 & st1 <= 47800000) | (st2 >= 47700000 & st2 <= 47800000)) | ((end1 >= 47700000 & end1 <= 47800000) | (end2 >= 47700000 & end2 <= 47800000)))
MAP_SMARCC %>% dplyr::select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")

ULK_SMARCC <- dat %>% filter(chrA == "chr3") %>%
  #ULK4 filter
  filter(((st1 >= 41700000 & st1 <= 42000000) | (st2 >= 41700000 & st2 <= 42000000)) | ((end1 >= 41700000 & end1 <= 42000000) | (end2 >= 41700000 & end2 <= 42000000)))%>%
  #SMARCC1 filter
  filter(((st1 >= 47700000 & st1 <= 47800000) | (st2 >= 47700000 & st2 <= 47800000)) | ((end1 >= 47700000 & end1 <= 47800000) | (end2 >= 47700000 & end2 <= 47800000)))
ULK_SMARCC %>% dplyr::select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")

#VSMC SNPs interactions
dat <- read.table(VSMC_dat_file)
colnames(dat) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
head(dat)

MAP_ULK <- dat %>% filter(chrA == "chr3") %>%
  #MAP4 filter
  filter(((st1 >= 47850000 & st1 <= 48100000) | (st2 >= 47850000 & st2 <= 48100000)) | ((end1 >= 47850000 & end1 <= 48100000) | (end2 >= 47850000 & end2 <= 48100000))) %>%
  #ULK4 filter
  filter(((st1 >= 41700000 & st1 <= 42000000) | (st2 >= 41700000 & st2 <= 42000000)) | ((end1 >= 41700000 & end1 <= 42000000) | (end2 >= 41700000 & end2 <= 42000000)))
MAP_ULK %>% dplyr::select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")

MAP_SMARCC <- dat %>% filter(chrA == "chr3") %>%
  #MAP4 filter
  filter(((st1 >= 47850000 & st1 <= 48100000) | (st2 >= 47850000 & st2 <= 48100000)) | ((end1 >= 47850000 & end1 <= 48100000) | (end2 >= 47850000 & end2 <= 48100000))) %>%
  #SMARCC1 filter
  filter(((st1 >= 47700000 & st1 <= 47800000) | (st2 >= 47700000 & st2 <= 47800000)) | ((end1 >= 47700000 & end1 <= 47800000) | (end2 >= 47700000 & end2 <= 47800000)))
MAP_SMARCC %>% dplyr::select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")

ULK_SMARCC <- dat %>% filter(chrA == "chr3") %>%
  #ULK4 filter
  filter(((st1 >= 41700000 & st1 <= 42000000) | (st2 >= 41700000 & st2 <= 42000000)) | ((end1 >= 41700000 & end1 <= 42000000) | (end2 >= 41700000 & end2 <= 42000000)))%>%
  #SMARCC1 filter
  filter(((st1 >= 47700000 & st1 <= 47800000) | (st2 >= 47700000 & st2 <= 47800000)) | ((end1 >= 47700000 & end1 <= 47800000) | (end2 >= 47700000 & end2 <= 47800000)))
ULK_SMARCC %>% dplyr::select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")

#MAP_ULK <- dat %>% filter(chrA == "chr3") %>%
#  #MAP4 filter
#  filter(((st1 >= 47850000 & st1 <= 48100000) | (st2 >= 47850000 & st2 <= 48100000)) | ((end1 >= 47850000 & end1 <= 48100000) | (end2 >= 47850000 & end2 <= 48100000))) %>%
#  #ULK4 filter
#  filter(((st1 >= 41700000 & st1 <= 42000000) | (st2 >= 41700000 & st2 <= 42000000)) | ((end1 >= 41700000 & end1 <= 42000000) | (end2 >= 41700000 & end2 <= 42000000)))
#MAP_ULK %>% dplyr::select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")
#
#MAP_SMARCC <- dat %>% filter(chrA == "chr3") %>%
#  #MAP4 filter
#  filter(((st1 >= 47850000 & st1 <= 48100000) | (st2 >= 47850000 & st2 <= 48100000)) | ((end1 >= 47850000 & end1 <= 48100000) | (end2 >= 47850000 & end2 <= 48100000))) %>%
#  #SMARCC1 filter
#  filter(((st1 >= 47700000 & st1 <= 47800000) | (st2 >= 47700000 & st2 <= 47800000)) | ((end1 >= 47700000 & end1 <= 47800000) | (end2 >= 47700000 & end2 <= 47800000)))
#MAP_SMARCC %>% dplyr::select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")
#
#ULK_SMARCC <- dat %>% filter(chrA == "chr3") %>%
#  #ULK4 filter
#  filter(((st1 >= 41700000 & st1 <= 42000000) | (st2 >= 41700000 & st2 <= 42000000)) | ((end1 >= 41700000 & end1 <= 42000000) | (end2 >= 41700000 & end2 <= 42000000)))%>%
#  #SMARCC1 filter
#  filter(((st1 >= 47700000 & st1 <= 47800000) | (st2 >= 47700000 & st2 <= 47800000)) | ((end1 >= 47700000 & end1 <= 47800000) | (end2 >= 47700000 & end2 <= 47800000)))
#ULK_SMARCC %>% dplyr::select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")

#combine all dfs into one df for arch plot
MAP_ULK$cellnum <- rep(1,nrow(MAP_ULK))
MAP_SMARCC$cellnum <- rep(2,nrow(MAP_SMARCC))
ULK_SMARCC$cellnum <- rep(3,nrow(ULK_SMARCC))
merged_df <- rbind(MAP_ULK,MAP_SMARCC,ULK_SMARCC)
colnames(merged_df) <- c("chr1","st1","end1","chr2","st2","end2","cell","zscore","cellnum")
#merged_df$cell <- as.factor(merged_df$cell)
summary(merged_df)

#arch plot with Sushi package
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

chromstart = 0
chromend = chrInf$size[which(chrInf$chrom == "chr3")]
# plot each cell separately
for (i in unique(merged_df$cell)){
  #i = "Cardiac_mesoderm_cell_day05_Zhang"
  print(i)
  #positive z-scores
  tdf <- merged_df %>% filter(cell == i) %>% filter(zscore >0)
  if (length(unique(tdf$cellnum))==1){
    #figure out the colour
    cn <- unique(tdf$cellnum)
    ccol <- ifelse(cn == 1, "blue",ifelse(cn == 2, "red","orange"))
    pdf(paste0(outprefix,"_cis_archplot_chr3_",i,"_positive_zscore.pdf"), width = 14, height = 8)
    pbpe = plotBedpe(tdf,"chr3",chromstart,chromend,
                     #  ylim=c(-1,2),
                     #color = ccol,
                     color = ccol,
                     heights = tdf$zscore,plottype="loops"
    )
    #plot(col = alpha(0.5))
    #x-axis
    labelgenome("chr3", chromstart,chromend,n=6,scale="Mb")
    #y-axis
    axis(side=2,las=2,tcl=.2, cex=2)
    #y-axis title
    mtext("Z-score (+)",side=2,line=3,cex=1,font=2)
    #title
    mtext(i,side=3,line=1,cex=2,font=2)
    dev.off()
  } else {
    pdf(paste0(outprefix,"_cis_archplot_chr3_",i,"_positive_zscore.pdf"), width = 14, height = 8)
    pbpe = plotBedpe(tdf,chrom,chromstart,chromend,
                     #  ylim=c(-1,2),
                     #color = "blue",
                     heights = tdf$zscore,plottype="loops",
                     colorby=tdf$cellnum,
                     colorbycol=SushiColors(length(unique(tdf$cellnum)))
    )
    #plot(col = alpha(0.5))
    #x-axis
    labelgenome("chr3", chromstart,chromend,n=6,scale="Mb")
    #y-axis
    axis(side=2,las=2,tcl=.2, cex=2)
    #y-axis title
    mtext("Z-score (+)",side=2,line=3,cex=1,font=2)
    #title
    mtext(i,side=3,line=1,cex=2,font=2)
    #legend
    legend("topright",inset =0.01,legend=c("MAP4 & ULK4","MAP4 & SMARCC1","ULK4 & SMARCC1"),
           col=SushiColors(3)(3),pch=19,bty='n',text.font=2)
    dev.off()
  } #ifelse
  #negative zscores
  tdf <- merged_df %>% filter(cell == i) %>% filter(zscore <0)
  if (length(unique(tdf$cellnum))==1){
    #figure out the colour
    cn <- unique(tdf$cellnum)
    ccol <- ifelse(cn == 1, "blue",ifelse(cn == 2, "red","orange"))
    pdf(paste0(outprefix,"_cis_archplot_chr3_",i,"_negative_zscore.pdf"), width = 14, height = 8)
    pbpe = plotBedpe(tdf,"chr3",chromstart,chromend,
                     #  ylim=c(-1,2),
                     #color = ccol,
                     flip = TRUE,
                     color = ccol,
                     heights = tdf$zscore,plottype="loops"
    )
    #plot(col = alpha(0.5))
    #x-axis
    labelgenome("chr3", chromstart,chromend,n=6,scale="Mb")
    #y-axis
    axis(side=2,las=2,tcl=.2, cex=2)
    #y-axis title
    mtext("Z-score (-)",side=2,line=3,cex=1,font=2)
    #title
    mtext(i,side=3,line=1,cex=2,font=2)
    dev.off()
  } else {
    pdf(paste0(outprefix,"_cis_archplot_chr3_",i,"_negative_zscore.pdf"), width = 14, height = 8)
    pbpe = plotBedpe(tdf,"chr3",chromstart,chromend,
                     #  ylim=c(-1,2),
                     #color = "blue",
                     flip = TRUE,
                     heights = tdf$zscore,plottype="loops",
                     colorby=tdf$cellnum,
                     colorbycol=SushiColors(length(unique(tdf$cellnum)))
    )
    #plot(col = alpha(0.5))
    #x-axis
    labelgenome("chr3", chromstart,chromend,n=6,scale="Mb")
    #y-axis
    axis(side=2,las=2,tcl=.2, cex=2)
    #y-axis title
    mtext("Z-score (-)",side=2,line=3,cex=1,font=2)
    #title
    mtext(i,side=3,line=1,cex=2,font=2)
    #legend
    legend("topright",inset =0.01,legend=c("MAP4 & ULK4","MAP4 & SMARCC1","ULK4 & SMARCC1"),
           col=SushiColors(3)(3),pch=19,bty='n',text.font=2)
    dev.off()
  } #ifelse
}

print("# DONE")