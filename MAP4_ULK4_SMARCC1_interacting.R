########################################
#plotting cis interactions in arch plot format for MAP4, SMARCC1, and ULK4 
###### DID NOT USE FOR ANALYSIS ####### not working perfectly, colours are not consistant when one of the three genes are missing
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: bed-like file for interaction data (tsv, first three columns are for the anchor region (chr, start, end) second three columns are for the target region (chr,start,end), last two columns are cell and zscore information)
#            outfile prefix (character)
########################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
VSMC_dat_file <- args[1]
outprefix <- args[2]

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
#VSMC_dat_file <- "VSMC_cis50Kb_SNPs_MERGED_overlapping_intearctions_for_circos.txt"
#outprefix <- "test_MAP_ULK_SMARC"

#VSMC SNPs interactions
dat <- read.table(VSMC_dat_file)
colnames(dat) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
head(dat)

MAP_ULK <- dat %>% filter(chrA == "chr3") %>%
  filter(
  #MAP4 filter in chrA bins
    ((st1 >= 47850000 & st1 < 48100000) | (end1 >= 47850000 & end1 <= 48100000)) & 
    #ULK4 in chrB bins
    ((st2 >= 41700000 & st2 < 42000000) | (end2 >= 41700000 & end2 <= 42000000)) | 
  #MAP4 filter in chrB bins
    ((st2 >= 47850000 & st2 < 48100000) | (end2 >= 47850000 & end2 <= 48100000)) & 
    #ULK4 in chrA bins
    ((st1 >= 41700000 & st1 < 42000000) | (end1 >= 41700000 & end1 <= 42000000))
           )
print("# MAP4 ULK4 num inters pos")
MAP_ULK %>% filter(zscore >0) %>% dplyr::select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")
print("# MAP4 ULK4 num inters neg")
MAP_ULK %>% filter(zscore <0) %>% dplyr::select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")

MAP_SMARCC <- dat %>% filter(chrA == "chr3") %>%
  filter(
    #MAP4 filter in chrA bins
    ((st1 >= 47850000 & st1 < 48100000) | (end1 >= 47850000 & end1 <= 48100000)) & 
      #SMARCC1 in chrB bins
      ((st2 >= 47700000 & st2 < 47800000) | (end2 >= 47700000 & end2 <= 47800000)) | 
      #MAP4 filter in chrB bins
      ((st2 >= 47850000 & st2 < 48100000) | (end2 >= 47850000 & end2 <= 48100000)) & 
      #SMARCC1 in chrA bins
      ((st1 >= 47700000 & st1 < 47800000) | (end1 >= 47700000 & end1 <= 47800000))
  )
print("# MAP4 SMARCC1 num inters pos")
MAP_SMARCC %>% filter(zscore >0) %>% dplyr::select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")
print("# MAP4 SMARCC1 num inters neg")
MAP_SMARCC %>% filter(zscore <0) %>% dplyr::select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")

ULK_SMARCC <- dat %>% filter(chrA == "chr3") %>%
  filter(
    #ULK filter in chrA bins
      ((st1 >= 41700000 & st1 < 42000000) | (end1 >= 41700000 & end1 <= 42000000)) & 
      #SMARCC1 in chrB bins
      ((st2 >= 47700000 & st2 < 47800000) | (end2 >= 47700000 & end2 <= 47800000)) | 
      #ULK filter in chrB bins
        ((st2 >= 41700000 & st2 < 42000000) | (end2 >= 41700000 & end2 <= 42000000)) & 
      #SMARCC1 in chrA bins
      ((st1 >= 47700000 & st1 < 47800000) | (end1 >= 47700000 & end1 <= 47800000))
  )
print("# ULK4 SMARCC1 num inters pos")
ULK_SMARCC %>% filter(zscore >0) %>% dplyr::select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")
print("# ULK4 SMARCC1 num inters neg")
ULK_SMARCC %>% filter(zscore <0) %>% dplyr::select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")

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
chromstart = 41650000
chromend = 48150000
# plot each cell separately
for (i in unique(merged_df$cell)){
  #i = "Cardiac_mesoderm_cell_day05_Zhang"
  print(i)
  #positive z-scores
  tdf <- merged_df %>% filter(cell == i) %>% filter(zscore >0)
  print(unique(tdf$cellnum))
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
  tdf$zscore <- abs(tdf$zscore)
  print(unique(tdf$cellnum))
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
                     heights =tdf$zscore,plottype="loops"
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
                     heights = abs(tdf$zscore),plottype="loops",
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
