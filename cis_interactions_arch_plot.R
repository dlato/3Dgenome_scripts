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
####interaction data
#library(circlize) # for circos
#options(scipen = 999)
#dat_file <- "cis_arch_plot_test_dat.txt"
#chrom = "chr10"
#outprefix = "test_cis_archplot"
#window_size <- 3000000
#SNP_bed = "ATF4.bed"

#SNP info (longest region for SNP example)
SNP_df <- read.table(SNP_bed, sep = "\t")
colnames(SNP_df) <- c("chrom","start","end")

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

dat <- read.table(dat_file, header = FALSE)
head(dat)
colnames(dat) <- c("chr1","st1","end1","chr2","st2","end2","cell","zscore")
dat$cell <- as.factor(dat$cell)
summary(dat)
head(dat)

#change cells to numbers for plotting
cdf <- data.frame(unique(dat$cell),seq(1,length(unique(dat$cell)),by=1))
colnames(cdf) <- c("c","cn")
cdf
dat$cellnum <- as.factor(cdf$cn[match(dat$cell,cdf$c)])
dat$cellnum <- as.numeric(as.character(dat$cellnum))
head(dat)
summary(dat)
#dat <- dat %>% filter(cell %in% c("Liver","Lung","IMR90","Adrenal.gland","Bladder_rep1","Bladder_rep2"))


#subdat <- dat %>% filter(cell %in% c("Liver","Lung","IMR90"))
#subdat$cellnum <- subdat$cell
#colnames(cdf) <- c("c","cn")
#subdat <- subdat %>%
#  mutate(cellnum = factor(cellnum, labels = c("1","2", "3")))
#subdat$cellnum <- as.numeric(as.character(subdat$cellnum))
#subdat
#plot arch plot for cis interactions
#pbpe = plotBedpe(subdat,chrom,chromstart,chromend,
#                   heights = subdat$zscore,plottype="loops",
#                   colorby=subdat$cellnum,
#                   colorbycol=SushiColors(length(unique(subdat$cellnum)))
#)

#plot arch plot for cis interactions
chromstart = 0
chromend = chrInf$size[which(chrInf$chrom == chrom)]
#pdf(paste0(outprefix,"_cis_archplot_",chrom,".pdf"), width = 14, height = 8)
#pbpe = plotBedpe(dat,chrom,chromstart,chromend,
#                   heights = dat$zscore,plottype="loops",
#                   colorby=dat$cellnum,
#                   colorbycol=SushiColors(length(unique(dat$cellnum)))
#)
##x-axis
#labelgenome(chrom, chromstart,chromend,n=6,scale="Mb")
##y-axis
#axis(side=2,las=2,tcl=.2)
##y-axis title
#mtext("Z-score",side=2,line=3,cex=1,font=2)
##legend
#legend("topright",inset =0.01,legend=cdf$c,
#       col=SushiColors(3)(3),pch=19,bty='n',text.font=2)
#dev.off()

# plot each cell separately
for (i in cdf$c){
  print(i)
  #i = "Adrenal.gland"
  #positive z-scores
  tdf <- dat %>% filter(cell == i) %>% filter(zscore > 0)
  pdf(paste0(outprefix,"_cis_archplot_",chrom,"_",i,"_positive_zscore.pdf"), width = 14, height = 8)
  pbpe = plotBedpe(tdf,chrom,chromstart,chromend,
                 #  ylim=c(-1,2),
                 color = opaque("black"),
                   heights = tdf$zscore,plottype="loops"
  #                 colorby=tdf$cellnum,
  #                 colorbycol=SushiColors(length(unique(tdf$cellnum)))
  )
  #plot(col = alpha(0.5))
  #x-axis
  labelgenome(chrom, chromstart,chromend,n=6,scale="Mb")
  #y-axis
  axis(side=2,las=2,tcl=.2, cex=2)
  #y-axis title
  mtext("Z-score (+)",side=2,line=3,cex=1,font=2)
  #title
  mtext(i,side=3,line=1,cex=2,font=2)
  dev.off()
  # negative z-scores
  tdf <- dat %>% filter(cell == i) %>% filter(zscore < 0)
  pdf(paste0(outprefix,"_cis_archplot_",chrom,"_",i,"_negative_zscore.pdf"), width = 14, height = 8)
  pbpe = plotBedpe(tdf,chrom,chromstart,chromend,
                   flip = TRUE,
                   color = opaque("black"),
                   heights = tdf$zscore,plottype="loops"
                   #                 colorby=tdf$cellnum,
                   #                 colorbycol=SushiColors(length(unique(tdf$cellnum)))
  )
  #x-axis
  labelgenome(chrom, chromstart,chromend,n=6,scale="Mb",side=3)
  #y-axis
  axis(side=2,las=2,tcl=.2)
  #y-axis title
  mtext("Z-score (-)",side=2,line=3,cex=1,font=2)
  #title
  mtext(i,side=1,line=1,cex=2,font=2)
  dev.off()
}

print("number of interactions with SNP per cell, positive and negative z-scores")
dat %>%
  dplyr::select(cell, zscore) %>%
  group_by(cell) %>%
  dplyr::summarise(numInters = n())

pz <- dat %>%
  filter(zscore > 0) %>%
  dplyr::select(cell, zscore) %>%
  group_by(cell) %>%
  dplyr::summarise(posZscore = n())
nz <- dat %>%
  filter(zscore < 0) %>%
  dplyr::select(cell, zscore) %>%
  group_by(cell) %>%
  dplyr::summarise(negZscore = n())
numInters_df <- merge(pz,nz, by = "cell", all = TRUE)
numInters_df <- numInters_df %>%
  replace(is.na(.), 0) %>%
  mutate(totalInters = rowSums(.[2:3]))

write.table(numInters_df, file = paste0(outprefix,"_number_of_interactions_per_cell_cis_archplot.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

print("#########")
print("# zoomed in archplot around SNP example")
print("#########")
#adjust chrom start and end to be window sized
chromstart = SNP_df$start[1] - window_size
chromend = SNP_df$end[1] + window_size
# plot each cell separately
for (i in cdf$c){
  print(i)
  #i = "Adrenal.gland"
  #positive z-scores
  tdf <- dat %>% filter(cell == i) %>% filter(zscore > 0)
  pdf(paste0(outprefix,"_cis_archplot_zoomed_",window_size,"_window_",chrom,"_",i,"_positive_zscore.pdf"), width = 14, height = 8)
  pbpe = plotBedpe(tdf,chrom,chromstart,chromend,
                   #  ylim=c(-1,2),
                   heights = tdf$zscore,plottype="loops"
                   #                 colorby=tdf$cellnum,
                   #                 colorbycol=SushiColors(length(unique(tdf$cellnum)))
  )
  #x-axis
  labelgenome(chrom, chromstart,chromend,n=6,scale="Mb")
  #y-axis
  axis(side=2,las=2,tcl=.2)
  #y-axis title
  mtext("Z-score (+)",side=2,line=3,cex=1,font=2)
  #title
  mtext(i,side=3,line=1,cex=2,font=2)
  dev.off()
  # negative z-scores
  tdf <- dat %>% filter(cell == i) %>% filter(zscore < 0)
  pdf(paste0(outprefix,"_cis_archplot_zoomed_",window_size,"_window_",chrom,"_",i,"_negative_zscore.pdf"), width = 14, height = 8)
  pbpe = plotBedpe(tdf,chrom,chromstart,chromend,
                   flip = TRUE,
                   heights = tdf$zscore,plottype="loops"
                   #                 colorby=tdf$cellnum,
                   #                 colorbycol=SushiColors(length(unique(tdf$cellnum)))
  )
  #x-axis
  labelgenome(chrom, chromstart,chromend,n=6,scale="Mb",side=3)
  #y-axis
  axis(side=2,las=2,tcl=.2)
  #y-axis title
  mtext("Z-score (-)",side=2,line=3,cex=1,font=2)
  #title
  mtext(i,side=1,line=1,cex=2,font=2)
  dev.off()
}



print("# number of sig inters per bin")
#count  each inter twice
tp_A <- dat %>% dplyr::select(chr1,st1,end1,cell,zscore) 
colnames(tp_A) <- c("chr","st","end","cell","zscore")
tp_B <- dat %>% dplyr::select(chr2,st2,end2,cell,zscore) 
colnames(tp_B) <- c("chr","st","end","cell","zscore")
tp_dat <- rbind(tp_A,tp_B)
tp_dat_sum = tp_dat %>% group_by(chr,st,cell) %>% dplyr::summarize(numInters=n(),.groups = "keep")
tp_dat_sum$chr <- gsub("chr", "", tp_dat_sum$chr)
summary(tp_dat_sum)
head(tp_dat_sum)
tpp <- tp_dat_sum
tp <- (ggplot(tpp, aes(st/1000000, y=1, fill = numInters))
       + geom_tile(aes(fill = numInters), width = 1, height = 1)
       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Total number of significant interactions", na.value = "grey")
       #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
       #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
       + labs(x = paste0("Chromosome ", unique(tpp$chr), " position [Mb]"),
              y = "",
              title = "Cis-chromosomal interactions overlapping with SNPs (significant)")
       + facet_grid(cell ~ .)
       #       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
       #         c("nonsig" = "Non-significant", "sig" = "Significant"))))
       #       + theme(axis.text.x = element_text(angle = 90))
       + expand_limits(x = c(chromstart/100000,chromend/100000))
       + scale_y_continuous(expand = c(0, 0))
       + scale_x_continuous(expand = c(0, 0))
       + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
               strip.background = element_rect(fill = "white"),
               panel.spacing = unit(0, "lines"),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
       + theme(panel.background = element_rect(fill = "grey85", colour = NA))
)
filename <- paste0(outprefix,"_numSig_inters_chrom",unique(tpp$chr),"_tickplot_zoomed_windos_",window_size,".pdf")
pdf(filename, width = 14, height = 4)
print(tp)
dev.off()

print("DONE")
