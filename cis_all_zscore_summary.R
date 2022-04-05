########################################
#summarizing all cis zscores
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: 3Dflow z-score output data (tsv) ** make sure it is ALL zscores! (not just sig)
#            3Dflow p-value output data (tsv) ** make sure it is ALL pvalues! (not just sig)
#            germlayer df (tsv, two columns: one with the cell name as it appears in the 3Dflow output files, and the second is the associated germ layer for that cell)
#            all interactions df (tsv, all possible cis interactions, created with another R script)
########################################

#read in command line args
options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
zdat_file <- args[1]
pdat_file <- args[2]
germlayer_file <- args[3]
allinters_file <- args[4]
#Atype <- args[4]

##########
library(dplyr)
library(tidyr)
#library(GenomicRanges)
library(ggplot2)
library(ggforce)#for ridgeline
library(ggridges)#for ridgeline
library(ggbiplot, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3")#for PCA
library(circlize, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2")#for circos
library(gtools)#for sorting alphanumeric
library(devtools)#for PCA
library(factoextra,lib="/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3")#for PCA
library(harrypotter, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3")
#install_github("vqv/ggbiplot")
##remotes::install_github("R-CoderDotCom/ridgeline@main")
#library(ridgeline)
##########

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
                  legend.text = element_text(size = 10),
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
######### for testing only ###########
##interaction data
#zdat_file <- "test_cis_zscore.txt"
#pdat_file <- "test_cis_pvalues.txt"
#allinters_file <- "all_cis_interactions_1Mb.txt"
#germlayer_file <- "germlayer_info.txt"
#gl_df <- read.table("germlayer_info.txt",sep = "\t", header = TRUE)
#library(harrypotter)
#library(gtools)
####################
zdat <- read.table(zdat_file, header = TRUE)
pdat <- read.table(pdat_file, header = TRUE)
print("reading all inters")
print(allinters_file)
allinters <- read.table(allinters_file, header = FALSE)
print("summary of ALL zscores per cell type")
summary(zdat)
print("summary of ALL pvalues per cell type")
summary(pdat)
#combine zscore and pvalue into one df
zdatL <- gather(zdat, key = "cell", value = "zscore", 2:length(colnames(zdat)))
pdatL <- gather(pdat, key = "cell", value = "pvalue", 2:length(colnames(pdat)))
dat <- merge(zdatL, pdatL, by=c("ID","cell"))
head(dat)
#read in data with all possible interactions
allinters <- read.table(allinters_file, header = FALSE)
colnames(allinters) <- c("chrA", "startA", "endA", "chrB", "startB", "endB")
head(allinters)
print("summary of ALL sig zscores per cell type")
summary(dat %>% filter(pvalue <= 0.05))
###
# ALL zscores
###
print("highest zscore info per cell: ALL zscores")
# sorting the cell and zscore.
sorted_dat <- dat[order(dat$cell, -dat$zscore),]
# removing duplicates
ans <- sorted_dat[!duplicated(sorted_dat$cell),]
ans
print("lowest zscore info per cell: ALL zscores")
# sorting the cell and zscore.
sorted_dat <- dat[order(dat$cell, dat$zscore),]
# removing duplicates
ans <- sorted_dat[!duplicated(sorted_dat$cell),]
ans
###
# SIG zscores
###
sigD <- dat %>% filter(pvalue <= 0.05)
print("highest zscore info per cell: SIG zscores")
# sorting the cell and zscore.
sorted_dat <- sigD[order(sigD$cell, -sigD$zscore),]
# removing duplicates
ans <- sorted_dat[!duplicated(sorted_dat$cell),]
ans
print("lowest zscore info per cell: SIG zscores")
# sorting the cell and zscore.
sorted_dat <- sigD[order(sigD$cell, sigD$zscore),]
# removing duplicates
ans <- sorted_dat[!duplicated(sorted_dat$cell),]
ans
###
# NON-SIG zscores
###
sigD <- dat %>% filter(pvalue >= 0.05)
print("highest zscore info per cell: NON-SIG zscores")
# sorting the cell and zscore.
sorted_dat <- sigD[order(sigD$cell, -sigD$zscore),]
# removing duplicates
ans <- sorted_dat[!duplicated(sorted_dat$cell),]
ans
print("lowest zscore info per cell: NON-SIG zscores")
# sorting the cell and zscore.
sorted_dat <- sigD[order(sigD$cell, sigD$zscore),]
# removing duplicates
ans <- sorted_dat[!duplicated(sorted_dat$cell),]
ans
dat$ID <- as.character(dat$ID)
#read in germlayer info file
gl_df <- read.table(germlayer_file, header = TRUE, sep = "\t")
colnames(gl_df) <- c("cell","germLayer")
summary(gl_df)
#order of germlayer
gl_ord <- c("ectoderm", "mesoderm", "endoderm", "bipotent", "ectoderm/mesoderm")
#gl_colours <- c("#0D3B66","#F4D35E","#F95738","#66999B","#EE964B")
#darker shades
gl_colours <- c("#071F36","#F2CB40","#ED2E07","#517A7B","#EA7E1F")

df <-dat
print("Total number of interactions in celltypes")
percellD <- df %>% filter(cell == "Aorta")
length(percellD$ID)
print("Total number of ALL possible interactions")
length(allinters$chrA)
print("percent of identified interactions across genome: all possible interactions")
(length(percellD$ID)/length(allinters$chrA)) *100

#split ID col
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
df$ID <- sub("B", "\\.B", as.character(df$ID))
dat2 <- df %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
#remove A and B from chrom names
dat2$chrA <- gsub("A", "", dat2$chrA)
dat2$chrB <- gsub("B", "", dat2$chrB)
#make col for sig and non sig zscores
dat2$sig <- dat2$pvalue
dat2 <- dat2 %>% mutate(
  sig = if_else(
    pvalue >= 0.05, 
    true      = "nonsig", 
    false     = "sig"
  )
)
head(dat2)
#################
print("# PCA of zscore")
#################
print("#PCA with cells as rows")
pca_dat <- zdat
pca_dat <- na.omit(pca_dat)
row.names(pca_dat) <- pca_dat$ID
pca_dat <- pca_dat %>% select(-ID)
pca_dat <- as.data.frame(t(pca_dat))
head(pca_dat)
commonInter.pca <- prcomp(pca_dat[,c(2:ncol(pca_dat))], center = TRUE,scale. = TRUE)
summary(commonInter.pca)
g <- (ggbiplot(commonInter.pca,
               obs.scale = 1,
               var.axes=FALSE,
               var.scale = 1,
               labels = row.names(pca_dat),
               #              groups = row.names(pca_dat),
               ellipse = TRUE,
               circle = TRUE,
               ellipse.prob = 0.68
)
+ labs(title = "All Cis-chromosomal Interactions z-scores")
)
pdf("zscore_PCAcellrows_interactions_all_cells.pdf", width = 14, height = 8)
g
dev.off()

print("#PCA with cells as cols")
pca_dat <- zdat
pca_dat <- na.omit(pca_dat)
row.names(pca_dat) <- pca_dat$ID
pca_dat <- pca_dat %>% select(-ID)
summary(pca_dat)
head(pca_dat)
commonInter.pca <- prcomp(pca_dat, center = TRUE,scale. = TRUE)
summary(commonInter.pca)
g <- (ggbiplot(commonInter.pca,
               obs.scale = 1,
               #              var.axes=FALSE,
               var.scale = 1,
               #              labels = row.names(pca_dat),
               groups = colnames(commonInter.pca),
               ellipse = TRUE,
               circle = TRUE,
               ellipse.prob = 0.68
)
+ labs(title = "All Cis-chromosomal Interactions z-scores")
)
pdf("zscore_PCAcellcols_interactions_all_cells.pdf", width = 14, height = 8)
g
dev.off()

#################

################
#ridgeline plot of all chroms and number of interactions per genomic region
#################
#wide to long format
r_dat <- dat2
head(r_dat)
tail(r_dat)
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
rev_chrs_len_ord <- rev(p_chr_ord)
chrs_len_ord <- p_chr_ord
chrInf
print("#circos plot of common interactions")
ucells <- unique(dat2$cell)
for(i in ucells) {
  print(i)
#testC <- "Aorta"
circos_df <- dat2 %>% filter(cell == i) %>% filter(pvalue <= 0.05)
circos_df <- circos_df %>% select(chrA, st1, end1, chrB, st2,end2)
circos_df$chrA <- gsub("chr","",circos_df$chrA)
circos_df$chrB <- gsub("chr","",circos_df$chrB)
head(circos_df)
nuc1 <- circos_df %>% select(chrA, st1, end1)
nuc2 <- circos_df %>% select(chrB, st2, end2)
nuc1$st1 <- as.numeric(as.character(nuc1$st1))
nuc1$end1 <- as.numeric(as.character(nuc1$end1))
nuc2$st2 <- as.numeric(as.character(nuc2$st2))
nuc2$end2 <- as.numeric(as.character(nuc2$end2))
summary(nuc1)
filename <- paste0("circos_significant_cis_inters_",i,".pdf")
pdf(filename, width = 14, height = 8)
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
rcols <- scales::alpha(ifelse(sign(nuc1$st1-nuc1$end1)!=sign(nuc2$st2-nuc2$end2),"black","black"),alpha=0.5)
circos.genomicLink(nuc1,nuc2, col=rcols, border=NA)
dev.off()
}#for

r_dat$chrA <- factor(r_dat$chrA, levels=rev_chrs_len_ord)
r_dat$chrB <- factor(r_dat$chrB, levels=rev_chrs_len_ord)
r_dat$st1 <- as.numeric(as.character(r_dat$st1))
r_dat$st2 <- as.numeric(as.character(r_dat$st2))
#adding missing interactions so ridgeline lengths are accurate for each chrom
#add ID column to all inters df
allinters$ID <- with(allinters, paste0("A",chrA, ".", startA,".", endA,".B", chrB,".", startB,".", endB))
#get which inters are missing from dataset
missInters <- as.data.frame(allinters$ID[! allinters$ID %in% unique(r_dat$ID)])
colnames(missInters) <- c("ID")
head(missInters)
#newdf with just missing inters
#split ID col
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
missInters <- missInters %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
#remove A and B from chrom names
missInters$chrA <- gsub("A", "", missInters$chrA)
missInters$chrB <- gsub("B", "", missInters$chrB)
missInters$cell <- rep("fake_cell", length(missInters$ID))
missInters$sig <- rep("fake_cell", length(missInters$ID))
missInters$zscore <- rep(0, length(missInters$ID))
missInters$pvalue <- rep(1, length(missInters$ID))
missInters$st1 <- as.numeric(missInters$st1)
missInters$end1 <- as.numeric(missInters$end1)
missInters$st2 <- as.numeric(missInters$st2)
missInters$end2 <- as.numeric(missInters$end2)
missInters$zscore <- as.numeric(missInters$zscore)

r_dat <- rbind(r_dat,missInters)

head(r_dat)
r_dat3 <- r_dat
#counting each interaction twice (once for each chrom in interaction)
anchD <- r_dat
anchD$AllChr <- anchD$chrA
anchD$AllSt <- anchD$st1
tarD <- r_dat
tarD$AllChr <- tarD$chrB
tarD$AllSt <- tarD$st2
r_dat2 <- rbind(anchD,tarD)
#scale genomic position by 10Mb
r_dat2$st1 <- r_dat2$st1/10000000
head(r_dat2)
chrInf2 <- chrInf
chrInf2$chrom <- factor(chrInf2$chrom, levels=rev_chrs_len_ord)
chrInf2$centromere <- chrInf2$centromere/10000000
#chrInf2$chrom <-gsub("chr","",chrInf2$chrom)
colnames(chrInf2) <- c("AllChr","centromere","chrClass")
#chrInf2 <- chrInf2 %>% filter(AllChr %in% r_dat2$AllChr)
#p <- (ggplot(r_dat2, aes(x = AllSt, y = AllChr))
#  + geom_density_ridges(scale = 2, alpha = 0.3, rel_min_height = 0.001) 
#  + geom_segment(data = chrInf2, aes(x=centromere, xend=centromere, 
#                                     y=as.numeric(AllChr), yend=as.numeric(AllChr) +0.9),
#                color = "red")
#  + labs(x="Genomic Position [10Mb]",
#         y="",
#         title = "Location of Trans-chromosomal Interactions")
##  + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
#  + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
#  + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
#)
#######################################
## NOTE: THE RIDGELINE PLOT HEIGHT CURRENTLY IS SCALED AMONG ALL Y-AXIS POINTS AND IS THE SAME FOR ALL Y-AXIS POINTS. IT DOES NOT REPRESENT ACTUAL AMOUNTS OF INTERACTIONS
## NOTE: THE PER CHROM RIDGELINE PLOT HAS FAKE DATA ADDED SO THAT THE LENGTH OF THE CHROMOSOMES IS ACCURATE, THIS COULD MAKE THE GRAPH LOOK FUNNY! 
#######################################
#pdf("chr_ridgeline_interactions_all_cells.pdf", width = 14, height = 8)
#p
#dev.off()
###test ridgeline with ACTUAL height representing ACTUAL values of the y-axis
##d <- data.frame(x = rep(1:5, 3), y = c(rep(0, 5), rep(1, 5), rep(3, 5)),
##                height = c(0, 1, 3, 4, 0, 1, 2, 3, 5, 4, 0, 5, 4, 4, 1))
##ggplot(d, aes(x, y, height = height, group = y)) + geom_ridgeline(fill="lightblue")
##

#ridgeline of all zscores (in all chroms) across cell types
#adding germlayer info
r_dat2$germL <- r_dat2$cell
r_dat2$germL <- as.factor(gl_df$germLayer[match(r_dat2$cell, gl_df$cell)])
#re-order based on gl_ord
r_dat2$germL <- factor(r_dat2$germL, levels=gl_ord)
r_dat2 <- r_dat2[order(r_dat2$germL),]
gl_cell_ord <- unique(r_dat2$cell)
gl_cell_ord
r_dat2$cell <- factor(r_dat2$cell, levels=rev(gl_cell_ord))
#r_dat2 <- r_dat2 %>% mutate(cell = factor(cell,levels=cell))
levels(r_dat2$germL)
levels(r_dat2$cell)
head(r_dat2)
r_dat2 <- r_dat2 %>% filter(cell != "fake_cell")
#remove NAs because these mean we did not have this information in the sequencing data (since we are reading in the all zscores df)
r_dat2 <- r_dat2 %>% filter(!is.na(zscore))

#adding germlayer info
r_dat3$germL <- r_dat3$cell
r_dat3$germL <- as.factor(gl_df$germLayer[match(r_dat3$cell, gl_df$cell)])
#re-order based on gl_ord
r_dat3$germL <- factor(r_dat3$germL, levels=gl_ord)
r_dat3 <- r_dat3[order(r_dat3$germL),]
gl_cell_ord <- unique(r_dat3$cell)
gl_cell_ord
r_dat3$cell <- factor(r_dat3$cell, levels=rev(gl_cell_ord))
levels(r_dat3$germL)
levels(r_dat3$cell)
head(r_dat3)
r_dat3 <- r_dat3 %>% filter(cell != "fake_cell")
#remove NAs because these mean we did not have this information in the sequencing data (since we are reading in the all zscores df)
r_dat3 <- r_dat3 %>% filter(!is.na(zscore))


p <- (ggplot(r_dat3, aes(x = zscore, y = cell, fill = germL))
      + stat_density_ridges(quantile_lines = TRUE, alpha = 0.3, scale=2, quantiles = 2, rel_min_height = 0.001)
      #+ geom_density_ridges(scale = 4, alpha = 0.3) 
      + labs(x="z-score",
             y="Cell",
             title = "z-scores of All Cis-chromosomal Interactions",
             fill = "Germ Layer")
      + scale_fill_manual(values = gl_colours)
#                                     gsub("F", "A", my_colors)))
      + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
      + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
      + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
      + facet_grid(.~ sig, labeller = labeller(sig= as_labeller(
        c("nonsig" = "Non-Significant", "sig" = "Significant"))))
)
pdf("zscore_ridgeline_common_interactions_all_cells.pdf", width = 14, height = 8)
p
dev.off()
##################
#print("#test between means of germ layers")
#print("#### ALL DATA ####")
#print("#Test each group for normality")
#print("sig = reject normality null")
##chrClass_dat$chrClass <- as.factor(chrClass_dat$chrClass)
#r_dat3 %>%
#  group_by(germL) %>%
#  summarise(W = shapiro.test(zscore)$statistic,
#            p.value = shapiro.test(zscore)$p.value)
#print("#Perform the Kruskal-Wallis test")
#print("sig = mean is diff btwn groups")
#kruskal.test(zscore ~ germL, data=r_dat3)
#print("# check which groups have sig diff")
#print("# perform pairwise wilcoxon test with FDR (Benjamini-Hochberg) correction")
#pairwise.wilcox.test(r_dat3$zscore, r_dat3$germL,
#                     p.adjust.method = "BH")
#print("#### SIG DATA ####")
#print("#Test each group for normality")
#print("sig = reject normality null")
##chrClass_dat$chrClass <- as.factor(chrClass_dat$chrClass)
#sig_rdat <- r_dat3 %>%
#  filter(pvalue <= 0.05) 
#sig_rdat %>%
#  group_by(germL) %>%
#  summarise(W = shapiro.test(zscore)$statistic,
#            p.value = shapiro.test(zscore)$p.value)
#print("#Perform the Kruskal-Wallis test")
#print("sig = mean is diff btwn groups")
#kruskal.test(zscore ~ germL, data=sig_rdat)
#print("# check which groups have sig diff")
#print("# perform pairwise wilcoxon test with FDR (Benjamini-Hochberg) correction")
#pairwise.wilcox.test(sig_rdat$zscore, sig_rdat$germL,
#                     p.adjust.method = "BH")
#print("#### NON-SIG DATA ####")
#print("#Test each group for normality")
#print("sig = reject normality null")
##chrClass_dat$chrClass <- as.factor(chrClass_dat$chrClass)
#sig_rdat <- r_dat3 %>%
#  filter(pvalue >= 0.05) 
#sig_rdat %>%
#  group_by(germL) %>%
#  summarise(W = shapiro.test(zscore)$statistic,
#            p.value = shapiro.test(zscore)$p.value)
#print("#Perform the Kruskal-Wallis test")
#print("sig = mean is diff btwn groups")
#kruskal.test(zscore ~ germL, data=sig_rdat)
#print("# check which groups have sig diff")
#print("# perform pairwise wilcoxon test with FDR (Benjamini-Hochberg) correction")
#pairwise.wilcox.test(sig_rdat$zscore, sig_rdat$germL,
#                     p.adjust.method = "BH")
##################
##common interactions z-scores broken down by chrom class
##################
#chrClass_dat <- r_dat2
##add metacentric chrom info to df
#chrClass_dat$chrClass <- chrClass_dat$AllChr
#chrClass_dat$chrClass <- chrInf$chrClass[match(chrClass_dat$AllChr, chrInf$chrom)]
#head(chrClass_dat)
#p <- (ggplot(chrClass_dat, aes(x = chrClass, y = zscore))
#      + geom_violin(fill="grey90", scale = "count")#"count" makes width of violins proportional to number of values
#      + geom_boxplot(fill="grey95",width = 0.1,outlier.size=4)
##      + stat_density_ridges(quantile_lines = TRUE, alpha = 0.3, scale=2, quantiles = 2, rel_min_height = 0.001)
##      #+ geom_density_ridges(scale = 4, alpha = 0.3) 
#      + labs(y="z-score",
#             x="Chromosome Class",
#             title = "z-scores of Common Trans-chromosomal Interactions")
##      #  + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
##      + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
##      + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
#)
#pdf("chrClass_violin_common_interactions_all_cells.pdf", width = 14, height = 8)
#p
#dev.off()
##with cell info
##p <- (ggplot(chrClass_dat, aes(y = cell, x = zscore))
##      + geom_violin(fill="grey90")
##      #      + stat_density_ridges(quantile_lines = TRUE, alpha = 0.3, scale=2, quantiles = 2, rel_min_height = 0.001)
##      #      #+ geom_density_ridges(scale = 4, alpha = 0.3) 
##      + facet_grid(.~ chrClass)
##      + labs(x="z-score",
##             y="Cell",
##             title = "z-scores of Common Trans-chromosomal Interactions")
##      #      #  + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
##      #      + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
##      #      + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
##)
#p <- (ggplot(chrClass_dat, aes(y = cell, x = zscore, fill = germL))
#      + stat_density_ridges(quantile_lines = TRUE, alpha = 0.3, scale=2, quantiles = 2, rel_min_height = 0.001)
#      #      + stat_density_ridges(quantile_lines = TRUE, alpha = 0.3, scale=2, quantiles = 2, rel_min_height = 0.001)
#      #      #+ geom_density_ridges(scale = 4, alpha = 0.3) 
#      + facet_grid(.~ chrClass)
#      + labs(x="z-score",
#             y="Cell",
#             title = "z-scores of Common Trans-chromosomal Interactions",
#             fill = "Germ Layer")
#      + scale_fill_manual(values = gl_colours)
#      #      #  + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
#      #      + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
#      #      + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
#)
#pdf("chrClass_cell_facet_ridgeline_common_interactions_all_cells.pdf", width = 14, height = 8)
#p
#dev.off()
##################

##################
##proportional interactions per chrom
##################
#prop_dat <- r_dat2
#prop_dat %>% filter(ID == "Achr10.0.1000000.Bchr18.79000000.80000000")
###re-format all interactions df
###split ID col
##dat$ID <- sub("B", "\\.B", as.character(dat$ID))
##dat <- dat %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
###remove A and B from chrom names
##dat$chrA <- gsub("A", "", dat$chrA)
##dat$chrB <- gsub("B", "", dat$chrB)
##dat_long <- gather(dat, cell, zscore, 8:ncol(dat2), factor_key=TRUE)
###counting each interaction twice (once for each chrom in interaction)
##anchD <- dat_long
##anchD$AllChr <- anchD$chrA
##anchD$AllSt <- anchD$st1
##tarD <- dat_long
##tarD$AllChr <- tarD$chrB
##tarD$AllSt <- tarD$st2
##dat_long2 <- rbind(anchD,tarD)
###scale genomic position by 10Mb
##dat_long2$AllSt <- as.numeric(as.character(dat_long2$AllSt))/10000000
###df with total number of interactions per chrom
##totInter <- dat_long2 %>%
##  select(AllChr, ID) %>%
##  group_by(AllChr) %>%
##  dplyr::summarise(n = n())
##totInter$AllChr <- factor(totInter$AllChr, levels=rev_chrs_len_ord)
#
##re-format all interactions df
#head(allinters)
##counting each interaction twice (once for each chrom in interaction)
#anchD <- allinters
#anchD$AllChr <- anchD$chrA
#anchD$AllSt <- anchD$startA
#tarD <- allinters
#tarD$AllChr <- tarD$chrB
#tarD$AllSt <- tarD$startB
#dat_long2 <- rbind(anchD,tarD)
##head(dat_long2)
##tail(dat_long2)
##allinters %>% filter(ID == "Achr1.0.1000000.Bchr2.0.1000000")
##dat_long2 %>% filter(ID == "Achr1.0.1000000.Bchr2.0.1000000")
##allinters %>% filter(ID == "Achr22.50000000.50818468.Bchr21.44000000.45000000")
##dat_long2 %>% filter(ID == "Achr22.50000000.50818468.Bchr21.44000000.45000000")
##scale genomic position by 10Mb
#dat_long2$AllSt <- as.numeric(as.character(dat_long2$AllSt))/10000000
##df with total number of interactions per chrom
#totInter <- dat_long2 %>%
#  select(AllChr, ID) %>%
#  group_by(AllChr) %>%
#  dplyr::summarise(n = n())
#totInter$AllChr <- factor(totInter$AllChr, levels=rev_chrs_len_ord)
#totInter
##df with total number of common interactions per chrom (ALL POSSIBLE INTERACTIONS)
##remove cell info (which makes interactions duplicated)
#prop_dat2 <- prop_dat %>% select(AllChr, ID)
#prop_dat2 %>% filter(ID == "Achr10.0.1000000.Bchr13.113000000.114000000")
#prop_dat2 <- unique(prop_dat2)
#prop_dat2 %>% filter(ID == "Achr10.0.1000000.Bchr13.113000000.114000000")
#commonInter <- prop_dat2 %>%
#  select(AllChr, ID) %>%
#  group_by(AllChr) %>%
#  dplyr::summarise(n = n())
#commonInter$AllChr <- factor(commonInter$AllChr, levels=rev_chrs_len_ord)
##df combining above 2 dfs
#prop_datAll <- totInter
#prop_datAll$commonInter <- prop_datAll$n
#colnames(prop_datAll) <- c("chrom", "totInter","commonInter")
#prop_datAll$commonInter <- as.numeric(commonInter$n[match(prop_datAll$chrom, commonInter$AllChr)])
#prop_datAll$commonInter[is.na(prop_datAll$commonInter)] <- 0
#prop_datAll$percent <- (prop_datAll$commonInter/prop_datAll$totInter)*100 
#prop_datAll$chrom <- as.factor(prop_datAll$chrom)
#prop_datAll$chrom <- as.factor(gsub("chr","", prop_datAll$chrom))
#chrs_len_ord_num <- gsub("chr","",chrs_len_ord)
#prop_datAll$chrom <- factor(prop_datAll$chrom, levels=chrs_len_ord_num)
#levels(prop_datAll$chrom)
###plot
##p <- (ggplot(prop_datAll, aes(y = percent,x=chrom))
##      + geom_bar(fill="grey90", color="black", stat = "identity")
##      + labs(y="Percent of Total Interactions [%]",
##             x="Chromosome",
##             title = "Percentage of Common Trans-chromosomal Interactions per Chromosome")
##)
##pdf("proportion_per_chrom_common_interactions_all_cells.pdf", width = 14, height = 8)
##p
##dev.off()
#prop_datAll2 <- prop_datAll %>% select(chrom, totInter, commonInter) %>% gather(key = "name", value = "num", 2:3)
#head(prop_datAll2)
#prop_datAll2$name <- factor(prop_datAll2$name, levels = c("totInter", "commonInter"))
##scaling number of interactions by 100,000
#prop_datAll2$num <- prop_datAll2$num / 100000
##plot
#p <- (ggplot(prop_datAll2, aes(y =num,x=chrom, fill = name))
#      + geom_bar(color="black",position = "dodge", stat = "identity")
#      + scale_fill_manual(values =c("#FAC9A1", "#013040"), labels= c("Total Interactions", "Common Interactions"))
#      + labs(y="Number of Interactions [100,000]",
#             x="Chromosome",
#             title = "Trans-chromosomal Interactions per Chromosome",
#             fill = "")
#)
#pdf("proportion_per_chrom_common_interactions_all_cells.pdf", width = 14, height = 8)
#p
#dev.off()
#
##proportional interactions per chrom pair
#pair_dat <- dat_long2
#head(pair_dat)
##new col for chrom pair
#pair_dat$pair <- paste0(pair_dat$chrA,pair_dat$chrB)
#pair_dat <- pair_dat %>% select(ID, pair) %>% unique()
##df with total number of interactions per chrom pair ALL POSSIBLE INTERACTIONS
#totInter <- pair_dat %>%
#  select(pair, ID) %>%
#  group_by(pair) %>%
#  dplyr::summarise(n = n())
##df with total number of common interactions per chrom pair
#commonInter <- prop_dat
#head(prop_dat)
##commonInter <- gather(dat2, cell, zscore, 8:ncol(dat2), factor_key=TRUE)
##new col for chrom pair
#commonInter$pair <- paste0(commonInter$chrA,commonInter$chrB)
#commonInter <- commonInter %>% select(ID, pair) %>% unique()
#head(commonInter)
#commonInter %>% filter(ID == "Achr10.0.1000000.Bchr13.113000000.114000000")
#commonInter <- commonInter %>%
#  select(pair, ID) %>%
#  group_by(pair) %>%
#  dplyr::summarise(n = n())
#commonInter
##df combining above 2 dfs
#prop_pairAll <- totInter
#prop_pairAll$commonInter <- prop_pairAll$n
#colnames(prop_pairAll) <- c("chrom", "totInter","commonInter")
#prop_pairAll$commonInter <- as.numeric(commonInter$n[match(prop_pairAll$chrom, commonInter$pair)])
#prop_pairAll$commonInter[is.na(prop_pairAll$commonInter)] <- 0
#prop_pairAll$percent <- (prop_pairAll$commonInter/prop_pairAll$totInter)*100 
#prop_pairAll$chrom <- as.factor(prop_pairAll$chrom)
#print("# chrom pairs proportional number of common interactions")
#prop_pairAll
#print("# chromosome pair(s) with highest (proportional) number of common interactions")
#prop_pairAll[which(prop_pairAll$percent == max(prop_pairAll$percent)),]
##################

#################
#heatmap of common interactions z-scores by cell type and chromosome
#################
r_dat2$AllChr <- gsub("chr","",r_dat2$AllChr)
r_dat2$AllChr <- as.factor(r_dat2$AllChr)
#calculate mean zscore per chrom per cell type so heat map is accutate
hm_dat = r_dat2 %>% group_by(AllChr,cell, sig) %>% dplyr::summarize(mzscore=mean(zscore, na.rm = TRUE))
head(hm_dat)
summary(hm_dat)
p_chr_ord2 <- gsub("chr","", p_chr_ord) 
hm_dat2 <- hm_dat %>% mutate(AllChr=factor(AllChr, levels=p_chr_ord2))
#hm <- (ggplot(r_dat2, aes(AllChr, cell))
#       + geom_tile(aes(fill = zscore), colour = "white")
hm <- (ggplot(hm_dat2, aes(AllChr, cell, fill = zscore))
       + geom_tile(aes(fill = mzscore), width = 1, height = 1)
+ scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score", na.value = "grey")
#       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
       #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
       + labs(x = "Chromosome",
              y = "Cell",
              title = "All Cis-chromosomal Interactions z-scores")
#       + facet_wrap(.~sig)
       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
         c("nonsig" = "Non-significant", "sig" = "Significant"))))
#       + theme(axis.text.x = element_text(angle = 90))
       + theme(axis.text=element_text(size=12),panel.background = element_rect(fill = "grey85", colour = NA))
)
pdf("zscore_mean_heatmap_interactions_chroms_all_cells.pdf", width = 14, height = 8)
hm
dev.off()
###########
# mean of ALL zscores (not separtated into sig and non-sig)
###########
hm_dat = r_dat2 %>% group_by(AllChr,cell) %>% dplyr::summarize(mzscore=mean(zscore, na.rm = TRUE))
hm_dat2 <- hm_dat %>% mutate(AllChr=factor(AllChr, levels=p_chr_ord2))
#hm <- (ggplot(r_dat2, aes(AllChr, cell))
#       + geom_tile(aes(fill = zscore), colour = "white")
hm <- (ggplot(hm_dat2, aes(AllChr, cell, fill = zscore))
       + geom_tile(aes(fill = mzscore), width = 1, height = 1)
+ scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score", na.value = "grey")
#       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
       #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
       + labs(x = "Chromosome",
              y = "Cell",
              title = "Cis-chromosomal Interactions (all) z-scores")
#       + facet_wrap(.~sig)
#       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
#         c("nonsig" = "Non-significant", "sig" = "Significant"))))
#       + theme(axis.text.x = element_text(angle = 90))
       + theme(axis.text=element_text(size=12),panel.background = element_rect(fill = "grey85", colour = NA))
)
pdf("zscore_mean_heatmap_interactions_chroms_all_cells_mean_of_ALL_zscores.pdf", width = 14, height = 8)
hm
dev.off()
###########
# mean of ALL zscores (not separtated into sig and non-sig)
###########
hm_dat = r_dat2 %>% filter(sig == "sig") %>% group_by(AllChr,cell) %>% dplyr::summarize(mzscore=mean(zscore, na.rm = TRUE))
hm_dat2 <- hm_dat2 %>% mutate(AllChr=factor(AllChr, levels=p_chr_ord2))
#hm <- (ggplot(r_dat2, aes(AllChr, cell))
#       + geom_tile(aes(fill = zscore), colour = "white")
hm <- (ggplot(hm_dat, aes(AllChr, cell, fill = zscore))
       + geom_tile(aes(fill = mzscore), width = 1, height = 1)
+ scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score", na.value = "grey")
#       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
       #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
       + labs(x = "Chromosome",
              y = "Cell",
              title = "Cis-chromosomal Interactions (significant) z-scores")
#       + facet_wrap(.~sig)
#       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
#         c("nonsig" = "Non-significant", "sig" = "Significant"))))
#       + theme(axis.text.x = element_text(angle = 90))
       + theme(axis.text=element_text(size=12),panel.background = element_rect(fill = "grey85", colour = NA))
)
pdf("zscore_mean_heatmap_interactions_chroms_all_cells_significant_interactions.pdf", width = 14, height = 8)
hm
dev.off()
#################
# total interactions per chrom pair vs all possible interactions all per cell 
#################
head(r_dat)
#make chrom pair col for all inters (sig and non-sig)
bar_df <- r_dat
#bar_df$chrPair <- paste0(bar_df$chrA, bar_df$chrB)
#head(bar_df)
chrpaircount <- bar_df %>% group_by(cell, chrA) %>% dplyr::summarise(n = n())
head(chrpaircount)
colnames(chrpaircount) <- c("cell","chrPair","obsCount")
#get cells to replicate for all inters chrom pair count
cells_obs <- unique(chrpaircount$cell)

allinters_chrpair <- allinters
#allinters_chrpair$chrPair <- paste0(allinters_chrpair$chrA, allinters_chrpair$chrB)
#head(allinters_chrpair)
allchrpaircount <- allinters_chrpair %>% group_by(chrA) %>% dplyr::summarise(n = n())
colnames(allchrpaircount) <- c("chrPair","allCount")
allchrpaircount <- allchrpaircount %>% slice(rep(1:n(), each = length(cells_obs)))
allcell <- rep(cells_obs, times = length(unique(allchrpaircount$chrPair)))
tail(allchrpaircount)
allchrpaircount <- cbind(allcell, allchrpaircount)
colnames(allchrpaircount) <- c("cell","chrPair","allCount")
head(allchrpaircount)
#combine dfs
chrpair_df <- merge(allchrpaircount, chrpaircount, by = c("cell","chrPair"), all = TRUE)
chrpair_df <- gather(chrpair_df, key="interaction", value="count", 3:4)
head(chrpair_df)

chrpair_df$chrPair = as.factor(chrpair_df$chrPair)
chrpair_df$chrPair = factor(chrpair_df$chrPair, levels = gtools::mixedsort(unique(chrpair_df$chrPair)))
chrpair_df <- chrpair_df %>% filter(cell != "fake_cell")
chrpair_df$chrPair <- gsub("chr","", chrpair_df$chrPair)
bp <- (ggplot(chrpair_df, aes(x=chrPair, y=count))
  +geom_bar(aes(fill = interaction), position = "dodge", stat="identity")
#  + scale_fill_hp(discrete = FALSE, option = "Always", name = "Mean z-score per chromosomal pair", na.value = "grey")
#  #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
#  #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
  + labs(x = "Chromosome",
         y = "Number of cis-chromosomal interactions",
         title = "Total Cis-chromosomal Interactions per Chromosome",
         fill = "")
  + facet_grid(cell ~ .)
#  #       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
#  #         c("nonsig" = "Non-significant", "sig" = "Significant"))))
#  #       + theme(axis.text.x = element_text(angle = 90))
  + scale_fill_manual(values =c("#FAC9A1", "#013040"), labels= c("Total Possible Interactions", "Common Interactions"))
  + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
          strip.background = element_rect(fill = "white"),
          panel.spacing = unit(0, "lines"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
#  + theme(axis.text=element_text(size=5))
   + theme(strip.text.y.right = element_text(angle = 0))
)
pdf("number_of_inters_chrPair_obs_exp_allInters_bargraph.pdf", width = 14, height = 8)
bp
dev.off()

#################
# total interactions per cell vs all possible interactions per cell 
#################
#make chrom pair col for all inters (sig and non-sig)
bar_df <- r_dat
head(bar_df)
bar_df <- na.omit(bar_df)
chrpaircount <- bar_df %>% group_by(cell) %>% dplyr::summarise(n = n())
head(chrpaircount)
colnames(chrpaircount) <- c("cell","obsCount")
#get cells to replicate for all inters chrom pair count
cells_obs <- unique(chrpaircount$cell)

allinters_tot <- length(allinters$chrA)
allrow <- c("All possible interactions", allinters_tot)
chrpaircount <- rbind(allrow,chrpaircount) %>% filter(cell != "fake_cell")
chrpaircount$obsCount <- as.numeric(as.character(chrpaircount$obsCount)) / 100000

bp <- (ggplot(chrpaircount, aes(y=obsCount, x=cell))
       +geom_bar(position = "dodge", stat="identity")
       + coord_flip()
       #  + scale_fill_hp(discrete = FALSE, option = "Always", name = "Mean z-score per chromosomal pair", na.value = "grey")
       #  #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
       #  #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
         + labs(y = "Number of Cis-chromosomal Interactions [100,000]",
                x = "",
                title = "Cis-chromosomal Interactions")
 #      + facet_grid(cell ~ .)
       #  #       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
       #  #         c("nonsig" = "Non-significant", "sig" = "Significant"))))
       #  #       + theme(axis.text.x = element_text(angle = 90))
       #  + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
       #          strip.background = element_rect(fill = "white"),
       #          panel.spacing = unit(0, "lines"),
       #          axis.text.y = element_blank(),
       #          axis.ticks.y = element_blank())
       #  + theme(axis.text=element_text(size=5))
       + scale_y_continuous(expand = c(0, 0), limits = c(0,max(chrpaircount$obsCount +0.2)))
       + theme(strip.text.y.right = element_text(angle = 0))
)
pdf("number_of_inters_obs_exp_allInters_bargraph.pdf", width = 14, height = 8)
bp
dev.off()
#bar graph with sig and non sig total interactions per cell
bar_df <- r_dat
head(bar_df)
bar_df <- na.omit(bar_df)
chrpaircount <- bar_df %>% group_by(cell, sig) %>% dplyr::summarise(n = n())
head(chrpaircount)
colnames(chrpaircount) <- c("cell","sig","obsCount")
#get cells to replicate for all inters chrom pair count
cells_obs <- unique(chrpaircount$cell)
chrpaircount <- chrpaircount %>% filter(cell != "fake_cell")

allinters_tot <- length(allinters$chrA)
allrow <- as.data.frame(t(c("All possible interactions", "all", allinters_tot)))
colnames(allrow) <- c("cell","sig","obsCount")
chrpaircount <- rbind(allrow,chrpaircount)
chrpaircount$obsCount <- as.numeric(as.character(chrpaircount$obsCount)) / 100000

bp <- (ggplot(chrpaircount, aes(y=obsCount, x=cell, fill = sig))
       +geom_bar(position = "dodge", stat="identity")
       + coord_flip()
       + scale_fill_manual(values = c("grey","#FAC9A1", "#013040"), labels= c("All possible interactions","Total", "Significant"))
       #  + scale_fill_hp(discrete = FALSE, option = "Always", name = "Mean z-score per chromosomal pair", na.value = "grey")
       #  #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
       #  #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
       + labs(y = "Number of Cis-chromosomal Interactions [100,000]",
              x = "",
              title = "Cis-chromosomal Interactions",
              fill = "")
       #      + facet_grid(cell ~ .)
       #  #       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
       #  #         c("nonsig" = "Non-significant", "sig" = "Significant"))))
       #  #       + theme(axis.text.x = element_text(angle = 90))
       #  + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
       #          strip.background = element_rect(fill = "white"),
       #          panel.spacing = unit(0, "lines"),
       #          axis.text.y = element_blank(),
       #          axis.ticks.y = element_blank())
         + theme(axis.text=element_text(size=12))
       + scale_y_continuous(expand = c(0, 0), limits = c(0,max(chrpaircount$obsCount +0.2)))
       + theme(strip.text.y.right = element_text(angle = 0))
)
pdf("number_of_inters_obs_exp_sig_allInters_bargraph.pdf", width = 14, height = 8)
bp
dev.off()
#################
# average z-score per chrom pair heatmap 
#################
head(r_dat)
hm_dat <- r_dat
hm_dat$chrPair <- paste0(hm_dat$chrA,hm_dat$chrB)
head(hm_dat)
#calculate mean zscore per chrom pair per cell type so heat map is accutate
#hm_dat2 = hm_dat %>% filter(cell != "fake_cell") %>% group_by(chrPair,cell,sig) %>% dplyr::summarize(mzscore=mean(zscore, na.rm = TRUE))
hm_dat2 = hm_dat %>% filter(cell != "fake_cell") %>% group_by(chrPair,cell) %>% dplyr::summarize(mzscore=mean(zscore, na.rm = TRUE))
hm_dat2$chrPair <- gsub("chr", "\\.chr", hm_dat2$chrPair)
coln <- c("tmp","chrA", "chrB")
hm_dat2 <- hm_dat2 %>% separate(chrPair, sep = "\\.", into = coln, remove = FALSE) %>% select(-tmp)
hm_dat2$chrA <- gsub("chr", "", hm_dat2$chrA)
hm_dat2$chrB <- gsub("chr", "", hm_dat2$chrB)
head(hm_dat2)
#  test_D <- hm_dat2 %>% filter(cell == "Hippocampus")
#head(test_D)
hm_dat2 <- hm_dat2 %>% mutate(chrA=factor(chrA, levels=p_chr_ord2))
hm_dat2 <- hm_dat2 %>% mutate(chrB=factor(chrB, levels=p_chr_ord2))
#hm <- (ggplot(hm_dat2, aes(chrA, chrB, fill = mzscore))
hm <- (ggplot(hm_dat2, aes(chrA, y = 1, fill = mzscore))
       + geom_tile(aes(fill = mzscore), width = 1, height = 1)
       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score per chromosomal pair", na.value = "grey")
       #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
       #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
       + labs(x = "Chromosome",
              y = "",
              title = "Trans-chromosomal Interactions (all) z-scores")
       + facet_wrap(cell ~ .)
#       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
#         c("nonsig" = "Non-significant", "sig" = "Significant"))))
       #       + theme(axis.text.x = element_text(angle = 90))
      + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
              strip.background = element_rect(fill = "white"),
              panel.spacing = unit(0, "lines"),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
       + theme(axis.text=element_text(size=5),panel.background = element_rect(fill = "grey85", colour = NA))
)
pdf("zscore_chrom_pair_mean_heatmap_AllInteractions.pdf", width = 14, height = 8)
hm
dev.off()
#########
# test one cell heat map
##########
hm_dat2 = hm_dat %>% filter(cell == "Aorta") %>% group_by(chrPair) %>% dplyr::summarize(mzscore=mean(zscore, na.rm = TRUE))
hm_dat2$chrPair <- gsub("chr", "\\.chr", hm_dat2$chrPair)
coln <- c("tmp","chrA", "chrB")
hm_dat2 <- hm_dat2 %>% separate(chrPair, sep = "\\.", into = coln, remove = FALSE) %>% select(-tmp)
hm_dat2$chrA <- gsub("chr", "", hm_dat2$chrA)
hm_dat2$chrB <- gsub("chr", "", hm_dat2$chrB)
head(hm_dat2)
#  test_D <- hm_dat2 %>% filter(cell == "Hippocampus")
#head(test_D)
hm_dat2 <- hm_dat2 %>% mutate(chrA=factor(chrA, levels=p_chr_ord2))
hm_dat2 <- hm_dat2 %>% mutate(chrB=factor(chrB, levels=p_chr_ord2))
hm <- (ggplot(hm_dat2, aes(chrA, y=1, fill = mzscore))
       + geom_tile(aes(fill = mzscore), width = 1, height = 1)
       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score per chromosomal pair", na.value = "grey")
       #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
       #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
       + labs(x = "Chromosome",
              y = "",
              title = "Trans-chromosomal Interactions (all) z-scores")
#       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
#         c("nonsig" = "Non-significant", "sig" = "Significant"))))
       #       + theme(axis.text.x = element_text(angle = 90))
      + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
              strip.background = element_rect(fill = "white"),
              panel.spacing = unit(0, "lines"))
)
pdf("TEST_Aorta_zscore_chrom_pair_mean_heatmap_AllInteractions.pdf", width = 14, height = 8)
hm
dev.off()

##############
# significant interactions
##############
hm_dat2 = hm_dat %>% filter(cell != "fake_cell") %>% filter(sig == "sig") %>% group_by(chrPair,cell) %>% dplyr::summarize(mzscore=mean(zscore, na.rm = TRUE))
hm_dat2$chrPair <- gsub("chr", "\\.chr", hm_dat2$chrPair)
coln <- c("tmp","chrA", "chrB")
hm_dat2 <- hm_dat2 %>% separate(chrPair, sep = "\\.", into = coln, remove = FALSE) %>% select(-tmp)
hm_dat2$chrA <- gsub("chr", "", hm_dat2$chrA)
hm_dat2$chrB <- gsub("chr", "", hm_dat2$chrB)
hm_dat2 <- hm_dat2 %>% mutate(chrA=factor(chrA, levels=p_chr_ord2))
hm_dat2 <- hm_dat2 %>% mutate(chrB=factor(chrB, levels=p_chr_ord2))
hm <- (ggplot(hm_dat2, aes(chrA, y=1, fill = mzscore))
       + geom_tile(aes(fill = mzscore), width = 1, height = 1)
       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score per chromosomal pair", na.value = "grey")
       #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
       #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
       + labs(x = "Chromosome",
              y = "",
              title = "Trans-chromosomal Interactions (significant) z-scores")
       + facet_wrap(cell ~ .)
       #       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
       #         c("nonsig" = "Non-significant", "sig" = "Significant"))))
       #       + theme(axis.text.x = element_text(angle = 90))
       + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
               strip.background = element_rect(fill = "white"),
               panel.spacing = unit(0, "lines"),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
       + theme(axis.text=element_text(size=5),panel.background = element_rect(fill = "grey85", colour = NA))
)
pdf("zscore_chrom_pair_mean_heatmap_sig_Interactions.pdf", width = 14, height = 8)
hm
dev.off()


print("##########")
print("# tickplot per chrom, each bin represented fact with cells")
print("# mean zscore per bin")
head(r_dat)
#count  each inter twice
tp_A <- r_dat %>% select(chrA,st1,end1,cell,zscore,pvalue) 
colnames(tp_A) <- c("chr","st","end","cell","zscore","pvalue")
tp_B <- r_dat %>% select(chrB,st2,end2,cell,zscore,pvalue) 
colnames(tp_B) <- c("chr","st","end","cell","zscore","pvalue")
#select only sig inters
tp_dat <- rbind(tp_A,tp_B) %>% filter(pvalue <=0.05)
tp_dat_sum = tp_dat %>% group_by(chr,st,cell) %>% dplyr::summarize(mzscore=mean(zscore, na.rm = TRUE))
tp_dat_sum$chr <- gsub("chr", "", tp_dat_sum$chr)
tp_dat_sum <- tp_dat_sum %>% mutate(chr=factor(chr, levels=p_chr_ord2))
#scale pts by 1Mb
tp_dat_sum$st <- tp_dat_sum$st / 1000000
summary(tp_dat_sum)
head(tp_dat_sum)
for(i in unique(tp_dat_sum$chr)) {
#i="X"
xmax <- chrInf$size[match(paste0("chr",i),chrInf$chrom)] / 1000000
tpp <- tp_dat_sum %>% filter(chr == i)
tp <- (ggplot(tpp, aes(st, y=1, fill = mzscore))
       + geom_tile(aes(fill = mzscore), colour = "white")
       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score per bin", na.value = "grey")
       #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
       #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
       + labs(x = paste0("Chromosome ", i, " position [Mb]"),
              y = "",
              title = "Cis-chromosomal interactions (significant) z-scores")
       + facet_grid(cell ~ .)
       #       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
       #         c("nonsig" = "Non-significant", "sig" = "Significant"))))
       #       + theme(axis.text.x = element_text(angle = 90))
      + expand_limits(x = c(0,xmax))
       + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
               strip.background = element_rect(fill = "white"),
               panel.spacing = unit(0, "lines"),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
       + theme(axis.text=element_text(size=5),panel.background = element_rect(fill = "grey85", colour = NA))
)
filename <- paste0("zscore_chrom",i,"_mean_tickplot_sig_Interactions.pdf")
pdf(filename, width = 14, height = 8)
print(tp)
dev.off()
}#for
print("##########")
print("##########")
print("# tickplot per chrom, each bin represented fact with cells")
print("# number of sig inters per bin")
head(r_dat)
#count  each inter twice
tp_A <- r_dat %>% select(chrA,st1,end1,cell,zscore,pvalue) 
colnames(tp_A) <- c("chr","st","end","cell","zscore","pvalue")
tp_B <- r_dat %>% select(chrB,st2,end2,cell,zscore,pvalue) 
colnames(tp_B) <- c("chr","st","end","cell","zscore","pvalue")
#select only sig inters
tp_dat <- rbind(tp_A,tp_B) %>% filter(pvalue <=0.05)
tp_dat_sum = tp_dat %>% group_by(chr,st,cell) %>% dplyr::summarize(numSig=n())
tp_dat_sum$chr <- gsub("chr", "", tp_dat_sum$chr)
tp_dat_sum
tp_dat_sum <- tp_dat_sum %>% mutate(chr=factor(chr, levels=p_chr_ord2))
#scale pts by 1Mb
tp_dat_sum$st <- tp_dat_sum$st / 1000000
summary(tp_dat_sum)
head(tp_dat_sum)
for(i in unique(tp_dat_sum$chr)) {
  #i="X"
  xmax <- chrInf$size[match(paste0("chr",i),chrInf$chrom)] / 1000000
  tpp <- tp_dat_sum %>% filter(chr == i)
  tp <- (ggplot(tpp, aes(st, y=1, fill = numSig))
         + geom_tile(aes(fill = numSig), colour = "white")
         + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Total number of significant interactions", na.value = "grey")
         #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
         #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
         + labs(x = paste0("Chromosome ", i, " position [Mb]"),
                y = "",
                title = "Cis-chromosomal interactions (significant)")
         + facet_grid(cell ~ .)
         #       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
         #         c("nonsig" = "Non-significant", "sig" = "Significant"))))
         #       + theme(axis.text.x = element_text(angle = 90))
         + expand_limits(x = c(0,xmax))
         + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
                 strip.background = element_rect(fill = "white"),
                 panel.spacing = unit(0, "lines"),
                 axis.text.y = element_blank(),
                 axis.ticks.y = element_blank())
         + theme(axis.text=element_text(size=5),panel.background = element_rect(fill = "grey85", colour = NA))
  )
  filename <- paste0("numSig_inters_chrom",i,"_tickplot_sig_Interactions.pdf")
  pdf(filename, width = 14, height = 8)
  print(tp)
  dev.off()
}#for
print("##########")

#################
# classic Hi-C heatmap, on steroids, presence/absence of sig and non-sig
#################

##################
##parallel sets
##################
##prep data for parallel sets plot
#ps_df <- dat2 %>% select(chrA, chrB)
#ps_df <- unique(ps_df)
#ps_df$chrA <- as.factor(ps_df$chrA)
#ps_df$chrB <- as.factor(ps_df$chrB)
##re-order chroms based on chrom len
#ps_df$chrA <- factor(ps_df$chrA, levels=rev_chrs_len_ord)
#ps_df$chrB <- factor(ps_df$chrB, levels=rev_chrs_len_ord)
#ps_df <- as.data.frame(cbind(as.character(ps_df$chrB),as.character(ps_df$chrA)))
#ps_df<- ps_df %>%
#  gather_set_data(1:2)
##ps_df
#
##plot parallel sets
#ps <- (ggplot(data =ps_df, aes(x, id=id, split = y, value = 1))
#       #  + geom_parallel_sets(aes(fill = U00096000))
#       + geom_parallel_sets(alpha = 0.8, fill ="#cbc0d3")
##       + scale_fill_manual(values = c("#2E294E","#BEBEBE"))
#       #  + geom_parallel_sets(aes(fill = U00096 ))
#       + xlab("") 
#       + ylab("")
#       + coord_flip()
##       + scale_x_discrete(expand = c(0,0))
##       + theme(legend.title=element_blank())
#      + geom_parallel_sets_axes(axis.width = 0.1, fill = "grey90", color = "black")
#      + geom_parallel_sets_labels(
#        color = 'black',
##        family = dviz_font_family,
#        size = 10/.pt,
#        angle = 0
#      )
#      + theme(axis.text.x = element_blank(),
#              axis.ticks.x = element_blank(),
#              axis.text.y = element_blank(),
#              axis.ticks.y = element_blank(),
#              axis.line = element_blank(),
#              panel.background = element_rect(color = "white"))
#)
#pdf("parallel_sets_common_interactions_all_cells.pdf", width = 14, height = 8)
#ps
#dev.off()
##################

