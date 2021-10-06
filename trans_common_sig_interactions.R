########################################
#validating known trans-chromosomal interactions
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: 3Dflow output data (tsv)
#            germlayer df (tsv)
#            all interactions file (tsv)
#            tissue/system df (tsv)
########################################

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
germlayer_file <- args[2]
allinters_file <- args[3]
tissue_file <- args[4]
#Atype <- args[4]

##########
library(dplyr)
library(tidyr)
#library(GenomicRanges)
library(ggplot2)
library(ggforce)#for ridgeline
library(ggridges)#for ridgeline
library(ggbiplot)#for PCA
library(devtools)#for PCA
library(factoextra, lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3/")#for PCA
library(harrypotter, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3") #for colours
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
                  legend.text = element_text(size = 18),
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
#interaction data
#Atype <- "1_vs_All"
#tissue_file <- "tissue_system_info.txt"
#dat_file <- "test_pairwise_dat.txt"
#allinters_file <- "all_trans_interactions_1Mb.txt"
#germlayer_file <- "germlayer_info.txt"
#library(factoextra)#for PCA
#library(harrypotter) #for colours

allinters <- read.table(allinters_file, header = FALSE)
colnames(allinters) <- c("chrA", "startA", "endA", "chrB", "startB", "endB")
head(allinters)
dat <- read.table(dat_file, header = TRUE)
print("summary of ALL sig zscores per cell type")
summary(dat)
dat$ID <- as.character(dat$ID)
print("#read in germlayer info file")
print(germlayer_file)
gl_df <- read.table(germlayer_file, header = TRUE, sep="\t")
print("colnames")
colnames(gl_df) <- c("cell","germLayer")
summary(gl_df)
#order of germlayer
gl_ord <- c("ectoderm", "mesoderm", "endoderm", "bipotent", "ectoderm/mesoderm")
#gl_colours <- c("#0D3B66","#F4D35E","#F95738","#66999B","#EE964B")
#darker shades
gl_colours <- c("#071F36","#F2CB40","#ED2E07","#517A7B","#EA7E1F")

#tissue info
ts_df <- read.table(tissue_file,sep = "\t", header = TRUE)
head(ts_df)

#remove interactions involving x and y chrs
#dat <- dat[grep("chrY", df$ID, invert=TRUE), ]
#dat2 <- dat[grep("chrY", df$ID), ]
#dat <- dat[grep("chrX", df$ID, invert=TRUE), ]

#select only rows with NO NAs in any cell type
df <- na.omit(dat)
print("Total number of common interactions")
length(df$ID)
print("Total number of ALL possible interactions")
length(allinters$chrA)
print("percent of ALL common interactions across genome: all possible interactions")
(length(df$ID)/length(allinters$chrA)) *100
print("percent of ALL common interactions across genome: interactions in dataset")
(length(df$ID)/length(dat$ID)) *100

#split ID col
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
df$ID <- sub("B", "\\.B", as.character(df$ID))
dat2 <- df %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
#remove A and B from chrom names
dat2$chrA <- gsub("A", "", dat2$chrA)
dat2$chrB <- gsub("B", "", dat2$chrB)

print("all chroms involved in common interactions")
Achrs <- unique(c(dat2$chrA, dat2$chrB))
Achrs
#################
# PCA 
#################
#PCA with cells as rows
pca_dat <- dat2 %>% select(-chrA, -st1, -end1, -chrB, -st2, -end2)
row.names(pca_dat) <- pca_dat$ID
pca_dat <- pca_dat %>% select(-ID)
pca_dat <- as.data.frame(t(pca_dat))
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
      + labs(title = "Common Trans-chromosomal Interactions")
)
pdf("zscore_PCAcellpts_common_interactions_all_cells.pdf", width = 14, height = 8)
g
dev.off()

#PCA with cells as cols
pca_dat <- as.data.frame(t(pca_dat))
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
+ labs(title = "Trans-chromosomal Interactions")
)
pdf("zscore_PCAcellarrs_common_interactions_all_cells.pdf", width = 14, height = 8)
g
dev.off()

#################

#################
#ridgeline plot of all chroms and number of interactions per genomic region
#################
#wide to long format
r_dat <- gather(dat2, cell, zscore, 8:ncol(dat2), factor_key=TRUE)
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
                      size = c(249300000,243200000,
                               198000000,191200000,
                               180900000,171100000,
                               159100000,155300000,
                               146400000,141200000,
                               135000000,135500000,
                               133900000,115200000,
                               107300000,102500000,
                               90400000,81200000,
                               78100000,63000000,
                               59100000,59400000,
                               51300000,48100000)
  
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
missInters$zscore <- rep(0, length(missInters$ID))
missInters$st1 <- as.numeric(missInters$st1)
missInters$end1 <- as.numeric(missInters$end1)
missInters$st2 <- as.numeric(missInters$st2)
missInters$end2 <- as.numeric(missInters$end2)
missInters$zscore <- as.numeric(missInters$zscore)

r_dat <- rbind(r_dat,missInters)

head(r_dat)
#counting each interaction twice (once for each chrom in interaction)
anchD <- r_dat
anchD$AllChr <- anchD$chrA
anchD$AllSt <- anchD$st1
tarD <- r_dat
tarD$AllChr <- tarD$chrB
tarD$AllSt <- tarD$st2
r_dat2 <- rbind(anchD,tarD)
#scale genomic position by 10Mb
r_dat2$AllSt <- r_dat2$AllSt/10000000
head(r_dat2)
chrInf2 <- chrInf
chrInf2$chrom <- factor(chrInf2$chrom, levels=rev_chrs_len_ord)
chrInf2$centromere <- chrInf2$centromere/10000000
chrInf2$size <- chrInf2$size/10000000
#chrInf2$chrom <-gsub("chr","",chrInf2$chrom)
colnames(chrInf2) <- c("AllChr","centromere","chrClass","size")
#chrInf2 <- chrInf2 %>% filter(AllChr %in% r_dat2$AllChr)
#with fake data which did not really work as expected
#p <- (ggplot(r_dat2, aes(x = AllSt, y = AllChr))
#without fake data
ndat2 <- r_dat2 %>% filter(cell != "fake_cell")
chsindf <- unique(ndat2$AllChr)
red_levels <- factor(chrs_len_ord[chrs_len_ord %in% chsindf])
ndat2$AllChr <- factor(ndat2$AllChr, levels = levels(red_levels))
chrInf_hm <- chrInf2 %>% filter(AllChr %in% chsindf)
chrInf_hm$AllChr <- factor(chrInf_hm$AllChr, levels = levels(red_levels))
p <- (ggplot(ndat2, aes(x = AllSt, y = AllChr))
  + geom_density_ridges(scale = 2, alpha = 0.3, rel_min_height = 0.001) 
#  + geom_segment(data = chrInf2, aes(x=centromere, xend=centromere, 
  + geom_segment(data = chrInf_hm, aes(x=centromere, xend=centromere, 
                                     y=as.numeric(AllChr), yend=as.numeric(AllChr) +0.9),
                color = "red")
#  + geom_segment(data = chrInf2, aes(x=size, xend=size, 
  + geom_segment(data = chrInf_hm, aes(x=size, xend=size, 
                                     y=as.numeric(AllChr), yend=as.numeric(AllChr) +0.9),
                color = "black")
  + labs(x="Genomic Position [10Mb]",
         y="",
         title = "Location of Common Trans-chromosomal Interactions")
  + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
#  + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
  + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
#  + theme(axis.text=element_text(size=12))
)
p
######################################
# NOTE: THE RIDGELINE PLOT HEIGHT CURRENTLY IS SCALED AMONG ALL Y-AXIS POINTS AND IS THE SAME FOR ALL Y-AXIS POINTS. IT DOES NOT REPRESENT ACTUAL AMOUNTS OF INTERACTIONS
# NOTE: THE PER CHROM RIDGELINE PLOT HAS FAKE DATA ADDED SO THAT THE LENGTH OF THE CHROMOSOMES IS ACCURATE, THIS COULD MAKE THE GRAPH LOOK FUNNY! 
######################################
pdf("chr_ridgeline_common_interactions_all_cells.pdf", width = 14, height = 8)
p
dev.off()
##test ridgeline with ACTUAL height representing ACTUAL values of the y-axis
#d <- data.frame(x = rep(1:5, 3), y = c(rep(0, 5), rep(1, 5), rep(3, 5)),
#                height = c(0, 1, 3, 4, 0, 1, 2, 3, 5, 4, 0, 5, 4, 4, 1))
#ggplot(d, aes(x, y, height = height, group = y)) + geom_ridgeline(fill="lightblue")
#
#calculate mean zscore per chrom per bin so that the heatmap is accurate
head(r_dat2)
hm_dat = ndat2 %>% group_by(AllChr, AllSt) %>% dplyr::summarize(mzscore=mean(zscore))
head(hm_dat)
hm_dat$AllChr <- factor(hm_dat$AllChr, levels=rev_chrs_len_ord)
hm_dat$AllChr <- gsub("chr", "", hm_dat$AllChr)
hm <- (ggplot(hm_dat, aes(x=AllSt,y=AllChr, fill = mzscore))
       #hm <- (ggplot(hm_dat, aes(AllChr, cell, fill = zscore))
       #       + geom_tile(aes(fill = mzscore), colour = "white")
       + scale_fill_hp(discrete = FALSE, option = "Always", name = "Mean z-score", na.value = "grey")
       + geom_tile(aes(fill = mzscore), colour = "white")
#       + scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
       + labs(x = "Genomic Postion [Mbp]",
              y = "Chromosome",
              title = "Common Trans-chromosomal Interactions z-scores")
       #       + facet_wrap(.~germL)
       + theme(axis.text=element_text(size=12))
)
pdf("zscore_mean_heatmap_common_interactions_chroms_all_cells.pdf", width = 14, height = 8)
hm
dev.off()



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
p <- (ggplot(r_dat2, aes(x = zscore, y = cell, fill = germL))
      + stat_density_ridges(quantile_lines = TRUE, alpha = 0.3, scale=2, quantiles = 2, rel_min_height = 0.001)
      #+ geom_density_ridges(scale = 4, alpha = 0.3) 
      + labs(x="z-score",
             y="Cell",
             title = "z-scores of Common Trans-chromosomal Interactions",
             fill = "Germ Layer")
      + scale_fill_manual(values = gl_colours)
#                                     gsub("F", "A", my_colors)))
      #  + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
      + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
      + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
  + theme(axis.text=element_text(size=12))
)
pdf("zscore_ridgeline_common_interactions_all_cells_germlayer.pdf", width = 14, height = 8)
p
dev.off()

# Tissue/system breakdown
t_dat2 <- r_dat2
head(t_dat2)
t_dat2$tissue <- ts_df$Tissue.System[match(t_dat2$cell, ts_df$X3Dflow.normalized_data.name)]
t_dat2 <- t_dat2[order(t_dat2$tissue),]
t_cell_ord <- unique(t_dat2$cell)
t_cell_ord
t_dat2$cell <- factor(t_dat2$cell, levels=rev(t_cell_ord))
p <- (ggplot(t_dat2, aes(x = zscore, y = cell, fill = tissue))
      + stat_density_ridges(quantile_lines = TRUE, alpha = 0.3, scale=2, quantiles = 2, rel_min_height = 0.001)
      #+ geom_density_ridges(scale = 4, alpha = 0.3) 
      + labs(x="z-score",
             y="Cell",
             title = "z-scores of Common Trans-chromosomal Interactions",
             fill = "Tissue")
#      + scale_fill_hp_d(option = "NewtScamander", name = "Tissue")
      #                                     gsub("F", "A", my_colors)))
      #  + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
      + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
      + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
  + theme(axis.text=element_text(size=12))
)
pdf("zscore_ridgeline_common_interactions_all_cells_tissue.pdf", width = 14, height = 8)
p
dev.off()
#################

#################
#common interactions z-scores broken down by chrom class
#################
chrClass_dat <- r_dat2
#add metacentric chrom info to df
chrClass_dat$chrClass <- chrClass_dat$AllChr
chrClass_dat$chrClass <- chrInf$chrClass[match(chrClass_dat$AllChr, chrInf$chrom)]
head(chrClass_dat)
p <- (ggplot(chrClass_dat, aes(x = chrClass, y = zscore))
      + geom_violin(fill="grey90", scale = "count")#"count" makes width of violins proportional to number of values
      + geom_boxplot(fill="grey95",width = 0.1,outlier.size=4)
#      + stat_density_ridges(quantile_lines = TRUE, alpha = 0.3, scale=2, quantiles = 2, rel_min_height = 0.001)
#      #+ geom_density_ridges(scale = 4, alpha = 0.3) 
      + labs(y="z-score",
             x="Chromosome Class",
             title = "z-scores of Common Trans-chromosomal Interactions")
#      #  + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
#      + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
#      + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
)
pdf("chrClass_violin_common_interactions_all_cells.pdf", width = 14, height = 8)
p
dev.off()
print("#test between means of chrom classes")
print("#Test each group for normality")
print("sig = reject normality null")
chrClass_dat$chrClass <- as.factor(chrClass_dat$chrClass)
chrClass_dat %>%
  group_by(chrClass) %>%
  summarise(W = shapiro.test(zscore)$statistic,
            p.value = shapiro.test(zscore)$p.value)
print("#Perform the Kruskal-Wallis test")
print("sig = mean is diff btwn groups")
kruskal.test(zscore ~ chrClass, data=chrClass_dat)
#with cell info
#p <- (ggplot(chrClass_dat, aes(y = cell, x = zscore))
#      + geom_violin(fill="grey90")
#      #      + stat_density_ridges(quantile_lines = TRUE, alpha = 0.3, scale=2, quantiles = 2, rel_min_height = 0.001)
#      #      #+ geom_density_ridges(scale = 4, alpha = 0.3) 
#      + facet_grid(.~ chrClass)
#      + labs(x="z-score",
#             y="Cell",
#             title = "z-scores of Common Trans-chromosomal Interactions")
#      #      #  + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
#      #      + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
#      #      + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
#)
p <- (ggplot(chrClass_dat, aes(y = cell, x = zscore, fill = germL))
      + stat_density_ridges(quantile_lines = TRUE, alpha = 0.3, scale=2, quantiles = 2, rel_min_height = 0.001)
      #      + stat_density_ridges(quantile_lines = TRUE, alpha = 0.3, scale=2, quantiles = 2, rel_min_height = 0.001)
      #      #+ geom_density_ridges(scale = 4, alpha = 0.3) 
      + facet_grid(.~ chrClass)
      + labs(x="z-score",
             y="Cell",
             title = "z-scores of Common Trans-chromosomal Interactions",
             fill = "Germ Layer")
      + scale_fill_manual(values = gl_colours)
      #      #  + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
      #      + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
      #      + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
  + theme(axis.text=element_text(size=12))
)
pdf("chrClass_cell_facet_ridgeline_common_interactions_all_cells_germlayer.pdf", width = 14, height = 8)
p
dev.off()
# Tissue/system
head(chrClass_dat)
t_dat2 <- chrClass_dat
t_dat2$tissue <- ts_df$Tissue.System[match(t_dat2$cell, ts_df$X3Dflow.normalized_data.name)]
t_dat2 <- t_dat2[order(t_dat2$tissue),]
t_cell_ord <- unique(t_dat2$cell)
t_dat2$cell <- factor(t_dat2$cell, levels=rev(t_cell_ord))
p <- (ggplot(t_dat2, aes(y = cell, x = zscore, fill = tissue))
      + stat_density_ridges(quantile_lines = TRUE, alpha = 0.3, scale=2, quantiles = 2, rel_min_height = 0.001)
      #      + stat_density_ridges(quantile_lines = TRUE, alpha = 0.3, scale=2, quantiles = 2, rel_min_height = 0.001)
      #      #+ geom_density_ridges(scale = 4, alpha = 0.3) 
      + facet_grid(.~ chrClass)
      + labs(x="z-score",
             y="Cell",
             title = "z-scores of Common Trans-chromosomal Interactions",
             fill = "Tissue")
#      + scale_fill_manual(values = gl_colours)
      #      #  + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
      #      + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
      #      + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
  + theme(axis.text=element_text(size=12))
)
pdf("chrClass_cell_facet_ridgeline_common_interactions_all_cells_tissue.pdf", width = 14, height = 8)
p
dev.off()
#################

#################
#proportional interactions per chrom
#################
prop_dat <- r_dat2
prop_dat %>% filter(ID == "Achr10.0.1000000.Bchr18.79000000.80000000")
##re-format all interactions df
##split ID col
#dat$ID <- sub("B", "\\.B", as.character(dat$ID))
#dat <- dat %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
##remove A and B from chrom names
#dat$chrA <- gsub("A", "", dat$chrA)
#dat$chrB <- gsub("B", "", dat$chrB)
#dat_long <- gather(dat, cell, zscore, 8:ncol(dat2), factor_key=TRUE)
##counting each interaction twice (once for each chrom in interaction)
#anchD <- dat_long
#anchD$AllChr <- anchD$chrA
#anchD$AllSt <- anchD$st1
#tarD <- dat_long
#tarD$AllChr <- tarD$chrB
#tarD$AllSt <- tarD$st2
#dat_long2 <- rbind(anchD,tarD)
##scale genomic position by 10Mb
#dat_long2$AllSt <- as.numeric(as.character(dat_long2$AllSt))/10000000
##df with total number of interactions per chrom
#totInter <- dat_long2 %>%
#  select(AllChr, ID) %>%
#  group_by(AllChr) %>%
#  dplyr::summarise(n = n())
#totInter$AllChr <- factor(totInter$AllChr, levels=rev_chrs_len_ord)

#re-format all interactions df
head(allinters)
#counting each interaction twice (once for each chrom in interaction)
anchD <- allinters
anchD$AllChr <- anchD$chrA
anchD$AllSt <- anchD$startA
tarD <- allinters
tarD$AllChr <- tarD$chrB
tarD$AllSt <- tarD$startB
dat_long2 <- rbind(anchD,tarD)
#head(dat_long2)
#tail(dat_long2)
#allinters %>% filter(ID == "Achr1.0.1000000.Bchr2.0.1000000")
#dat_long2 %>% filter(ID == "Achr1.0.1000000.Bchr2.0.1000000")
#allinters %>% filter(ID == "Achr22.50000000.50818468.Bchr21.44000000.45000000")
#dat_long2 %>% filter(ID == "Achr22.50000000.50818468.Bchr21.44000000.45000000")
#scale genomic position by 10Mb
dat_long2$AllSt <- as.numeric(as.character(dat_long2$AllSt))/10000000
#df with total number of interactions per chrom
totInter <- dat_long2 %>%
  select(AllChr, ID) %>%
  group_by(AllChr) %>%
  dplyr::summarise(n = n())
totInter$AllChr <- factor(totInter$AllChr, levels=rev_chrs_len_ord)
totInter
#df with total number of common interactions per chrom (ALL POSSIBLE INTERACTIONS)
#remove cell info (which makes interactions duplicated)
prop_dat2 <- prop_dat %>% select(AllChr, ID)
prop_dat2 %>% filter(ID == "Achr10.0.1000000.Bchr13.113000000.114000000")
prop_dat2 <- unique(prop_dat2)
prop_dat2 %>% filter(ID == "Achr10.0.1000000.Bchr13.113000000.114000000")
commonInter <- prop_dat2 %>%
  select(AllChr, ID) %>%
  group_by(AllChr) %>%
  dplyr::summarise(n = n())
commonInter$AllChr <- factor(commonInter$AllChr, levels=rev_chrs_len_ord)
#df combining above 2 dfs
prop_datAll <- totInter
prop_datAll$commonInter <- prop_datAll$n
colnames(prop_datAll) <- c("chrom", "totInter","commonInter")
prop_datAll$commonInter <- as.numeric(commonInter$n[match(prop_datAll$chrom, commonInter$AllChr)])
prop_datAll$commonInter[is.na(prop_datAll$commonInter)] <- 0
prop_datAll$percent <- (prop_datAll$commonInter/prop_datAll$totInter)*100 
prop_datAll$chrom <- as.factor(prop_datAll$chrom)
prop_datAll$chrom <- as.factor(gsub("chr","", prop_datAll$chrom))
chrs_len_ord_num <- gsub("chr","",chrs_len_ord)
prop_datAll$chrom <- factor(prop_datAll$chrom, levels=chrs_len_ord_num)
levels(prop_datAll$chrom)
##plot
#p <- (ggplot(prop_datAll, aes(y = percent,x=chrom))
#      + geom_bar(fill="grey90", color="black", stat = "identity")
#      + labs(y="Percent of Total Interactions [%]",
#             x="Chromosome",
#             title = "Percentage of Common Trans-chromosomal Interactions per Chromosome")
#)
#pdf("proportion_per_chrom_common_interactions_all_cells.pdf", width = 14, height = 8)
#p
#dev.off()
prop_datAll2 <- prop_datAll %>% select(chrom, totInter, commonInter) %>% gather(key = "name", value = "num", 2:3)
head(prop_datAll2)
prop_datAll2$name <- factor(prop_datAll2$name, levels = c("totInter", "commonInter"))
#scaling number of interactions by 100,000
prop_datAll2$num <- prop_datAll2$num / 100000
#plot
p <- (ggplot(prop_datAll2, aes(y =num,x=chrom, fill = name))
      + geom_bar(color="black",position = "dodge", stat = "identity")
      + scale_fill_manual(values =c("#FAC9A1", "#013040"), labels= c("Total Interactions", "Common Interactions"))
      + labs(y="Number of Interactions [100,000]",
             x="Chromosome",
             title = "Trans-chromosomal Interactions per Chromosome",
             fill = "")
)
pdf("proportion_per_chrom_common_interactions_all_cells.pdf", width = 14, height = 8)
p
dev.off()

#proportional interactions per chrom pair
pair_dat <- dat_long2
head(pair_dat)
#new col for chrom pair
pair_dat$pair <- paste0(pair_dat$chrA,pair_dat$chrB)
pair_dat <- pair_dat %>% select(ID, pair) %>% unique()
#df with total number of interactions per chrom pair ALL POSSIBLE INTERACTIONS
totInter <- pair_dat %>%
  select(pair, ID) %>%
  group_by(pair) %>%
  dplyr::summarise(n = n())
#df with total number of common interactions per chrom pair
commonInter <- prop_dat
head(prop_dat)
#commonInter <- gather(dat2, cell, zscore, 8:ncol(dat2), factor_key=TRUE)
#new col for chrom pair
commonInter$pair <- paste0(commonInter$chrA,commonInter$chrB)
commonInter <- commonInter %>% select(ID, pair) %>% unique()
head(commonInter)
commonInter %>% filter(ID == "Achr10.0.1000000.Bchr13.113000000.114000000")
commonInter <- commonInter %>%
  select(pair, ID) %>%
  group_by(pair) %>%
  dplyr::summarise(n = n())
commonInter
#df combining above 2 dfs
prop_pairAll <- totInter
prop_pairAll$commonInter <- prop_pairAll$n
colnames(prop_pairAll) <- c("chrom", "totInter","commonInter")
prop_pairAll$commonInter <- as.numeric(commonInter$n[match(prop_pairAll$chrom, commonInter$pair)])
prop_pairAll$commonInter[is.na(prop_pairAll$commonInter)] <- 0
prop_pairAll$percent <- (prop_pairAll$commonInter/prop_pairAll$totInter)*100 
prop_pairAll$chrom <- as.factor(prop_pairAll$chrom)
print("# chrom pairs proportional number of common interactions")
prop_pairAll
print("# chromosome pair(s) with highest (proportional) number of common interactions")
prop_pairAll[which(prop_pairAll$percent == max(prop_pairAll$percent)),]
#################

#################
#heatmap of common interactions z-scores by cell type and chromosome
#################
r_dat2$AllChr <- gsub("chr","",r_dat2$AllChr)
r_dat2$AllChr <- as.factor(r_dat2$AllChr)
chrs_ord <- gsub("chr", "", chrs_len_ord)
r_dat2$AllChr <- factor(r_dat2$AllChr, levels=chrs_ord)
#calculate mean zscore per chrom per cell type so heat map is accutate
hm_dat = r_dat2 %>% group_by(AllChr,cell) %>% dplyr::summarize(mzscore=mean(zscore))
hm <- (ggplot(r_dat2, aes(AllChr, cell))
#hm <- (ggplot(hm_dat, aes(AllChr, cell, fill = zscore))
#       + geom_tile(aes(fill = mzscore), colour = "white")
       + geom_tile(aes(fill = zscore), colour = "white")
       + scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
       + labs(x = "Chromosome",
              y = "Cell",
              title = "Common Trans-chromosomal Interactions z-scores")
#       + facet_wrap(.~germL)
       + theme(axis.text=element_text(size=12))
)
pdf("zscore_mean_heatmap_common_interactions_chroms_all_cells.pdf", width = 14, height = 8)
hm
dev.off()
#################
# average z-score per chrom pair heatmap 
#################
head(r_dat)
hm_dat <- r_dat
hm_dat$chrPair <- paste0(hm_dat$chrA,hm_dat$chrB)
head(hm_dat)
print("chrom pairs involved in common interactions")
unique(hm_dat$chrPair)
#calculate mean zscore per chrom pair per cell type so heat map is accutate
hm_dat2 = hm_dat %>% filter(cell != "fake_cell") %>% group_by(chrPair,cell) %>% dplyr::summarize(mzscore=mean(zscore, na.rm = TRUE))
hm_dat2$chrPair <- gsub("chr", "\\.chr", hm_dat2$chrPair)
coln <- c("tmp","chrA", "chrB")
hm_dat2 <- hm_dat2 %>% separate(chrPair, sep = "\\.", into = coln, remove = FALSE) %>% select(-tmp)
hm_dat2$chrA <- gsub("chr", "", hm_dat2$chrA)
hm_dat2$chrB <- gsub("chr", "", hm_dat2$chrB)
hm_dat2$chrA <- factor(hm_dat2$chrA, levels=chrs_ord)
hm_dat2$chrB <- factor(hm_dat2$chrB, levels=chrs_ord)
head(hm_dat2)
hm <- (ggplot(hm_dat2, aes(chrA, chrB, fill = mzscore))
       + geom_tile(aes(fill = mzscore), colour = "white")
       + scale_fill_hp(discrete = FALSE, option = "Always", name = "Mean z-score per chromosomal pair", na.value = "grey")
       #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
       #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
       + labs(x = "Chromosome",
              y = "",
              title = "Common Trans-chromosomal Interactions z-scores")
       + facet_wrap(cell ~ .)
       #       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
       #         c("nonsig" = "Non-significant", "sig" = "Significant"))))
       #       + theme(axis.text.x = element_text(angle = 90))
       + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
               strip.background = element_rect(fill = "white"),
               panel.spacing = unit(0, "lines"),
#               axis.text.y = element_blank(),
#               axis.ticks.y = element_blank(),
               axis.text=element_text(size=10))
)
pdf("zscore_chrom_pair_mean_heatmap_common_interactions.pdf", width = 14, height = 8)
hm
dev.off()



#################
#parallel sets
#################
#prep data for parallel sets plot
ps_df <- dat2 %>% select(chrA, chrB)
ps_df <- unique(ps_df)
ps_df$chrA <- as.factor(ps_df$chrA)
ps_df$chrB <- as.factor(ps_df$chrB)
#re-order chroms based on chrom len
ps_df$chrA <- factor(ps_df$chrA, levels=rev_chrs_len_ord)
ps_df$chrB <- factor(ps_df$chrB, levels=rev_chrs_len_ord)
ps_df <- as.data.frame(cbind(as.character(ps_df$chrB),as.character(ps_df$chrA)))
ps_df<- ps_df %>%
  gather_set_data(1:2)
#ps_df

#plot parallel sets
ps <- (ggplot(data =ps_df, aes(x, id=id, split = y, value = 1))
       #  + geom_parallel_sets(aes(fill = U00096000))
       + geom_parallel_sets(alpha = 0.8, fill ="#cbc0d3")
#       + scale_fill_manual(values = c("#2E294E","#BEBEBE"))
       #  + geom_parallel_sets(aes(fill = U00096 ))
       + xlab("") 
       + ylab("")
       + coord_flip()
#       + scale_x_discrete(expand = c(0,0))
#       + theme(legend.title=element_blank())
      + geom_parallel_sets_axes(axis.width = 0.1, fill = "grey90", color = "black")
      + geom_parallel_sets_labels(
        color = 'black',
#        family = dviz_font_family,
        size = 10/.pt,
        angle = 0
      )
      + theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line = element_blank(),
              panel.background = element_rect(color = "white"))
)
pdf("parallel_sets_common_interactions_all_cells.pdf", width = 14, height = 8)
ps
dev.off()
#################

