########################################
#validating known trans-chromosomal interactions
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: 3Dflow output data (tsv)
#            germlayer df (tsv)
########################################

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
germlayer_file <- args[2]
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
library(harrypotter)
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
dat <- read.table("test_pairwise_dat.txt", header = TRUE)
dat <- read.table(dat_file, header = TRUE)
dat <- as.data.frame(dat)
print("summary of ALL sig zscores per cell type")
summary(dat)
dat$ID <- as.character(dat$ID)
#read in germlayer info file
gl_df <- read.table("germlayer_info.txt",sep = "\t", header = TRUE)
gl_df <- read.table(germlayer_file, header = TRUE)
colnames(gl_df) <- c("cell","germLayer")
summary(gl_df)
#order of germlayer
gl_ord <- c("ectoderm", "mesoderm", "endoderm", "bipotent", "ectoderm/mesoderm")
#gl_colours <- c("#0D3B66","#F4D35E","#F95738","#66999B","#EE964B")
#darker shades
gl_colours <- c("#071F36","#F2CB40","#ED2E07","#517A7B","#EA7E1F")

#remove interactions involving x and y chrs
#dat <- dat[grep("chrY", df$ID, invert=TRUE), ]
#dat2 <- dat[grep("chrY", df$ID), ]
#dat <- dat[grep("chrX", df$ID, invert=TRUE), ]

#select only rows with NO NAs in any cell type
df <- na.omit(dat)
print("percent of ALL common interactions across genome")
(length(df$ID)/length(dat$ID)) *100
head(df)



#split ID col
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
df$ID <- sub("B", "\\.B", as.character(df$ID))
dat2 <- df %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
#remove A and B from chrom names
dat2$chrA <- gsub("A", "", dat2$chrA)
dat2$chrB <- gsub("B", "", dat2$chrB)

#################
# PCA 
#################
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
pdf("zscore_PCA_common_interactions_all_cells.pdf", width = 14, height = 8)
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
#chrInf2$chrom <-gsub("chr","",chrInf2$chrom)
colnames(chrInf2) <- c("AllChr","centromere","chrClass")
#chrInf2 <- chrInf2 %>% filter(AllChr %in% r_dat2$AllChr)
p <- (ggplot(r_dat2, aes(x = AllSt, y = AllChr))
  + geom_density_ridges(scale = 2, alpha = 0.3) 
  + geom_segment(data = chrInf2, aes(x=centromere, xend=centromere, 
                                     y=as.numeric(AllChr), yend=as.numeric(AllChr) +0.9),
                color = "red")
  + labs(x="Genomic Position [10Mb]",
         y="",
         title = "Location of Common Trans-chromosomal Interactions")
#  + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
  + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
  + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
)
######################################
# NOTE: THE RIDGELINE PLOT HEIGHT CURRENTLY IS SCALED AMONG ALL Y-AXIS POINTS AND IS THE SAME FOR ALL Y-AXIS POINTS. IT DOES NOT REPRESENT ACTUAL AMOUNTS OF INTERACTIONS
######################################
pdf("chr_ridgeline_common_interactions_all_cells.pdf", width = 14, height = 8)
p
dev.off()
##test ridgeline with ACTUAL height representing ACTUAL values of the y-axis
#d <- data.frame(x = rep(1:5, 3), y = c(rep(0, 5), rep(1, 5), rep(3, 5)),
#                height = c(0, 1, 3, 4, 0, 1, 2, 3, 5, 4, 0, 5, 4, 4, 1))
#ggplot(d, aes(x, y, height = height, group = y)) + geom_ridgeline(fill="lightblue")
#




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
)
pdf("zscore_ridgeline_common_interactions_all_cells.pdf", width = 14, height = 8)
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
)
pdf("chrClass_cell_facet_ridgeline_common_interactions_all_cells.pdf", width = 14, height = 8)
p
dev.off()
#################

#################
#proportional interactions per chrom
#################
prop_dat <- r_dat2
#re-format all interactions df
#split ID col
dat$ID <- sub("B", "\\.B", as.character(dat$ID))
dat <- dat %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
#remove A and B from chrom names
dat$chrA <- gsub("A", "", dat$chrA)
dat$chrB <- gsub("B", "", dat$chrB)
dat_long <- gather(dat, cell, zscore, 8:ncol(dat2), factor_key=TRUE)
#counting each interaction twice (once for each chrom in interaction)
anchD <- dat_long
anchD$AllChr <- anchD$chrA
anchD$AllSt <- anchD$st1
tarD <- dat_long
tarD$AllChr <- tarD$chrB
tarD$AllSt <- tarD$st2
dat_long2 <- rbind(anchD,tarD)
#scale genomic position by 10Mb
dat_long2$AllSt <- as.numeric(as.character(dat_long2$AllSt))/10000000
#df with total number of interactions per chrom
totInter <- dat_long2 %>%
  select(AllChr, ID) %>%
  group_by(AllChr) %>%
  dplyr::summarise(n = n())
totInter$AllChr <- factor(totInter$AllChr, levels=rev_chrs_len_ord)
#df with total number of common interactions per chrom
commonInter <- prop_dat %>%
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
#plot
p <- (ggplot(prop_datAll, aes(y = percent,x=chrom))
      + geom_bar(fill="grey90", color="black", stat = "identity")
      + labs(y="Percent of Total Interactions [%]",
             x="Chromosome",
             title = "Percentage of Common Trans-chromosomal Interactions per Chromosome")
)
pdf("proportion_per_chrom_common_interactions_all_cells.pdf", width = 14, height = 8)
p
dev.off()

#proportional interactions per chrom pair
pair_dat <- dat_long
#new col for chrom pair
pair_dat$pair <- paste0(pair_dat$chrA,pair_dat$chrB)
#df with total number of interactions per chrom pair
totInter <- pair_dat %>%
  select(pair, ID) %>%
  group_by(pair) %>%
  dplyr::summarise(n = n())
#df with total number of common interactions per chrom pair
commonInter <- gather(dat2, cell, zscore, 8:ncol(dat2), factor_key=TRUE)
#new col for chrom pair
commonInter$pair <- paste0(commonInter$chrA,commonInter$chrB)
commonInter <- commonInter %>%
  select(pair, ID) %>%
  group_by(pair) %>%
  dplyr::summarise(n = n())
#df combining above 2 dfs
prop_pairAll <- totInter
prop_pairAll$commonInter <- prop_pairAll$n
colnames(prop_pairAll) <- c("chrom", "totInter","commonInter")
prop_pairAll$commonInter <- as.numeric(commonInter$n[match(prop_pairAll$chrom, commonInter$pair)])
prop_pairAll$commonInter[is.na(prop_pairAll$commonInter)] <- 0
prop_pairAll$percent <- (prop_pairAll$commonInter/prop_pairAll$totInter)*100 
prop_pairAll$chrom <- as.factor(prop_pairAll$chrom)
print("# chromosome pair(s) with highest (proportional) number of interactions")
prop_pairAll[which(prop_pairAll$percent == max(prop_pairAll$percent)),]
#################

#################
#heatmap of common interactions z-scores by cell type and chromosome
#################
r_dat2$AllChr <- gsub("chr","",r_dat2$AllChr)
r_dat2$AllChr <- as.factor(r_dat2$AllChr)
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
#       + theme(axis.text.x = element_text(angle = 90))
)
pdf("zscore_heatmap_common_interactions_chroms_all_cells.pdf", width = 14, height = 8)
hm
dev.off()
#################

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

