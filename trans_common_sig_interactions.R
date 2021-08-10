########################################
#validating known trans-chromosomal interactions
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: 3Dflow output data
########################################

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
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

# PCA 
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



#ridgeline plot of all chroms and number of interactions per genomic region
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
p <- (ggplot(r_dat2, aes(x = AllSt, y = AllChr))
  + stat_density_ridges(scale = 2, alpha = 0.3,) 
  + labs(x="Genomic Position [10Mb]",
         y="",
         title = "Location of Common Trans-chromosomal Interactions")
#  + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
  + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
  + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
)
pdf("chr_ridgeline_common_interactions_all_cells.pdf", width = 14, height = 8)
p
dev.off()

#ridgeline of all zscores (in all chroms) across cell types
head(r_dat2)
p <- (ggplot(r_dat2, aes(x = zscore, y = cell))
      + stat_density_ridges(quantile_lines = TRUE, alpha = 0.3, scale=2, quantiles = 2, rel_min_height = 0.001)
      #+ geom_density_ridges(scale = 4, alpha = 0.3) 
      + labs(x="z-score",
             y="Cell",
             title = "z-scores of Common Trans-chromosomal Interactions")
      #  + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
      + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
      + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
)
pdf("zscore_ridgeline_common_interactions_all_cells.pdf", width = 14, height = 8)
p
dev.off()

#heatmap of common interactions z-scores by cell type and chromosome
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
#       + theme(axis.text.x = element_text(angle = 90))
)
pdf("zscore_heatmap_common_interactions_chroms_all_cells.pdf", width = 14, height = 8)
hm
dev.off()

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


