########################################
#Originally used to look replicates in 3Dflow pipeline. I think a simple correlation of the zscores or interaction frequencies would be better (seems like what others are doing).
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: 3Dflow z-score output data (tsv)
#            3Dflow p-value output data (tsv)
#            output path
########################################

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
zscoreFile <- args[1]
pvaluFile <- args[2]
output <- args[3]

##########
library(dplyr)
library(tidyr)
library(VIM, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3")#for KNN imputation
library(class)#for KNN imputation (picking k)
#library(GenomicRanges)
library(ggplot2)
library(ggforce)#for ridgeline
library(ggridges)#for ridgeline
library(ggbiplot, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3")#for PCA
library(devtools)#for PCA
library(inauguration, lib="/hpf/largeprojects/pmaass/Daniella/R-lib")
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
#library(inauguration)
#library(VIM)
#library(ggbiplot)
#zscoreFile <- "test_1vsAll_dat.txt"
#pvaluFile <- "test_1vsAll_pvalues.txt"
datz <- read.table(zscoreFile, header = TRUE)
#datz <- as.data.frame(datz)
datp <- read.table(pvaluFile, header = TRUE)
#datp <- as.data.frame(datp)
head(datp)
head(datz)
#combine zscore and pvals into one df
datzL <- datz %>% gather(key="cell", value = "zscore", 2:length(colnames(datz)))
datpL <- datp %>% gather(key="cell", value = "pvalue", 2:length(colnames(datp)))
head(datzL)
head(datpL)
dat <- merge(datzL, datpL, by = c("ID", "cell"))
print("summary of ALL sig zscores per cell type")
summary(dat)
dat$ID <- as.character(dat$ID)
#add in column to indicate which reps are part of the same cell
dat$cell_noreps <- as.character(dat$cell)
dat$cell_noreps[grepl("Cardiomyocites_Progenitors",dat$cell)]<-as.character("Cardiomyocites_Progenitors")
dat$cell_noreps[grepl("Cardiomyocites_primitive",dat$cell)]<-"Cardiomyocites_Primitive"
dat$cell_noreps[grepl("Cardiomyocites_ventricular",dat$cell)]<-"Cardiomyocites_Ventricular"
dat$cell_noreps[grepl("Astrocytes_cerebellum",dat$cell)]<-"Astrocytes_Cerebellum"
head(dat)
##read in germlayer info file
#gl_df <- read.table("germlayer_info.txt",sep = "\t", header = TRUE)
#gl_df <- read.table(germlayer_file, header = TRUE)
#colnames(gl_df) <- c("cell","germLayer")
#summary(gl_df)
##order of germlayer
#gl_ord <- c("ectoderm", "mesoderm", "endoderm", "bipotent", "ectoderm/mesoderm")
##gl_colours <- c("#0D3B66","#F4D35E","#F95738","#66999B","#EE964B")
##darker shades
#gl_colours <- c("#071F36","#F2CB40","#ED2E07","#517A7B","#EA7E1F")

#remove interactions involving x and y chrs
#dat <- dat[grep("chrY", df$ID, invert=TRUE), ]
#dat2 <- dat[grep("chrY", df$ID), ]
#dat <- dat[grep("chrX", df$ID, invert=TRUE), ]

##select only rows with NO NAs in any cell type
#df <- na.omit(dat)
#print("percent of ALL common interactions across genome")
#(length(df$ID)/length(dat$ID)) *100
#head(df)



#split ID col
df <- dat
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
df$ID <- sub("B", "\\.B", as.character(df$ID))
dat2 <- df %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
#remove A and B from chrom names
dat2$chrA <- gsub("A", "", dat2$chrA)
dat2$chrB <- gsub("B", "", dat2$chrB)
head(dat2)
dat2$cell <- as.factor(dat2$cell)
dat2$cell_noreps <- as.factor(dat2$cell_noreps)
#NAs are a problem with PCA. need to be removed or dealt with
#https://www.edureka.co/blog/knn-algorithm-in-r/
#print("#using KNN imputation to fill in the missing values (NAs)")
##optimizing value for k (using zscore col)
#set.seed(123)
#datsimple <- dat2 %>% select(zscore,pvalue)
##datsimple <- dat2 %>% select(pvalue,zscore)
#datsimple <- na.omit(datsimple)
#datsub <- dat2 %>% select(zscore,pvalue) %>% select(-zscore)
##datsub <- dat2 %>% select(pvalue,zscore) %>% select(-pvalue)
#datsub <- na.omit(datsub)
#dat.d <- sample(1:nrow(datsub),size=nrow(datsub)*0.7,replace = FALSE) #random selection of 70% data.
#train.loan <- datsimple[dat.d,] # 70% training data
#test.loan <- datsimple[-dat.d,] # remaining 30% test data
##Creating seperate dataframe for 'zscore' feature which is our target.
#train.loan_labels <- datsimple[dat.d,1]
#test.loan_labels <-datsimple[-dat.d,1]
#i=1
#k.optm=1
#k.optmvec <- c()
#for (i in 1:100){
#   knn.mod <- as.numeric(as.character(class::knn(train=train.loan, test=test.loan, cl=train.loan_labels, k=i)))
#   k.optm[i] <- 100 * sum(test.loan_labels == knn.mod)/NROW(test.loan_labels)
#   k=i
#   k.optmvec <- c(k.optmvec, k.optm[i])
#   #cat(k,'=',k.optm[i],'')
#}
#optK <- which(k.optmvec == max(k.optmvec))
#summary(dat2)
#dat.knn <- VIM::kNN(dat2, variable = c("zscore", "pvalue"), k=50)
##dat.knn <- VIM::kNN(dat2, variable = c("zscore", "pvalue"), k=optK)
#dat.knn <- dat.knn %>% select(-zscore_imp, -pvalue_imp)
#summary(dat.knn)
#
##################
#print("# PCA with z-score and p-value as cols ")
##################
##NAs are a problem with PCA. need to be removed.
#pca_dat <- dat.knn %>% select(-chrA, -st1, -end1, -chrB, -st2, -end2, -ID, -cell)
#head(pca_dat)
##row.names(pca_dat) <- pca_dat$ID
##pca_dat <- pca_dat %>% select(-ID)
##pca_dat <- as.data.frame(t(pca_dat))
#commonInter.pca <- prcomp(pca_dat[,-3], center = TRUE,scale. = TRUE)
#summary(commonInter.pca)
#print("graph 1")
#g <- (ggbiplot(commonInter.pca,
#              obs.scale = 1,
##              var.axes=FALSE,
#              var.scale = 1,
##              labels = row.names(pca_dat),
#              groups = dat2$cell,
##              ellipse = TRUE,
#              circle = TRUE,
##              ellipse.prob = 0.68
#      )
#      + labs(title = "Trans-chromosomal Interactions")
#)
#pdf("zscore_PCA_trans_interactions_zscore_pvalue_replicates.pdf", width = 14, height = 8)
#g
#dev.off()
##################
##################
#print("# PCA with cells as rows (zscore only) ")
##################
#datW <- dat.knn %>% select(-pvalue, -cell_noreps) %>% spread(key = "cell", "zscore")
#pca_dat <- datW %>% select(-chrA, -st1, -end1, -chrB, -st2, -end2)
#head(pca_dat)
#row.names(pca_dat) <- pca_dat$ID
#pca_dat <- pca_dat %>% select(-ID)
#pca_dat <- as.data.frame(t(pca_dat))
#head(pca_dat)
#commonInter.pca <- prcomp(pca_dat, center = TRUE,scale. = TRUE)
#summary(commonInter.pca)
#print("graph 2")
#g <- (ggbiplot(commonInter.pca,
#              obs.scale = 1,
#              var.axes=FALSE,
#              var.scale = 1,
#              labels = row.names(pca_dat),
##              groups = colnames(pca_dat),
#              ellipse = TRUE,
#              circle = TRUE,
#              ellipse.prob = 0.68
#      )
#      + labs(title = "Trans-chromosomal Interactions")
#)
#pdf("zscore_PCA_common_interactions_zscore_cells_as_cols_replicates.pdf", width = 14, height = 8)
#g
#dev.off()
#################
# REMOVING NAs PCAs
#################
dat3 <- na.omit(dat2)
#################
print("# REMOVING NAs, PCA with z-score and p-value as cols ")
#################
#NAs are a problem with PCA. need to be removed.
pca_dat <- dat3 %>% select(-chrA, -st1, -end1, -chrB, -st2, -end2, -ID, -cell)
head(pca_dat)
#row.names(pca_dat) <- pca_dat$ID
#pca_dat <- pca_dat %>% select(-ID)
#pca_dat <- as.data.frame(t(pca_dat))
commonInter.pca <- prcomp(pca_dat[,-3], center = TRUE,scale. = TRUE)
summary(commonInter.pca)
g <- (ggbiplot(commonInter.pca,
               obs.scale = 1,
               #              var.axes=FALSE,
               var.scale = 1,
               #              labels = row.names(pca_dat),
               groups = dat3$cell,
               ellipse = TRUE,
               circle = TRUE,
               ellipse.prob = 0.68
)
+ labs(title = "Trans-chromosomal Interactions")
)
pdf("zscore_PCA_trans_interactions_zscore_pvalue_replicates_NArm.pdf", width = 14, height = 8)
g
dev.off()
##################
##Had to remove this section because when removing NAs for the PCA there are no points left (i.e. no interactions that are sig in all the cells
#################
print("# REMOVING NAs, PCA with cells as rows (zscore only) ")
#################
datW <- dat2 %>% select(-pvalue, -cell_noreps) %>% spread(key = "cell", "zscore")
head(datW)
datW <- na.omit(datW)
pca_dat <- datW %>% select(-chrA, -st1, -end1, -chrB, -st2, -end2)
print("head pca_dat")
head(pca_dat)
row.names(pca_dat) <- pca_dat$ID
pca_dat <- pca_dat %>% select(-ID)
pca_dat <- as.data.frame(t(pca_dat))
head(pca_dat)
commonInter.pca <- prcomp(pca_dat, center = TRUE,scale. = TRUE)
summary(commonInter.pca)
g <- (ggbiplot(commonInter.pca,
               obs.scale = 1,
               var.axes=FALSE,
               var.scale = 1,
               labels = row.names(pca_dat),
               #              groups = colnames(pca_dat),
               ellipse = TRUE,
               circle = TRUE,
               ellipse.prob = 0.68
)
+ labs(title = "Trans-chromosomal Interactions")
)
pdf("zscore_PCA_common_interactions_zscore_cells_as_cols_replicates_NArm.pdf", width = 14, height = 8)
g
dev.off()
#################




#################
#ridgeline plot of all chroms and number of interactions per genomic region
#################
#wide to long format
head(dat2)
r_dat <- dat2
#r_dat <- gather(dat2, cell, zscore, 8:ncol(dat2), factor_key=TRUE)
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
##chrInf2 <- chrInf2 %>% filter(AllChr %in% r_dat2$AllChr)
#p <- (ggplot(r_dat2, aes(x = AllSt, y = AllChr))
#  + geom_density_ridges(scale = 2, alpha = 0.3) 
#  + geom_segment(data = chrInf2, aes(x=centromere, xend=centromere, 
#                                     y=as.numeric(AllChr), yend=as.numeric(AllChr) +0.9),
#                color = "red")
#  + labs(x="Genomic Position [10Mb]",
#         y="",
#         title = "Location of Common Trans-chromosomal Interactions")
##  + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
#  + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
#  + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
#)
#######################################
## NOTE: THE RIDGELINE PLOT HEIGHT CURRENTLY IS SCALED AMONG ALL Y-AXIS POINTS AND IS THE SAME FOR ALL Y-AXIS POINTS. IT DOES NOT REPRESENT ACTUAL AMOUNTS OF INTERACTIONS
#######################################
#pdf("chr_ridgeline_common_interactions_all_cells.pdf", width = 14, height = 8)
#p
#dev.off()
###test ridgeline with ACTUAL height representing ACTUAL values of the y-axis
##d <- data.frame(x = rep(1:5, 3), y = c(rep(0, 5), rep(1, 5), rep(3, 5)),
##                height = c(0, 1, 3, 4, 0, 1, 2, 3, 5, 4, 0, 5, 4, 4, 1))
##ggplot(d, aes(x, y, height = height, group = y)) + geom_ridgeline(fill="lightblue")
##




##ridgeline of all zscores (in all chroms) across cell types
##adding germlayer info
#r_dat2$germL <- r_dat2$cell
#r_dat2$germL <- as.factor(gl_df$germLayer[match(r_dat2$cell, gl_df$cell)])
##re-order based on gl_ord
#r_dat2$germL <- factor(r_dat2$germL, levels=gl_ord)
#r_dat2 <- r_dat2[order(r_dat2$germL),]
#gl_cell_ord <- unique(r_dat2$cell)
#gl_cell_ord
#r_dat2$cell <- factor(r_dat2$cell, levels=rev(gl_cell_ord))
##r_dat2 <- r_dat2 %>% mutate(cell = factor(cell,levels=cell))
#levels(r_dat2$germL)
#levels(r_dat2$cell)
head(r_dat2)
p <- (ggplot(r_dat2, aes(x = zscore, y = cell, fill = cell_noreps))
      + stat_density_ridges(quantile_lines = TRUE, alpha = 0.3, scale=2, quantiles = 2, rel_min_height = 0.001)
      #+ geom_density_ridges(scale = 4, alpha = 0.3) 
      + labs(x="z-score",
             y="Cell",
             title = "z-scores of Trans-chromosomal Interactions",
             fill = "")
      + scale_fill_manual(values = inauguration("inauguration_2021"))
#                                     gsub("F", "A", my_colors)))
      #  + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
      + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
      + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
)
pdf("zscore_ridgeline_interactions_all_cells_reps.pdf", width = 14, height = 8)
p
dev.off()
#################

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
#
##################
##proportional interactions per chrom
##################
#prop_dat <- r_dat2
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
##df with total number of common interactions per chrom
#commonInter <- prop_dat %>%
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
#
##proportional interactions per chrom pair
#pair_dat <- dat_long
##new col for chrom pair
#pair_dat$pair <- paste0(pair_dat$chrA,pair_dat$chrB)
##df with total number of interactions per chrom pair
#totInter <- pair_dat %>%
#  select(pair, ID) %>%
#  group_by(pair) %>%
#  dplyr::summarise(n = n())
##df with total number of common interactions per chrom pair
#commonInter <- gather(dat2, cell, zscore, 8:ncol(dat2), factor_key=TRUE)
##new col for chrom pair
#commonInter$pair <- paste0(commonInter$chrA,commonInter$chrB)
#commonInter <- commonInter %>%
#  select(pair, ID) %>%
#  group_by(pair) %>%
#  dplyr::summarise(n = n())
##df combining above 2 dfs
#prop_pairAll <- totInter
#prop_pairAll$commonInter <- prop_pairAll$n
#colnames(prop_pairAll) <- c("chrom", "totInter","commonInter")
#prop_pairAll$commonInter <- as.numeric(commonInter$n[match(prop_pairAll$chrom, commonInter$pair)])
#prop_pairAll$commonInter[is.na(prop_pairAll$commonInter)] <- 0
#prop_pairAll$percent <- (prop_pairAll$commonInter/prop_pairAll$totInter)*100 
#prop_pairAll$chrom <- as.factor(prop_pairAll$chrom)
#print("# chromosome pair(s) with highest (proportional) number of interactions")
#prop_pairAll[which(prop_pairAll$percent == max(prop_pairAll$percent)),]
##################

#################
print("#heatmap of common interactions z-scores by cell type and chromosome")
#################
r_dat2$AllChr <- gsub("chr","",r_dat2$AllChr)
r_dat2$AllChr <- as.factor(r_dat2$AllChr)
head(r_dat2)
#calculate mean zscore per chrom per cell type so heat map is accutate
hm_dat = r_dat2 %>% group_by(AllChr,cell) %>% dplyr::summarize(mzscore=mean(zscore, na.rm = TRUE))
head(hm_dat)
tail(hm_dat)
hm <- (ggplot(hm_dat, aes(AllChr, cell, fill = zscore))
       + geom_tile(aes(fill = mzscore), colour = "white")
#hm <- (ggplot(r_dat2, aes(AllChr, cell))
#       + geom_tile(aes(fill = zscore), colour = "white")
       + scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
       + labs(x = "Chromosome",
              y = "Cell",
              title = "Trans-chromosomal Interactions z-scores")
#       + facet_wrap(.~germL)
#       + theme(axis.text.x = element_text(angle = 90))
)
pdf("zscore_heatmap_interactions_chroms_all_cells_reps.pdf", width = 14, height = 8)
#png("zscore_heatmap_interactions_chroms_all_cells_reps.png")
hm
dev.off()
#################


#################
# violin and box plot
#################
head(dat2)
v <- (ggplot(dat2, aes(y=cell, x=zscore, fill = cell_noreps))
      + geom_violin(scale = "count")#height of violin proportional to data
      + geom_boxplot(alpha = 0.4, outlier.shape = NA)
      + scale_fill_manual(values = inauguration("inauguration_2021"))
      + labs(x = "z-score",
         y = "Cell",
         title = "Trans-chromosomal Interactions z-scores",
         fill = "")
)
pdf("zscore_violin_interactions_reps.pdf", width = 14, height = 8)
v
dev.off()

#################
print("# Pearson Correlation ")
#################
#pcDat <- dat2 %>% filter(cell_noreps != "Astrocytes_Cerebellum")
pcDat <- dat2
cell_nr <- as.character(unique(pcDat$cell_noreps))
cell_nr
cor_df <- data.frame(cell=character(),
                     correlation=numeric(),
                     cor_pvalue=numeric(),
                     value=character(),
                     stringsAsFactors=FALSE)
#i=as.character(cell_nr[6])
#i
for(i in cell_nr){
    tmpD <- pcDat %>% filter(cell_noreps == i)
    print(i)
    #zscore correlation
    tmp_zW <- tmpD %>% select(ID,cell, zscore) %>% spread(key=cell, value = zscore)
    result = cor.test(tmp_zW[,2],tmp_zW[,3] , method = "pearson")
    cor_val <- round(result$estimate,digits = 2)
    cor_pval <- round(result$p.value,digits = 2)
    cor_df <- rbind(cor_df, c(i,cor_val,cor_pval,"zscore"))
    #pvalue correlation
    tmp_zW <- tmpD %>% select(ID,cell, pvalue) %>% spread(key=cell, value = pvalue)
    result = cor.test(tmp_zW[,2],tmp_zW[,3] , method = "pearson")
    print(result)
    cor_val <- round(result$estimate,digits = 2)
    cor_pval <- round(result$p.value,digits = 2)
    cor_df <- rbind(cor_df, c(i,cor_val,cor_pval,"pvalue"))
}    
colnames(cor_df) <- c("cell","cor_val","cor_pval","value")
print("#correlation test results")
cor_df


###############
print("# checking correlation when adding constant to all zscores to make them all positive")
#adding 50 to all zscores so they are all positive
pcDat <- dat2
head(pcDat)
pcDat$zscore <- pcDat$zscore + 50
head(pcDat)
#################
print("# Pearson Correlation ")
#################
cell_nr <- as.character(unique(pcDat$cell_noreps))
cell_nr
cor_df <- data.frame(cell=character(),
                     correlation=numeric(),
                     cor_pvalue=numeric(),
                     value=character(),
                     stringsAsFactors=FALSE)
#i=as.character(cell_nr[6])
#i
for(i in cell_nr){
    tmpD <- pcDat %>% filter(cell_noreps == i)
    print(i)
    #zscore correlation
    tmp_zW <- tmpD %>% select(ID,cell, zscore) %>% spread(key=cell, value = zscore)
    result = cor.test(tmp_zW[,2],tmp_zW[,3] , method = "pearson")
    cor_val <- round(result$estimate,digits = 2)
    cor_pval <- round(result$p.value,digits = 2)
    cor_df <- rbind(cor_df, c(i,cor_val,cor_pval,"zscore"))
    #pvalue correlation
    tmp_zW <- tmpD %>% select(ID,cell, pvalue) %>% spread(key=cell, value = pvalue)
    result = cor.test(tmp_zW[,2],tmp_zW[,3] , method = "pearson")
    print(result)
    cor_val <- round(result$estimate,digits = 2)
    cor_pval <- round(result$p.value,digits = 2)
    cor_df <- rbind(cor_df, c(i,cor_val,cor_pval,"pvalue"))
}
colnames(cor_df) <- c("cell","cor_val","cor_pval","value")
print("#correlation test results")
cor_df

####################
#cn <- gsub("_", " ", i)
#corTit <- paste(cn,"Replicates\n correlation coefficient =",cor_val)
#corTit
#p <- (ggplot(tmp_zW, aes(x=zscore,y=zscoreDup))
#      #      +  geom_dotplot(binwidth = 0.5,stackdir = "centerwhole", dotsize = 12, alpha = 0.2)
#      #+ geom_point(alpha=0.4,size=3)
#      + geom_smooth(method="lm")
#      + labs(x = "z-score B",
#             y = "z-score A",
#             title = corTit)
#      # change ylim
#      #      + ylim(-0.25, 0.25)
#      # reference line at 0
#      #     + geom_hline(yintercept = 0, color="red")
#      # remove x-axis ticks and text
#      #      + theme(axis.ticks.x = element_blank(),
#      #             axis.text.x = element_blank())
#)
#p

#[grepl("Cardiomyocites_primitive",pcDat$cell)

