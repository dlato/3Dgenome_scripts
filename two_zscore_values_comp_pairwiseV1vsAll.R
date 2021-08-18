# test zscore difference between 1vsAll and pairwise 3Dflow calculation

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
oneVall <- args[1]
pairwise <- args[2]
#install.packages("dplyr", repos="http://cran.r-project.org", lib="/hpf/largeprojects/pmaass/3D-flow/scripts/src/R-lib")
#install.packages("FSA", repos="http://cran.r-project.org", lib="/hpf/largeprojects/pmaass/3D-flow/scripts/src/R-lib")


#library(msir)
#library(mclust)
library(dplyr)
library(ggplot2)
#library(FSA)

#library(msir)
#library(mclust)
#library(dplyr)
#library(ggplot2)
#library(FSA)



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
                  legend.title = element_blank(),
                  legend.text = element_text(size = 18),
                  #change the colour of facet label background
                  strip.background = element_rect(fill = "#E6E1EA"),
                  #remove space between facest
                  panel.spacing.x=unit(0, "lines"),
                  #                  legend.key = element_blank(),
                  legend.background=element_blank(),
                  #legend background
                  legend.key = element_rect(fill = NA),
                  #                  legend.position="none")
                  legend.position="top")
)
#########################################################################

#read in data
#oneVall <- "test_1vsAll_dat.txt"
#pairwise <- "test_pairwise_dat.txt"
#one vs All
o_df<-read.table(oneVall, header = TRUE)
head(o_df)
summary(o_df)
#subset df
#o_hESC <- o_df %>% select(ID,Thymus) 
o_hESC <- o_df %>% select(hESC_Dekker) 
colnames(o_hESC) <- c("ID", "hESC_o")
#pairwise
p_df<-read.table(pairwise, header = TRUE)
head(p_df)
summary(p_df)
#subset df
#p_hESC <- p_df %>% select(ID,Thymus) 
p_hESC <- p_df %>% select(hESC_Dekker)
colnames(p_hESC) <- c("ID", "hESC_p")
head(p_hESC)

#merge dfs by ID
m_df <- merge(o_hESC, p_hESC, by = 'ID')
head(m_df)
#select chrs 12 and 17 interactions
cistr_sox_df <- m_df[grep("chr12", m_df$ID), ]
#cistr_sox_df <- cistr_sox_df[grep("chr10", cistr_sox_df$ID), ]
cistr_sox_df <- cistr_sox_df[grep("chr17", cistr_sox_df$ID), ]
head(cistr_sox_df)

#plot matching zscores in pairwise and 1vsAll
p<-(ggplot(data=cistr_sox_df, aes(x=hESC_o, y=hESC_p)) 
    + geom_point(size=3, alpha=.3)
    + labs(x = "1vsAll z-score",
           y = "pairwise z-score",
           title = "hESC Dekker Interactions between chr 12 and chr 17")
    + geom_hline(yintercept = 0, color = "red", linetype = "dashed")
    + geom_vline(xintercept = 0, color = "red", linetype = "dashed")
)
pdf("1vsAll_pairwise_zscore_comparison_hESCDekker_chr12chr17.pdf", width = 14, height = 8)
p
dev.off()




###pairwise
##colnames(interactions)=c("index","id","chr1","st1","end1","chr2","st2","end2","dist","freq","mean","sd","zscore","pvalue")
##interactions <- interactions[,-1]
##1vsAll
#colnames(interactions)=c("id","chr1","st1","end1","chr2","st2","end2","dist","freq","mean","sd","zscore","pvalue","anchor_chrom")
#interactions <- interactions[,-14]
#interactions <- interactions[,c("id","chr1","st1","end1","chr2","st2","end2","dist","freq","mean","sd","pvalue","zscore")]
##interactions <- rbind(interactions,interactions[3,])
##interactions[3,]
##interactions[3,13] = 16
#id2 <- rownames(interactions)
#interactions$uID <- id2
#summary(interactions)
#interactions$interID <- do.call(paste0, interactions[c("chr1","st1","end1", "chr2", "st2", "end2")])
##interactions[which(interactions$interID == "chr1001000000chr121000000011000000"),]
##interactions[which(interactions$interID == "chr2189000000190000000chr1001000000"),]
#print("summary of ALL zscores (duplicated)")
#summary(interactions)
#print("summary of ALL zscores (duplicated)")
#outRange <- interactions %>% filter(zscore >= 3 | zscore <=-3)
#summary(outRange)
#print("percent of z-scores outside of -3 - +3 range")
#(length(outRange$zscore) / length(interactions$zscore)) * 100
##order df by interID
#interactions <- arrange(interactions, group_by = interID)
#interactions$interID <- as.factor(interactions$interID)
#interactions_orig <- interactions
##add duplicated rows as new col
##remove duplicates from orig df
#interactions <- interactions_orig[order(interactions_orig$interID, interactions_orig$zscore, decreasing=FALSE),]
#interactions <- interactions[!duplicated(interactions$interID),]
#interactions <- interactions[order(interactions$interID, interactions$zscore),]
##interactions[which(interactions$interID == "chr1001000000chr121000000011000000"),]
##keep highest zscore in separate df
## Reverse sort
#z <- interactions_orig[order(interactions_orig$interID, interactions_orig$zscore, decreasing=TRUE),]
## Keep only the first row for each duplicate 
#z <- z[!duplicated(z$interID),]
#z <- z[order(z$interID, z$zscore),]
##z[which(z$interID == "chr1001000000chr121000000011000000"),]
##z[which(z$interID == "chr2189000000190000000chr1001000000"),]
##rename zscore col
#colnames(z)[colnames(z) == 'zscore'] <- 'zscoreDup'
#z$zscoreDup <- as.numeric(z$zscoreDup)
##add duplicate as second col
#m1 <- left_join(interactions, select(z,  -id, -chr1, -st1, -end1, -chr2, -st2, -end2, -dist, -freq, -mean, -sd, -pvalue, -uID))
#head(m1)
#summary(m1)
##m1[which(m1$interID == "chr1001000000chr121000000011000000"),]
##m1[which(m1$interID == "chr2189000000190000000chr1001000000"),]
#
##calculate difference between two zscore values per interaction
##diff between two columns
#df <- m1%>%
#  group_by(interID)%>%
#  mutate(zDiff=abs(zscore-zscoreDup))
#head(df)
#summary(df)
##df[which(df$interID == "chr1001000000chr121000000011000000"),]
##df[which(df$interID == "chr2189000000190000000chr1001000000"),]
#high_diff <- df %>%
#          filter(zDiff >= 1.5)
#summary(high_diff)
#print("percent of zscores with zDiff >= 1.5")
#(length(high_diff$zscore) / length(df$zscore)) * 100
#
##where are these high zDiffs located?
#high_diff$chrID <- do.call(paste0, high_diff[c("chr1","chr2")])
##head(as.data.frame(high_diff))
##plot
##tit <- "Cardiomyocytes_Primitive_Rep_1_1Mb"
#mytitle <- gsub("_"," ",tit)
#p<-(ggplot(data=high_diff, aes(x=chrID)) 
#    + geom_histogram(stat = "count")
#    + labs(x = "",
#           y = "count",
#           title = paste(mytitle,"\n |Difference in z-scores| >= 1.5"))
#    + theme(axis.text.x = element_text(angle = 90))
#)
##p
#pdf("two_zscore_value_comparison_high_zdiff_histogram_per_chrom_pair.pdf", width = 14, height = 8)
#p
#dev.off()
#
##high_df_chroms <- as.data.frame(c(high_diff$chr1, high_diff$chr2))
##colnames(high_df_chroms) <- c("chroms")
##high_df_chroms$chroms <- as.factor(high_df_chroms$chroms)
##head(high_df_chroms)
##
##p<-(ggplot(data=high_df_chroms, aes(x=chroms)) 
##    + geom_histogram(stat = "count")
##    + labs(x = "chromosome",
##           y = "count",
##           title = paste(mytitle,"\n |Difference in z-scores| >= 1.5"))
##    + theme(axis.text.x = element_text(angle = 45))
##)
##calculating percentage of interactions involving each chrom with zDiff > 1.5
#uniq_chroms <- unique(c(df$chr1, df$chr2))
#uniq_chroms
##create empty df
#df2 <- data.frame(chroms=character(),
#                  perc=integer(),
#                  stringsAsFactors=FALSE)
#for (i in 1:length(uniq_chroms)){
#  sub_d <- as.data.frame(df[grep(uniq_chroms[i], df$id),])
#  high_zD_df <- as.data.frame(sub_d %>% filter(zDiff >= 1.5))
#  perc_val <- (length(high_zD_df$zDiff) / length(sub_d$id)) *100
#  df2[nrow(df2) + 1, ] <- c(uniq_chroms[i],perc_val)
#}
#colnames(df2) <- c("chroms","perc")
#p<-(ggplot(data=df2, aes(x=chroms, y=perc)) 
#    + geom_point(size =3)
#    + labs(x = "chromosome",
#           y = "Percent [%]",
#           title = paste(mytitle,"\n |Difference in z-scores| >= 1.5"))
#    + theme(axis.text.x = element_text(angle = 45))
#)
##p
#pdf("two_zscore_value_comparison_high_zdiff_histogram_per_chrom.pdf", width = 14, height = 8)
#p
#dev.off()
#
#
#######################
## Exploring chrY
#######################
#chrY_df <- high_diff[grep("chrY", high_diff$chrID),]
#head(as.data.frame(chrY_df))
##create new col to get chrY starts if chrY is anchor or target
#chrY_df <- chrY_df %>%
#           mutate(genomePos = ifelse(chr1 == "chrY", st1, st2))
##scale genomePos to be on the Mb scale
#chrY_df$genomePos <- chrY_df$genomePos / 1000000
###################save file##########################
#write.table(chrY_df,"chrY.high.zDiff.zscores.trans.txt", sep="\t", quote=FALSE, row.name=FALSE,col.name=FALSE)
#p<-(ggplot(data=chrY_df, aes(x=genomePos, y=zDiff)) 
#    + geom_point(size=3, alpha=.5)
#    + labs(x = "Chromosome Y Position [Mb]",
#           y = "|Difference in z-scores| \n >= 1.5",
#           title = mytitle)
#    + scale_x_continuous(limits = c(0, 58))
##    + theme(axis.text.x = element_text(angle = 45))
#)
##p
#pdf("two_zscore_value_comparison_high_zdiff_chrY_high_zDiff.pdf", width = 14, height = 8)
#p
#dev.off()
#
##all chrY interactions
#chrY_all <- df
#chrY_all$chrID <- do.call(paste0, chrY_all[c("chr1","chr2")])
#chrY_all <- chrY_all[grep("chrY", chrY_all$chrID),]
#chrY_all <- chrY_all %>%
#           mutate(genomePos = ifelse(chr1 == "chrY", st1, st2))
##new column for high zDiff
#chrY_all <- chrY_all %>%
#           mutate(zDiff_col = ifelse(zDiff >= 1.5, "High", "Low"))
#chrY_all$zDiff_col <- as.factor(chrY_all$zDiff_col)
#head(as.data.frame(chrY_all))
#summary(as.data.frame(chrY_all))
#levels(chrY_all$zDiff_col)
##scale genomePos to be on the Mb scale
#chrY_all$genomePos <- chrY_all$genomePos / 1000000
#p<-(ggplot(data=chrY_all, aes(x=genomePos, y=zDiff, color = zDiff_col)) 
#    + geom_point(size=3, alpha=.5)
#    + geom_hline(yintercept =1.5, colour = "red", linetype = 2)
#    + labs(x = "Chromosome Y Position [Mb]",
#           y = "|Difference in z-scores| \n >= 1.5",
#           title = mytitle)
#    + scale_color_manual(values = c("High" = "red", "Low" = "grey"))
#    + scale_x_continuous(limits = c(0, 58))
#        + theme(legend.position = "none")
#)
##p
#pdf("two_zscore_value_comparison_high_zdiff_chrY_all.pdf", width = 14, height = 8)
#p
#dev.off()
#######################
## Exploring chrYchrX
#######################
#chrYchrX_all <- df
#chrYchrX_all$chrID <- do.call(paste0, chrYchrX_all[c("chr1","chr2")])
#chrYchrX_all <- chrYchrX_all %>% filter(chrID == "chrXchrY" | chrID == "chrYchrX")
#chrYchrX_all$chrY_pos <- chrYchrX_all$st2
#chrYchrX_all$chrX_pos <- chrYchrX_all$st1
##new column for high zDiff
#chrYchrX_all <- chrYchrX_all %>%
#  mutate(zDiff_col = ifelse(zDiff >= 1.5, "High", "Low"))
#chrYchrX_all$zDiff_col <- as.factor(chrYchrX_all$zDiff_col)
#head(as.data.frame(chrYchrX_all))
#levels(chrYchrX_all$zDiff_col)
#chrYchrX_all <- chrYchrX_all %>% filter(zDiff_col == "High")
##scale genomePos to be on the Mb scale
#chrYchrX_all$chrY_pos <- chrYchrX_all$chrY_pos / 1000000
#chrYchrX_all$chrX_pos <- chrYchrX_all$chrX_pos / 1000000
###################save file##########################
#write.table(chrYchrX_all,"chrYchrX.high.zDiff.zscores.trans.txt", sep="\t", quote=FALSE, row.name=FALSE,col.name=FALSE)
#p<-(ggplot(data=chrYchrX_all, aes(x=chrX_pos, y=chrY_pos, color = zDiff_col)) 
#    + geom_point(size=3, alpha=.5)
#    + labs(x = "Chromosome X Position [Mb]",
#           y = "Chromosome Y Position [Mb]",
#           title = mytitle)
#    + scale_color_manual(values = c("High" = "red", "Low" = "grey"))
#    + theme(legend.position = "none")
#)
##p
#pdf("two_zscore_value_comparison_high_zdiff_chrYchrX_all.pdf", width = 14, height = 8)
#p
#dev.off()
#
## make col with cell name
#df$cell <- rep("CardPrimRep1",length(df$zDiff))
##title
#p <- (ggplot(df, aes(x=cell,y=zDiff))
#      #      +  geom_dotplot(binwidth = 0.5,stackdir = "centerwhole", dotsize = 12, alpha = 0.2)
#      + geom_violin(alpha=0.8)
##      + geom_boxplot(width=0.1, color="black",outlier.shape = NA)
#      + geom_boxplot(width=0.1, color="black")
#      + labs(x = "",
#             y = "|Difference in pairwise z-scores|",
#             title = mytitle)
#      # change ylim
#      #      + ylim(-0.25, 0.25)
#      # reference line at 0
# #     + geom_hline(yintercept = 0, color="red")
#      # remove x-axis ticks and text
#      + theme(axis.ticks.x = element_blank(),
#             axis.text.x = element_blank())
#)
##p
#pdf("two_zscore_value_comparison.pdf", width = 14, height = 8)
#p
#dev.off()
#
#print("CORRELATION between zscore values")
##correlation graph
#result = cor.test(df$zscore,df$zscoreDup , method = "pearson")
#result
#cor_val <- round(result$estimate,digits = 2)
#df$cell <- rep("CardPrimRep1",length(df$zDiff))
#corTit <- paste(mytitle,"\n correlation coefficient =",cor_val)
#p <- (ggplot(df, aes(x=zscore,y=zscoreDup))
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
##p
#pdf("two_zscore_value_comparison_correlation.pdf", width = 14, height = 8)
#p
#dev.off()
#
##mean comparison
##add duplicated rows as new col
##remove duplicates from orig df
#interactions <- interactions_orig[order(interactions_orig$interID, interactions_orig$mean, decreasing=FALSE),]
#interactions <- interactions[!duplicated(interactions$interID),]
#interactions <- interactions[order(interactions$interID, interactions$mean),]
##interactions[which(interactions$interID == "chr1001000000chr121000000011000000"),]
##keep highest mean in separate df
## Reverse sort
#z <- interactions_orig[order(interactions_orig$interID, interactions_orig$mean, decreasing=TRUE),]
## Keep only the first row for each duplicate 
#z <- z[!duplicated(z$interID),]
#z <- z[order(z$interID, z$mean),]
##z[which(z$interID == "chr1001000000chr121000000011000000"),]
##z[which(z$interID == "chr2189000000190000000chr1001000000"),]
##rename mean col
#colnames(z)[colnames(z) == 'mean'] <- 'meanDup'
#z$meanDup <- as.numeric(z$meanDup)
##add duplicate as second col
#m1 <- left_join(interactions, select(z,  -id, -chr1, -st1, -end1, -chr2, -st2, -end2, -dist, -freq,-zscore,  -sd, -pvalue, -uID))
#head(m1)
#summary(m1)
##m1[which(m1$interID == "chr1001000000chr121000000011000000"),]
##m1[which(m1$interID == "chr2189000000190000000chr1001000000"),]
#
##calculate difference between two mean values per interaction
##diff between two columns
#df <- m1%>%
#  group_by(interID)%>%
#  mutate(mDiff=abs(mean-meanDup))
#head(df)
#summary(df)
##df[which(df$interID == "chr1001000000chr121000000011000000"),]
##df[which(df$interID == "chr2189000000190000000chr1001000000"),]
#
#
#
#
## make col with cell name
#df$cell <- rep("CardPrimRep1",length(df$mDiff))
#p <- (ggplot(df, aes(x=cell,y=mDiff))
#      #      +  geom_dotplot(binwidth = 0.5,stackdir = "centerwhole", dotsize = 12, alpha = 0.2)
#      + geom_violin(alpha=0.8)
##      + geom_boxplot(width=0.1, color="black",outlier.shape = NA)
#      + geom_boxplot(width=0.1, color="black")
#      + labs(x = "",
#             y = "|Difference in 1vsAll LOESS mean|",
#             title = mytitle)
#      # change ylim
#      #      + ylim(-0.25, 0.25)
#      # reference line at 0
#      #     + geom_hline(yintercept = 0, color="red")
#      # remove x-axis ticks and text
#            + theme(axis.ticks.x = element_blank(),
#                   axis.text.x = element_blank())
#)
##p
#pdf("two_zscore_mean_value_comparison.pdf", width = 14, height = 8)
#p
#dev.off()
#
##sd comparison
##add duplicated rows as new col
##remove duplicates from orig df
#interactions <- interactions_orig[order(interactions_orig$interID, interactions_orig$sd, decreasing=FALSE),]
#interactions <- interactions[!duplicated(interactions$interID),]
#interactions <- interactions[order(interactions$interID, interactions$sd),]
##interactions[which(interactions$interID == "chr1001000000chr121000000011000000"),]
##keep highest sd in separate df
## Reverse sort
#z <- interactions_orig[order(interactions_orig$interID, interactions_orig$sd, decreasing=TRUE),]
## Keep only the first row for each duplicate 
#z <- z[!duplicated(z$interID),]
#z <- z[order(z$interID, z$sd),]
##z[which(z$interID == "chr1001000000chr121000000011000000"),]
##z[which(z$interID == "chr2189000000190000000chr1001000000"),]
##rename sd col
#colnames(z)[colnames(z) == 'sd'] <- 'sdDup'
#z$sdDup <- as.numeric(z$sdDup)
##add duplicate as second col
#m1 <- left_join(interactions, select(z,  -id, -chr1, -st1, -end1, -chr2, -st2, -end2, -dist, -freq,-zscore,  -mean, -pvalue, -uID))
#head(m1)
#summary(m1)
##m1[which(m1$interID == "chr1001000000chr121000000011000000"),]
##m1[which(m1$interID == "chr2189000000190000000chr1001000000"),]
#
##calculate difference between two sd values per interaction
##diff between two columns
#df <- m1%>%
#  group_by(interID)%>%
#  mutate(sDiff=abs(sd-sdDup))
#head(df)
#summary(df)
##df[which(df$interID == "chr1001000000chr121000000011000000"),]
##df[which(df$interID == "chr2189000000190000000chr1001000000"),]
#
#
#
## make col with cell name
#df$cell <- rep("CardPrimRep1",length(df$sDiff))
#
#p <- (ggplot(df, aes(x=cell,y=sDiff))
#      #      +  geom_dotplot(binwidth = 0.5,stackdir = "centerwhole", dotsize = 12, alpha = 0.2)
#      + geom_violin(alpha=0.8)
##      + geom_boxplot(width=0.1, color="black",outlier.shape = NA)
#      + geom_boxplot(width=0.1, color="black")
#      + labs(x = "",
#             y = "|Difference in 1vsAll LOESS sd|",
#             title = mytitle)
#      # change ylim
#      #      + ylim(-0.25, 0.25)
#      # reference line at 0
#      #     + geom_hline(yintercept = 0, color="red")
#      # remove x-axis ticks and text
#            + theme(axis.ticks.x = element_blank(),
#                   axis.text.x = element_blank())
#)
##p
#pdf("two_zscore_sd_value_comparison.pdf", width = 14, height = 8)
#p
#dev.off()
#