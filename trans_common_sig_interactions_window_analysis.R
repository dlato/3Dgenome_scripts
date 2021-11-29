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
#            bin size (bp)
########################################

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
germlayer_file <- args[2]
allinters_file <- args[3]
tissue_file <- args[4]
bin_size <- args[5]
#Atype <- args[4]

##########
library(tidyr)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(ggforce)#for ridgeline
library(ggridges)#for ridgeline
library(ggbiplot)#for PCA
library(devtools)#for PCA
.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2")
library(circlize,lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2/") # for circos
library(regioneR,lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2/")#for permutation
library(factoextra, lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3/")#for PCA
library(harrypotter, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3") #for colours
library(regioneR)
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
##interaction data
Atype <- "1_vs_All"
tissue_file <- "tissue_system_info.txt"
dat_file <- "test_pairwise_dat.txt"
allinters_file <- "all_trans_interactions_1Mb.txt"
germlayer_file <- "germlayer_info.txt"
bin_size <- 1000000
library(factoextra)#for PCA
library(harrypotter) #for colours

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
print("all chrom pairs involved in ALL SIG interactions")
#split ID col
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
tmpD <- dat
tmpD$ID <- sub("B", "\\.B", as.character(tmpD$ID))
tmpD <- tmpD %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
#remove A and B from chrom names
tmpD$chrA <- gsub("A", "", tmpD$chrA)
tmpD$chrB <- gsub("B", "", tmpD$chrB)
tmpD$chrPair <- paste(tmpD$chrA,tmpD$chrB)
unique(tmpD$chrPair)

#remove NAs and therefore only look at common inters
df <- dat
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
print("all chrom pairs involved in common interactions")
tmpD <- dat2
tmpD$chrPair <- paste(dat2$chrA,dat2$chrB)
unique(tmpD$chrPair)

print("#########################")
print("# calculating number of cells that have the interaction as common")
print("#########################")
tot_cells_per_inter <- dat2
rownames(tot_cells_per_inter) <- tot_cells_per_inter$ID
tot_cells_per_inter <- tot_cells_per_inter %>% select(-ID,-chrA,-st1,-end1,-chrB,-st2,-end2)
tot_cells_per_inter <- as.data.frame(t(as.data.frame(t(apply(tot_cells_per_inter, MARGIN = 1, FUN = function(x) length(x[!is.na(x)]) )))))
tot_cells_per_inter$ID <- row.names(tot_cells_per_inter)
row.names(tot_cells_per_inter) <- NULL
colnames(tot_cells_per_inter) <- c("totCells","ID")
head(tot_cells_per_inter)

print("# histogram of how many cells most interactions have in common")
p <- (ggplot(data=tot_cells_per_inter, aes(totCells)) 
  + geom_histogram()
  + labs(y="Number of Interactions",
         x="Number of Cells With Common Interaction",
         title = "Trans-chromosomal Common Interactions",
         fill = "")
  + scale_x_continuous(expand = c(0, 0))
  + scale_y_continuous(expand = c(0, 0))
)
pdf("interactions_num_common_cells_histogram.pdf", width = 14, height = 8)
p
dev.off()

print("# common interactions along chrom length")
tmpD <- dat2
#new col for the chr start and end
tmpD$chrAID <- paste0(tmpD$chrA,".",tmpD$st1,".",tmpD$end1)
tmpD$chrBID <- paste0(tmpD$chrB,".",tmpD$st2,".",tmpD$end2)
tmpD <- tmpD %>% select(8:ncol(tmpD)) %>% select(chrAID, chrBID, everything())
tmpA <- tmpD %>% select(-chrBID)
colnames(tmpA)[1] <- "ID" 
tmpB <- tmpD %>% select(-chrAID)
colnames(tmpB)[1] <- "ID" 
tmpD <- rbind(tmpA,tmpB)
head(tmpD)
#calculate number of cells with interaction for each row
tot_cells_per_inter <- tmpD
tot_cells_per_inter <- as.data.frame(t(as.data.frame(t(apply(tot_cells_per_inter, MARGIN = 1, FUN = function(x) length(x[!is.na(x)]) )))))
tot_cells_per_inter <- bind_cols(tmpD$ID,tot_cells_per_inter$V1)
colnames(tot_cells_per_inter) <- c("ID","totCells")
head(tot_cells_per_inter)
#split ID col
colN <- c("chr","start","end")
tmpD <- tot_cells_per_inter %>% separate(ID, sep = "\\.", into = colN, remove = FALSE)
tmpD$start <- as.numeric(tmpD$start) / 1000000
head(tmpD)
summary(tmpD)
C=10
uchroms <- unique(tmpD$chr)
for(i in uchroms) {
  print(i)
  chr_df <- tmpD %>% filter(chr == i)
  C = gsub("chr","",i)
p <- (ggplot(data=chr_df, aes(x = start, y = totCells)) 
      #+ geom_boxplot()
      + geom_smooth(method = "loess", formula = y~x, colour="black")
      + labs(x=paste("Chromosome",C,"position [Mbp]"),
             y="Number of Cells With Common Interaction",
             title = "Trans-chromosomal Common Interactions",
             fill = "")
      + scale_x_continuous(expand = c(0, 0))
      + scale_y_continuous(expand = c(0, 0))
)
filename <- paste0("interactions_num_common_cells_chr",C,"_line.pdf")
pdf(filename, width = 14, height = 8)
print(p)
dev.off()
p <- (ggplot(data=chr_df, aes(x = start, y = totCells, group = start)) 
      + geom_boxplot()
      + labs(x=paste("Chromosome",C,"position [Mbp]"),
             y="Number of Cells With Common Interaction",
             title = "Trans-chromosomal Common Interactions",
             fill = "")
      + scale_x_continuous(expand = c(0, 0))
      + scale_y_continuous(expand = c(0, 0))
)
filename <- paste0("interactions_num_common_cells_chr",C,"_boxplot.pdf")
pdf(filename, width = 14, height = 8)
print(p)
dev.off()
}#for
######### EDIT HERE #################


print("# histogram of how many cells most interactions have in common")
p <- (ggplot(data=tot_cells_per_inter, aes(totCells)) 
      + geom_histogram()
      + labs(y="Number of Interactions",
             x="Number of Cells With Common Interaction",
             title = "Trans-chromosomal Common Interactions",
             fill = "")
      + scale_x_continuous(expand = c(0, 0))
      + scale_y_continuous(expand = c(0, 0))
)
pdf("interactions_num_common_cells_histogram.pdf", width = 14, height = 8)
p
dev.off()

#split ID col
tmpD <- tmpD %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
#remove A and B from chrom names
tmpD$chrA <- gsub("A", "", tmpD$chrA)
tmpD$chrB <- gsub("B", "", tmpD$chrB)
#count  each inter twice
tp_A <- tmpD %>% select(chrA,st1,end1,totCells) 
colnames(tp_A) <- c("chr","st","end","totCells")
tp_B <- tmpD %>% select(chrB,st2,end2,totCells) 
colnames(tp_B) <- c("chr","st","end","totCells")
tp_dat <- rbind(tp_A,tp_B)




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
                title = "Trans-chromosomal interactions (significant) z-scores")
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
         + theme(axis.text=element_text(size=18),panel.background = element_rect(fill = "grey85", colour = NA))
  )
  filename <- paste0("zscore_chrom",i,"_mean_tickplot_sig_Interactions.pdf")
  pdf(filename, width = 14, height = 8)
  print(tp)
  dev.off()
}#for