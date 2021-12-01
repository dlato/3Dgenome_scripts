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
#C=10
uchroms <- unique(tmpD$chr)
for(i in uchroms) {
#  i="chr10"
  print(i)
  chr_df <- tmpD %>% filter(chr == i)
  meandf <- chr_df %>% group_by(chr,start) %>%  dplyr::summarize(mNumCells=mean(totCells, na.rm = TRUE))
  C = gsub("chr","",i)
p <- (ggplot(data=meandf, aes(x = start, y = mNumCells)) 
      #+ geom_boxplot()
      + geom_point(size =3)
      +geom_line(size=1)
#      + geom_smooth(method = "loess", formula = y~x, colour="black")
      + labs(x=paste("Chromosome",C,"position [Mbp]"),
             y="Mean Number of Cells With Common Interaction",
             title = "Trans-chromosomal Common Interactions",
             fill = "")
      + scale_x_continuous(expand = c(0, 0))
#      + scale_y_continuous(expand = c(0, 0))
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
warnings()


print('##############################')
print("# COMMON INTERACTIONS USING SLIDING WINDOW")
print('##############################')
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
head(dat2)
window_sz <- 5000000
dat2$st1 <- as.numeric(dat2$st1)
dat2$st2 <- as.numeric(dat2$st2)
dat2$end1 <- as.numeric(dat2$end1)
dat2$end2 <- as.numeric(dat2$end2)
#initiate new df for common inters within window
common_window_df <- dat2[1,]
common_window_df <- common_window_df[-1,]
for(i in 1:length(dat2$ID)) {
#  i=1
  r=dat2[i,]
  rstA <- r$st1
  rendA <- r$end1
  rstB <- r$st1
  rendB <- r$end1
  #if start is less than window size, alter range start to be beginning of chrom
  if(r$st1 <= window_sz){
    rstA <- 0
  } else {
    rstA <- rstA - window_sz
  }
  if(r$st2 <= window_sz){
    rstB <- 0
  } else {
    rstB <- rstB - window_sz
  }
  #get end of chrom for each chrom
Aend <- chrInf$size[match(r$chrA, chrInf$chrom)]
Bend <- chrInf$size[match(r$chrB, chrInf$chrom)]
  #if end is within window size, alter range end to be end of chrom
  if(r$end1 >= Aend - window_sz){
    rendA <- Aend
  } else {
    rendA <- rendA + window_sz
  }
  if(r$end2 >= Bend - window_sz){
    rendB <- Bend
  } else {
    rendB <- rendB + window_sz
  }
  #grab the rows that fall within the range on BOTH chroms
  tmp_df <- dat2 %>% filter(chrA == r$chrA & chrB == r$chrB) %>%
                     #chrA needs to have both start and end within range 9bc we want the whole bin to be in the range)
                     filter(st1 >= rstA & st1 <= rendA & end1 >= rstA & end1 <= rendA) %>%
                     #chrB needs to have both start and end within range
                     filter(st2 >= rstA & st2 <= rendA & end2 >= rstA & end2 <= rendA)
  csum <- colSums(is.na(tmp_df))
  if (length(class(which(csum == nrow(tmp_df)))) == 0) {
   # print("#common interaction (one sig zscore in all cells)")
    common_window_df <- rbind(common_window_df, r)
  }
}#for

print("# DONE")
