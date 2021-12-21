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
#            window size (bp) *will calculate window size input in BOTH directions around bin. i.e. bin 5-10 would be 0-15 with a window size of 5
#            full path and name of output file
########################################

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
germlayer_file <- args[2]
allinters_file <- args[3]
tissue_file <- args[4]
bin_size <- args[5]
window_sz <- as.numeric(args[6])
outfile <- args[7]
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
###interaction data
#Atype <- "1_vs_All"
#tissue_file <- "tissue_system_info.txt"
#dat_file <- "test_pairwise_dat.txt"
#allinters_file <- "all_trans_interactions_1Mb.txt"
#germlayer_file <- "germlayer_info.txt"
#bin_size <- 1000000
#window_sz <- 5000000
#outfile <- "test_common_inters_df.txt"
#library(factoextra)#for PCA
#library(harrypotter) #for colours
#library(circlize) #for circos

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
print("window size")
window_sz
class(window_sz)
type(window_sz)
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
head(dat2)
dat2$st1 <- as.numeric(dat2$st1)
dat2$st2 <- as.numeric(dat2$st2)
dat2$end1 <- as.numeric(dat2$end1)
dat2$end2 <- as.numeric(dat2$end2)
chrInf$size <- as.numeric(chrInf$size)
print("#initiate new df for common inters within window")
common_window_df <- dat2[1,]
common_window_df <- common_window_df[-1,]
print("window_sz")
print(window_sz)
print(type(window_sz))
print(class(window_sz))
print("for loop")
for(i in 1:length(dat2$ID)) {
#  i=32
print ("--------------")
  r=dat2[i,]
  print(r)
  rstA <- r$st1
  rendA <- r$end1
  rstB <- r$st2
  rendB <- r$end2
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
print("Aend")
  print(Aend)
print("Bend")
  print(Bend)
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
  print(rstA)
  print(rendA)
  print(rstB)
  print(rendB)
  #grab the rows that fall within the range on BOTH chroms
  tmp_df <- dat2 %>% filter(chrA == r$chrA & chrB == r$chrB) %>%
                     #chrA needs to have both start and end within range 9bc we want the whole bin to be in the range)
                     filter(st1 >= rstA & st1 <= rendA & end1 >= rstA & end1 <= rendA) %>%
                     #chrB needs to have both start and end within range
                     filter(st2 >= rstB & st2 <= rendB & end2 >= rstB & end2 <= rendB)
#  tmp_df
  csum <- colSums(is.na(tmp_df))
  if (length(which(csum == nrow(tmp_df))) == 0) {
#    print("#common interaction (one sig zscore in all cells)")
    common_window_df <- rbind(common_window_df, r)
  }
}#for
print("#save df")
options(scipen=999)
type(outfile)
class(outfile)
outfile
write.table(common_window_df %>% select(chrA, st1, end1,chrB,st2,end2), file = as.character(outfile), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

print("###############################")
print("# graphs and tests for common inters with window")
print("###############################")
dat2 <- common_window_df
print("#circos plot of common interactions")
circos_df <- dat2 %>% select(chrA, st1, end1, chrB, st2,end2)
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
pdf("circos_common_inters.pdf", width = 14, height = 8)
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
#rcols <- scales::alpha(ifelse(sign(nuc1$st1-nuc1$end1)!=sign(nuc2$st2-nuc2$end2),"#f46d43","#66c2a5"),alpha=0.4)
#rcols <- scales::alpha(ifelse(sign(nuc1$st1-nuc1$end1)!=sign(nuc2$st2-nuc2$end2),"black","red"))
#circos.genomicLink(nuc1,nuc2,col=rcols,border=NA)
circos.genomicLink(nuc1,nuc2)
dev.off()

rev_chrs_len_ord <- rev(p_chr_ord)
chrs_len_ord <- p_chr_ord
chrInf
#wide to long format
r_dat <- gather(dat2, cell, zscore, 8:ncol(dat2), factor_key=TRUE)
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
r_dat3 <- r_dat
#counting each interaction twice (once for each chrom in interaction)
anchD <- r_dat
anchD$AllChr <- anchD$chrA
anchD$AllSt <- anchD$st1
tarD <- r_dat
tarD$AllChr <- tarD$chrB
tarD$AllSt <- tarD$st2
r_dat2 <- rbind(anchD,tarD)
pq_DAT <- r_dat2
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
p <- (ggplot(r_dat2, aes(x = AllSt, y = AllChr))
      + geom_density_ridges(scale = 2, alpha = 0.3, rel_min_height = 0.001) 
      #  + geom_segment(data = chrInf2, aes(x=centromere, xend=centromere, 
      + geom_segment(data = chrInf2, aes(x=centromere, xend=centromere, 
                                         y=as.numeric(AllChr), yend=as.numeric(AllChr) +0.9),
                     color = "red")
      #  + geom_segment(data = chrInf2, aes(x=size, xend=size, 
      + geom_segment(data = chrInf2, aes(x=size, xend=size, 
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
pdf("chr_ridgeline_common_interactions_all_cells_FAKE_DATA.pdf", width = 14, height = 8)
p
dev.off()
#without fake data
ndat2 <- r_dat2 %>% filter(cell != "fake_cell")
chsindf <- unique(ndat2$AllChr)
red_levels <- factor(rev_chrs_len_ord[rev_chrs_len_ord %in% chsindf])
ndat2$AllChr <- factor(ndat2$AllChr, levels =red_levels)
chrInf_hm <- chrInf2 %>% filter(AllChr %in% chsindf)
chrInf_hm$AllChr <- factor(chrInf_hm$AllChr, levels =red_levels)
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
print("# max genomic position for each chromosome (common interactions)")
print(as_tibble(ndat2 %>% group_by(AllChr) %>% dplyr::summarise(max_pos = max(AllSt))), n = 100)
print("# chromosome centromere and size info")
chrInf_hm
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
#``
print("# differnce in interactions on P and Q arms")
pq_dat <- pq_DAT %>% filter(cell != "fake_cell")
summary(pq_dat)
pq_tot_df <- data.frame(cell=character(),
                        AllChr=character(),
                        arm=character(),
                        n=integer(),
                        stringsAsFactors=FALSE)
for (i in unique(pq_dat$AllChr)){
  #create df with interactions classified in each arm for each chrom
  print(i)
  chr_df <- pq_dat %>% filter(AllChr == i)
  print(head(chr_df))
  cd <- chrInf %>% filter(chrom == i)
  c <- cd$centromere[1] 
  #find which arm is the short one and which is the long one
  lB <- cd$size[1] - cd$centromere[1]
  print(lB)
  l_df <- c(c,lB)
  s <- which.min(l_df)
  print(s)
  A <- ifelse(s == 1, "p","q")
  B <- ifelse(s == 1, "q","p")
  print(A)
  print(B)
  chr_df <- chr_df %>%  mutate(arm = if_else(AllSt <= c,A,B))
  print(head(chr_df))
  #add total number of p and q interactions on each arm to a new df
  appD <- chr_df %>% group_by(cell,AllChr, arm) %>% dplyr::summarise(n = n())
  print(appD)
  pq_tot_df <- rbind(appD,pq_tot_df)
  #  pq_perm_df <- rbind(pq_perm_df,chr_df)
}
head(pq_tot_df)
head(pq_dat)

#get df with if the arm before centromere is p or q
chr_pq <- data.frame(chrom=character(),
                     centromere=character(),
                     arm=character(),
                     stringsAsFactors=FALSE)
for (i in unique(chrInf$chrom)){
  #create df with interactions classified in each arm for each chrom
  cd <- chrInf %>% filter(chrom == i)
  c <- cd$centromere[1] 
  #find which arm is the short one and which is the long one
  lB <- cd$size[1] - cd$centromere[1]
  l_df <- c(c,lB)
  s <- which.min(l_df)
  A <- ifelse(s == 1, "p","q")
  B <- ifelse(s == 1, "q","p")
  tmprow <- c(i,c,A)
  #add total number of p and q interactions on each arm to a new df
  chr_pq <- rbind(chr_pq,tmprow)
}
colnames(chr_pq) <- c("chrom", "centromere", "arm")
head(chr_pq)
#pq dataframe for permutation test
pq_perm_df <- data.frame(matrix(ncol = 11, nrow = 0))
pq_dat2 <- unique(pq_dat %>% select(-AllChr, -AllSt))
for (r in 1:nrow(pq_dat2)){
  #get chromA of row
  rdf <- pq_dat2[r,]
  cA <- rdf$chrA[1]
  endA <- rdf$end1[1]
  centdf <- chr_pq %>% filter(chrom == cA)
  centA <- as.numeric(as.character(centdf$centromere[1]))
  aA <- centdf$arm[1]
  aB <- ifelse(aA == "p","q","p")
  #create df with interactions classified in each arm for each chrom
  rdf$armA[1] <- if_else(endA <= centA,aA,aB)
  #interaction chrom B
  cB <- rdf$chrB[1]
  endB <- rdf$end2[1]
  centdf <- chr_pq %>% filter(chrom == cB)
  centB <- as.numeric(as.character(centdf$centromere[1]))
  aA <- centdf$arm[1]
  aB <- ifelse(aA == "p","q","p")
  #create df with interactions classified in each arm for each chrom
  rdf$armB[1] <- if_else(endB <= centB,aA,aB)
  pq_perm_df <- rbind(pq_perm_df,rdf)
}
head(pq_perm_df)

#convert common inters to GRanges object
#get both chromosomes as diff rows for each interaction
firstchr <- unique(pq_dat %>% select(chrA, st1, end1, ID))
secondchr <- unique(pq_dat %>% select(chrB, st2, end2, ID))
firstchr$end1 <- as.numeric(as.character(firstchr$end1))
secondchr$end2 <- as.numeric(as.character(secondchr$end2))
colnames(firstchr) <- c("chrom","start","end", "ID")
colnames(secondchr) <- c("chrom","start","end", "ID")
commoninters_perm <- rbind(firstchr, secondchr)
head(commoninters_perm)
commonintersGRange <- toGRanges(commoninters_perm)
#df with just p arm ranges for each chrom
parm_perm <- chrInf %>% select(chrom, centromere)
parm_perm$start <- rep(0,nrow(parm_perm))
parm_perm <- parm_perm %>% select(chrom,start,centromere)
parmGRanges <- toGRanges(parm_perm)

#df with just p arm ranges for each chrom
qarm_perm <- chrInf %>% select(chrom, centromere,size)
qarmGRanges <- toGRanges(qarm_perm)

set.seed(369)
print('# permutation test to see if common interactions overlap more with p arm regions')
pt <- overlapPermTest(commonintersGRange, parmGRanges, ntimes=10000, genome="hg38", count.once=TRUE)
pt
pdf("commonInters_vs_pArm_permutation.pdf", width = 14, height = 8)
plot(pt)
dev.off()
print('# permutation test to see if common interactions overlap more with q arm regions')
pt <- overlapPermTest(commonintersGRange, qarmGRanges, ntimes=10000, genome="hg38", count.once=TRUE)
pt
pdf("commonInters_vs_qArm_permutation.pdf", width = 14, height = 8)
plot(pt)
dev.off()

print("### boxplot/violin for pq arms")
#get dataframe with total number of interactions per chrom arm
tot_inter_arm <- data.frame(chrom=character(),
                            inter=numeric(),
                            arm=character(),
                            stringsAsFactors=FALSE)
for (i in unique(chrInf$chrom)){
  #i="chr21"
  chr_all <- allinters %>% filter(chrA == i | chrB == i)
  A_chr_all <- toGRanges(chr_all %>% filter(chrA == i) %>% select(chrA, startA, endA))
  B_chr_all <- toGRanges(chr_all %>% filter(chrB == i) %>% select(chrB, startB, endB))
  chrInf_chr <- chrInf %>% filter(chrom == i)
  chrInf_chr$start <- 0
  p_range <- toGRanges(chrInf_chr %>% select(chrom, start, centromere))
  q_range <- toGRanges(chrInf_chr %>% select(chrom, centromere, size))
  pA <- findOverlaps(p_range,A_chr_all)
  pB <- findOverlaps(p_range,B_chr_all)
  qA <- findOverlaps(q_range,A_chr_all)
  qB <- findOverlaps(q_range,B_chr_all)
  q_tot <- sum(length(qA),length(qB))
  p_tot <- sum(length(pA),length(pB))
  tmprow_p <- c(i,p_tot,"p")
  tmprow_q <- c(i,q_tot,"q")
  #add total number of p and q interactions on each arm to a new df
  tot_inter_arm <- rbind(tot_inter_arm,tmprow_p,tmprow_q)
}#for
colnames(tot_inter_arm) <- c("chr", "tot_inters", "arm")
print("# total possible inters per chrom arm (head)")
head(tot_inter_arm)
A_dat <- pq_perm_df %>% select(chrA,st1,end1,armA)
colnames(A_dat) <- c("chr","start","end","arm")
B_dat <- pq_perm_df %>% select(chrB,st2,end2,armB)
colnames(B_dat) <- c("chr","start","end","arm")
pq_box_dat <- rbind(A_dat,B_dat)
pq_sum <- pq_box_dat %>%  group_by(chr, arm) %>% dplyr::summarize(common=n())
pq_sum
#bin_df <- chrInf %>% select(chrom,size)
#summary(bin_df)
#bin_df$numBins <- ceiling(as.numeric(as.character(bin_df$size)) / as.numeric(as.character(bin_size)))
#bin_df <- bin_df %>% select(-size)
#colnames(bin_df) <- c("chr","numBins")
box_df <-merge(pq_sum,tot_inter_arm, by = c("chr","arm"))
box_df$prop <- box_df$common / as.numeric(as.character(box_df$tot_inters))
box_df
set.seed(369)
p <- (ggplot(box_df, aes(x = arm, y = prop))
      #+ geom_violin(fill="grey90", scale = "count")#"count" makes width of violins proportional to number of values
      + geom_boxplot(fill="grey95",width = 0.3,outlier.shape = NA, alpha = 0.5)
      + geom_jitter(fill="grey90", alpha = 0.4, width = 0.15, size =4)
      #      + stat_density_ridges(quantile_lines = TRUE, alpha = 0.3, scale=2, quantiles = 2, rel_min_height = 0.001)
      #      #+ geom_density_ridges(scale = 4, alpha = 0.3) 
      + labs(y="Proportion of Common Trans-chromosomal Interactions",
             x="Chromosome Arm",
             title = "Proportion of Common Trans-chromosomal Interactions per Chromosome Arm")
      #      #  + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
      #      + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
      #      + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
)
pdf("PQ_arm_violin_common_interactions.pdf", width = 14, height = 8)
p
dev.off()

print("#calculate mean zscore per chrom per bin so that the heatmap is accurate")
head(r_dat2)
hm_dat = ndat2 %>% group_by(AllChr, AllSt) %>% dplyr::summarize(mzscore=mean(zscore))
head(hm_dat)
chrs_ord <- gsub("chr", "", chrs_len_ord)
hm_dat$AllChr <- factor(hm_dat$AllChr, levels=chrs_len_ord)
hm_dat$AllChr <- gsub("chr", "", hm_dat$AllChr)
hm <- (ggplot(hm_dat, aes(x=AllSt,ordered(AllChr, levels=rev(chrs_ord)), fill = mzscore))
       #hm <- (ggplot(hm_dat, aes(AllChr, cell, fill = zscore))
       #       + geom_tile(aes(fill = mzscore), colour = "white")
       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score", na.value = "grey")
       + geom_tile(aes(fill = mzscore), width = 1, height = 1)
       #       + scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
       + labs(x = "Genomic Postion [Mbp]",
              y = "Chromosome",
              title = "Common Trans-chromosomal Interactions z-scores")
       + theme(panel.background = element_rect(fill = "grey85", colour = NA))
       #       + facet_wrap(.~germL)
       #       + theme(axis.text=element_text(size=12))
)
pdf("zscore_mean_heatmap_common_interactions_chroms_all_cells.pdf", width = 14, height = 8)
hm
dev.off()



print("#ridgeline of all zscores (in all chroms) across cell types")
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

#adding germlayer info
r_dat3$germL <- r_dat3$cell
r_dat3$germL <- as.factor(gl_df$germLayer[match(r_dat3$cell, gl_df$cell)])
#re-order based on gl_ord
r_dat3$germL <- factor(r_dat3$germL, levels=gl_ord)
r_dat3 <- r_dat3[order(r_dat3$germL),]
gl_cell_ord <- unique(r_dat3$cell)
gl_cell_ord
r_dat3$cell <- factor(r_dat3$cell, levels=rev(gl_cell_ord))
#r_dat2 <- r_dat2 %>% mutate(cell = factor(cell,levels=cell))
levels(r_dat3$germL)
levels(r_dat3$cell)
head(r_dat3)
r_dat3 <- r_dat3 %>% filter(cell != "fake_cell")
#remove NAs because these mean we did not have this information in the sequencing data (since we are reading in the all zscores df)
r_dat3 <- r_dat3 %>% filter(!is.na(zscore))


p <- (ggplot(r_dat3, aes(x = zscore, y = cell, fill = germL))
      #p <- (ggplot(r_dat2, aes(x = zscore, y = cell, fill = germL))
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

##################
#print("#test between means of germ layers")
#print("#### common interactions ####")
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

print("# Tissue/system breakdown")
t_dat2 <- r_dat3
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
print("#common interactions z-scores broken down by chrom class")
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
#print("#test between means of chrom classes")
#print("#Test each group for normality")
#print("sig = reject normality null")
#chrClass_dat$chrClass <- as.factor(chrClass_dat$chrClass)
#chrClass_dat %>%
#  group_by(chrClass) %>%
#  summarise(W = shapiro.test(zscore)$statistic,
#            p.value = shapiro.test(zscore)$p.value)
#print("#Perform the Kruskal-Wallis test")
#print("sig = mean is diff btwn groups")
#kruskal.test(zscore ~ chrClass, data=chrClass_dat)
#print("# check which groups have sig diff")
#print("# perform pairwise wilcoxon test with FDR (Benjamini-Hochberg) correction")
#pairwise.wilcox.test(chrClass_dat$zscore, chrClass_dat$chrClass,
#                 p.adjust.method = "BH")

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
      + scale_fill_manual(values =c("#FAC9A1", "#013040"), labels= c("Total Possible Interactions", "Common Interactions"))
      + labs(y="Number of Interactions [100,000]",
             x="Chromosome",
             title = "Trans-chromosomal Interactions per Chromosome",
             fill = "")
)
pdf("proportion_per_chrom_common_interactions_all_cells.pdf", width = 14, height = 8)
p
dev.off()
print("#### chromosome with the highest proportional number of common interactions")
prop_datAll[which(prop_datAll$percent == max(prop_datAll$percent)),]
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
       + geom_tile(aes(fill = zscore), width = 1, height = 1)
       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score", na.value = "grey")
       #       + scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
       + labs(x = "Chromosome",
              y = "Cell",
              title = "Common Trans-chromosomal Interactions z-scores")
       #       + facet_wrap(.~germL)
       + theme(axis.text=element_text(size=12),
               panel.background = element_rect(fill = "grey85", colour = NA))
)
pdf("zscore_mean_heatmap_common_interactions_chroms_per_cell.pdf", width = 14, height = 8)
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
otherDir <- hm_dat2 %>% select(chrPair,chrB,chrA,cell,mzscore)
colnames(otherDir) <- c("chrPair", "chrA","chrB", "cell", "mzscore")
head(otherDir)
hm_bothdirs <- bind_rows(hm_dat2,otherDir)
hm_bothdirs$chrA <- factor(hm_bothdirs$chrA, levels=chrs_ord)
hm_bothdirs$chrB <- factor(hm_bothdirs$chrB, levels=chrs_ord)
hm <- (ggplot(hm_bothdirs, aes(chrA, ordered(chrB, levels=rev(chrs_ord)), fill = mzscore))
       + geom_tile(aes(fill = mzscore), width = 1, height = 1)
       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score per chromosomal pair", na.value = "grey")
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
               panel.background = element_rect(fill = "grey85", colour = NA),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank(),
               axis.text=element_text(size=5))
)
pdf("zscore_chrom_pair_mean_heatmap_common_interactions.pdf", width = 14, height = 8)
hm
dev.off()

print("# boxplot of chrom pair mean zscores (common inters)")
bpD <- hm_bothdirs %>% select(chrPair,cell,mzscore) %>% distinct()
bpD$chrPair <- gsub("\\.","", bpD$chrPair)
set.seed(369)
p <- (ggplot(bpD, aes(y =mzscore,x=reorder(chrPair,-mzscore, na.rm = TRUE)))
      + geom_boxplot()
      + geom_jitter(width = 0.3, alpha = 0.4)
      #+ scale_colour_manual(values =c("#FAC9A1", "#013040"), labels= c("All", "Significant"))
      + labs(y="Mean z-score per Chromosome Pair",
             x="Chromosome Pair",
             title = "Common Trans-chromosomal Interactions Across Cell Types",
             colour = "Interactions")
      + theme(axis.text.x = element_text(angle = 90))
)
pdf("zscore_chromPair_common_inters_all_cells_line.pdf", width = 14, height = 8)
p
dev.off()
#lineplot
set.seed(369)
p <- (ggplot(bpD, aes(y =mzscore,x=reorder(chrPair,-mzscore, na.rm = TRUE), group = 1))
      + geom_point(alpha = 0.4)
      + geom_smooth(method = "loess", formula = y~x, colour = 'red')
      #+ geom_boxplot()
      #+ geom_jitter(width = 0.3, alpha = 0.4)
      #+ scale_colour_manual(values =c("#FAC9A1", "#013040"), labels= c("All", "Significant"))
      + labs(y="Mean z-score per Chromosome Pair",
             x="Chromosome Pair",
             title = "Common Trans-chromosomal Interactions Across Cell Types",
             colour = "Interactions")
      + theme(axis.text.x = element_text(angle = 90))
)
pdf("zscore_chromPair_common_inters_all_cells_boxplot.pdf", width = 14, height = 8)
p
dev.off()
print("# top 10 chrom pair ordered highest to lowest MEAN of mean zscore per chrom pair (across all cells)")
meanprop <- bpD %>% group_by(chrPair) %>% dplyr::summarize(mProp=mean(mzscore, na.rm = TRUE))
meanprop[order(meanprop$mProp, decreasing = TRUE),]


print("# DONE")
