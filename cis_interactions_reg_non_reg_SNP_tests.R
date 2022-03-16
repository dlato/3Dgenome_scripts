########################################
#validating known trans-chromosomal interactions
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: reg SNPs filtered interactions file (tab separated, first 6 columns are a bed-like format, last few columns are other info)
#            non-reg SNPs filtered interactions file (tab separated, first 6 columns are a bed-like format, last few columns are other info)
#            reg SNPs file (tab separated)
#            non-reg SNPs file (tab separated)
#            bin size (bp)
#            output file prefix (character)
#            cell order for graphs (file with name of cells as it matches the data frame, one cell per line)
########################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
regSNPs_intersfile <- args[1]
nonRegSNPs_intersfile <- args[2]
regSNPs_file <- args[3]
nonRegSNPs_file <- args[4]
bin_size <- args[5]
outprefix <- args[6]
cellsfile <- args[7]

##########
library(tidyr)
library(forcats)
library(dplyr)
library(ggplot2)
#library(nortest, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.3.6.1") #for normality test with large sample size
#library(hexbin)
library(ggforce)#for ridgeline
library(ggridges)#for ridgeline
.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.3.6.1")
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
                 # legend.title = element_blank(),
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


print("#read in files")
##interaction data
#regSNPs_intersfile <- "VSMC_cis50Kb_SNPs_MERGED_overlapping_intearctions_for_circos.txt"
#nonRegSNPs_intersfile <- "nonRegSNPs_all_cells_cis50Kb_SNPs_MERGED_overlapping_intearctions_for_circos.txt"
#regSNPs_file <- "VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt"
#nonRegSNPs_file <- "sentinel_snps_not_reg_for_HiC_all.txt"
#bin_size <- 50000
#outprefix <- "test_cis_SNP_blood_pressure"
##cellsfile <- "cell_subset_CM.txt"
#cellsfile <- "cell_subset2.txt"
#library(harrypotter)
#library(factoextra)
#library(hexbin)
#library(ggVennDiagram)
#library(ggupset)

cells_sub <- read.table(cellsfile)
cells_sub

#read inters files
regSNPs_inters <- read.table(regSNPs_intersfile, header = F)
colnames(regSNPs_inters) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
#filter by desired cells
regSNPs_inters <- regSNPs_inters %>% filter(cell %in% cells_sub$V1)
head(regSNPs_inters)
nonRegSNPs_inters <- read.table(nonRegSNPs_intersfile, header = F)
colnames(nonRegSNPs_inters) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
#filter by desired cells
nonRegSNPs_inters <- nonRegSNPs_inters %>% filter(cell %in% cells_sub$V1)
head(nonRegSNPs_inters)
#merge reg and non-reg inters
SNP_inters <- full_join(regSNPs_inters,nonRegSNPs_inters, by = c("chrA","st1","end1","chrB","st2","end2","cell","zscore")) %>%
              distinct()
head(SNP_inters)

#num inters per bin
print("#counting each interaction twice (once for each chrom in interaction)")
anchD <- SNP_inters
anchD$AllChr <- anchD$chrA
anchD$AllSt <- anchD$st1
anchD$AllEnd <- anchD$end1
tarD <- SNP_inters
tarD$AllChr <- tarD$chrB
tarD$AllSt <- tarD$st2
tarD$AllEnd <- tarD$end2
SNP_inters_bin <- rbind(anchD,tarD) %>% select(cell,AllChr,AllSt, zscore)
#separating pos and neg zscores
SNP_inters_bin_pos <- SNP_inters_bin %>%
  group_by(AllChr,AllSt,cell) %>%
  filter(zscore >0) %>%
  summarise(mPosZscore = mean(zscore), .groups = "keep")
head(SNP_inters_bin_pos)
SNP_inters_bin_neg <- SNP_inters_bin %>%
  group_by(AllChr,AllSt,cell) %>%
  filter(zscore <0) %>%
  summarise(mNegZscore = mean(zscore), .groups = "keep")
head(SNP_inters_bin_neg)
SNP_inters_bin_all <- merge(SNP_inters_bin_pos, SNP_inters_bin_neg, by = c("AllChr","AllSt","cell"))
head(SNP_inters_bin_all)

SNP_inters_count <- rbind(anchD,tarD) %>% select(cell,AllChr,AllSt) %>%
  group_by(AllChr,AllSt,cell) %>%
  summarise(numInters = n(), .groups = "keep")
head(SNP_inters_count)

#read SNP files
regSNPs <- read.table(regSNPs_file, sep = "\t", header = T)
head(regSNPs)
print("#bins that SNPs are in")
regSNPs$bin <- plyr::round_any(regSNPs$start_b38, as.numeric(as.character(bin_size)), f = floor)
regSNPs_count <- regSNPs %>% dplyr::select(chr_b38,start_b38,end_b38,bin) %>%
  dplyr::group_by(chr_b38,bin) %>%
  dplyr::summarise(reg = n(), .groups = "keep")
head(regSNPs_count)
#nonreg
nonRegSNPs <- read.table(nonRegSNPs_file, sep = "\t", header = T)
head(nonRegSNPs)
print("#bins that SNPs are in")
nonRegSNPs$bin <- plyr::round_any(nonRegSNPs$start_b38, as.numeric(as.character(bin_size)), f = floor)
nonRegSNPs_count <- nonRegSNPs %>% dplyr::select(chr_b38,start_b38,end_b38,bin) %>%
  dplyr::group_by(chr_b38,bin) %>%
  dplyr::summarise(nonReg = n(), .groups = "keep")
head(nonRegSNPs_count)
#merge reg and non reg SNPs counts per bin
SNPs_count <- full_join(regSNPs_count,nonRegSNPs_count, by = c("chr_b38","bin"))
#add constant to reg and non-reg cols
SNPs_count <- SNPs_count %>% mutate(reg = reg +1) %>% mutate(nonReg = nonReg +1) %>%
              replace(is.na(.),1)
#calculate ratio of reg:non-reg SNPs per bin
SNPs_count <- SNPs_count %>% mutate(ratio = reg/nonReg)
head(SNPs_count)
summary(SNPs_count)

#combine ratio with num of interactions
colnames(SNP_inters_count) <- c("chr_b38","bin","cell","numInters")
count_SNPs_inters <- merge(SNPs_count,SNP_inters_count, by = c("chr_b38","bin"))
head(count_SNPs_inters)
summary(count_SNPs_inters)


#plot and correlation per cell
for (c in cells_sub$V1){
  #c = "Primitive_cardiomyocyte_day15_Zhang"
  print("##############")
  print(c)
  print("##############")
  tdat <- count_SNPs_inters %>% filter(cell == c)
  #line plot of data
  p <- (ggplot(tdat, aes(y=numInters, x=ratio))
        + geom_point(alpha=0.5)
        + geom_smooth(method = 'loess', formula = y ~ x)
        #+ geom_text(stat='count', aes(label=after_stat(count)), vjust=-1,hjust=-1,size=3,angle = 90)
        #+ geom_text(stat='count', aes(label=after_stat(count)),hjust=-0.2,size=3,angle = 90)
        #+ scale_x_upset(n_intersections = 90)
        #+ scale_x_upset(order_by = "freq")
        #+ scale_y_continuous(breaks = NULL, name = "", lim = c(0,30000))
        + labs(y = "# of cis-interactions",
               x = "# of regulatory SNPs : # of non-regulatory SNPs",
               title = c)
  )
  pdf(paste0(outprefix,"_",c,"_num_inters_line_pt_plot.pdf"), width = 14, height = 4)
  print(p)
  dev.off()
  print("# pearson correlation test between # inters and # of SNPs ratio")
  print(cor.test(tdat$ratio,tdat$numInters,method="pearson"))
}

#combine ratio with mean zscore of interactions
colnames(SNP_inters_bin_all) <- c("chr_b38","bin","cell","mPosZscore","mNegZscore")
count_SNPs_inters_zscore <- merge(SNPs_count,SNP_inters_bin_all, by = c("chr_b38","bin"))
count_SNPs_inters_zscore <- count_SNPs_inters_zscore %>% gather(key = "Zsign", value = "zscore", 7:8)
#re-order levels for pos and neg zscores
count_SNPs_inters_zscore <- count_SNPs_inters_zscore %>%
  mutate(Zsign = fct_relevel(Zsign,"mPosZscore", "mNegZscore"))
head(count_SNPs_inters_zscore)
summary(count_SNPs_inters_zscore)

#plot and correlation per cell
for (c in cells_sub$V1){
  #c = "Primitive_cardiomyocyte_day15_Zhang"
  #c = "OmniC_pooled_M_0d"
  print("##############")
  print(c)
  print("##############")
  tdat <- count_SNPs_inters_zscore %>% filter(cell == c)
  head(tdat)
#line plot of data
p <- (ggplot(tdat, aes(y=zscore, x=ratio))
      + geom_point(alpha=0.5)
      + geom_smooth(method = 'loess', formula = y ~ x)
      + facet_grid(Zsign ~ ., scales = "free")
      #+ geom_text(stat='count', aes(label=after_stat(count)), vjust=-1,hjust=-1,size=3,angle = 90)
      #+ geom_text(stat='count', aes(label=after_stat(count)),hjust=-0.2,size=3,angle = 90)
      #+ scale_x_upset(n_intersections = 90)
      #+ scale_x_upset(order_by = "freq")
      #+ scale_y_continuous(breaks = NULL, name = "", lim = c(0,30000))
      + labs(y = "z-score",
             x = "# of regulatory SNPs : # of non-regulatory SNPs",
             title = c)
)
pdf(paste0(outprefix,"_",c,"_zscore_line_pt_plot.pdf"), width = 14, height = 4)
print(p)
dev.off()
print("# pearson correlation test between # inters and # of SNPs ratio POS zscores")
pdat <- tdat %>% filter(Zsign == "mPosZscore")
print(cor.test(pdat$ratio,pdat$zscore,method="pearson"))
ndat <- tdat %>% filter(Zsign == "mNegZscore")
if (nrow(tdat) != 0){
print("# pearson correlation test between # inters and # of SNPs ratio NEG zscores")
print(cor.test(ndat$ratio,ndat$zscore,method="pearson"))
}#if
}


#regSNPs_inters <- read.table(regSNPs_intersfile, header = F)
#colnames(regSNPs_inters) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
##filter by desired cells
#regSNPs_inters <- regSNPs_inters %>% filter(cell %in% cells_sub$V1)
#head(regSNPs_inters)
##num inters per bin
#print("#counting each interaction twice (once for each chrom in interaction)")
#anchD <- regSNPs_inters
#anchD$AllChr <- anchD$chrA
#anchD$AllSt <- anchD$st1
#anchD$AllEnd <- anchD$end1
#tarD <- regSNPs_inters
#tarD$AllChr <- tarD$chrB
#tarD$AllSt <- tarD$st2
#tarD$AllEnd <- tarD$end2
#regSNPs_inters_bin <- rbind(anchD,tarD) %>% select(cell,AllChr,AllSt, zscore)
##separating pos and neg zscores
#regSNPs_inters_bin_pos <- regSNPs_inters_bin %>%
#  group_by(AllChr,AllSt,cell) %>%
#  filter(zscore >0) %>%
#  summarise(mPosZscore = mean(zscore), .groups = "keep")
#head(regSNPs_inters_bin_pos)
#regSNPs_inters_bin_neg <- regSNPs_inters_bin %>%
#  group_by(AllChr,AllSt,cell) %>%
#  filter(zscore <0) %>%
#  summarise(mNegZscore = mean(zscore), .groups = "keep")
#head(regSNPs_inters_bin_neg)
#regSNPs_inters_bin_all <- merge(regSNPs_inters_bin_pos, regSNPs_inters_bin_neg, by = c("AllChr","AllSt","cell"))
#regSNPs_inters_bin_all$SNP <- rep("reg",nrow(regSNPs_inters_bin_all))
#head(regSNPs_inters_bin_all)
#
#regSNPs_inters_count <- rbind(anchD,tarD) %>% select(cell,AllChr,AllSt) %>%
#                        group_by(AllChr,AllSt,cell) %>%
#                        summarise(numRegInters = n(), .groups = "keep")
#head(regSNPs_inters_count)
#
#nonRegSNPs_inters <- read.table(nonRegSNPs_intersfile, header = F)
#colnames(nonRegSNPs_inters) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
##filter by desired cells
#nonRegSNPs_inters <- nonRegSNPs_inters %>% filter(cell %in% cells_sub$V1)
#head(nonRegSNPs_inters)
##merge reg and non-reg inters
#SNP_inters <- full_join(regSNPs_inters,nonRegSNPs_inters, by = c("chrA","st1","end1","chrB","st2","end2","cell","zscore"))
#head(SNP_inters)
##num inters per bin
#print("#counting each interaction twice (once for each chrom in interaction)")
#anchD <- nonRegSNPs_inters
#anchD$AllChr <- anchD$chrA
#anchD$AllSt <- anchD$st1
#anchD$AllEnd <- anchD$end1
#tarD <- nonRegSNPs_inters
#tarD$AllChr <- tarD$chrB
#tarD$AllSt <- tarD$st2
#tarD$AllEnd <- tarD$end2
#nonRegSNPs_inters_bin <- rbind(anchD,tarD) %>% select(cell,AllChr,AllSt, zscore)
##separating pos and neg zscores
#nonRegSNPs_inters_bin_pos <- nonRegSNPs_inters_bin %>%
#  group_by(AllChr,AllSt,cell) %>%
#  filter(zscore >0) %>%
#  summarise(mPosZscore = mean(zscore), .groups = "keep")
#head(nonRegSNPs_inters_bin_pos)
#nonRegSNPs_inters_bin_neg <- nonRegSNPs_inters_bin %>%
#  group_by(AllChr,AllSt,cell) %>%
#  filter(zscore <0) %>%
#  summarise(mNegZscore = mean(zscore), .groups = "keep")
#head(nonRegSNPs_inters_bin_neg)
#nonRegSNPs_inters_bin_all <- merge(nonRegSNPs_inters_bin_pos, nonRegSNPs_inters_bin_neg, by = c("AllChr","AllSt","cell"))
#nonRegSNPs_inters_bin_all$SNP <- rep("nonReg",nrow(nonRegSNPs_inters_bin_all))
#head(nonRegSNPs_inters_bin_all)
#
##merge reg and non-reg inters
#SNP_inters_bin <- rbind(regSNPs_inters_bin_all,nonRegSNPs_inters_bin_all)
#head(SNP_inters_bin)
#
#
#nonRegSNPs_inters_count <- rbind(anchD,tarD) %>% select(cell,AllChr,AllSt) %>%
#  group_by(AllChr,AllSt,cell) %>%
#  summarise(numNonRegInters = n(), .groups = "keep")
#head(nonRegSNPs_inters_count)
##merge num inters dfs
#SNP_inters_count <- full_join(regSNPs_inters_count,nonRegSNPs_inters_count, by = c("AllChr","AllSt","cell"))
#wide to long format
#SNP_inters_count <- gather(SNP_inters_count, key="SNP", value="numInters", 4:5)
#head(SNP_inters_count)
#
##read SNP files
#regSNPs <- read.table(regSNPs_file, sep = "\t", header = T)
#head(regSNPs)
#print("#bins that SNPs are in")
#regSNPs$bin <- plyr::round_any(regSNPs$start_b38, as.numeric(as.character(bin_size)), f = floor)
#regSNPs_count <- regSNPs %>% dplyr::select(chr_b38,start_b38,end_b38,bin) %>%
#                 dplyr::group_by(chr_b38,bin) %>%
#                 dplyr::summarise(reg = n(), .groups = "keep")
#head(regSNPs_count)
##nonreg
#nonRegSNPs <- read.table(nonRegSNPs_file, sep = "\t", header = T)
#head(nonRegSNPs)
#print("#bins that SNPs are in")
#nonRegSNPs$bin <- plyr::round_any(nonRegSNPs$start_b38, as.numeric(as.character(bin_size)), f = floor)
#nonRegSNPs_count <- nonRegSNPs %>% dplyr::select(chr_b38,start_b38,end_b38,bin) %>%
#                 dplyr::group_by(chr_b38,bin) %>%
#                 dplyr::summarise(nonReg = n(), .groups = "keep")
#head(nonRegSNPs_count)
##merge reg and non reg SNPs counts per bin
#SNPs_count <- merge(regSNPs_count,nonRegSNPs_count, by = c("chr_b38","bin"))
##calculate ratio of reg:non-reg SNPs per bin
#SNPs_count <- SNPs_count %>% mutate(ratio = reg/nonReg)
#head(SNPs_count)
#summary(SNPs_count)
#
##combine ratio with num of interactions
#colnames(SNP_inters_count) <- c("chr_b38","bin","cell","SNP","numInters")
#count_SNPs_inters <- merge(SNPs_count,SNP_inters_count, by = c("chr_b38","bin"))
#head(count_SNPs_inters)
#summary(count_SNPs_inters)
#
##line plot of data
#p <- (ggplot(count_SNPs_inters, aes(y=numInters, x=ratio, fill = SNP, colour = SNP))
#       + geom_point()
#       + geom_smooth(method = 'loess', formula = y ~ x)
#       #+ geom_text(stat='count', aes(label=after_stat(count)), vjust=-1,hjust=-1,size=3,angle = 90)
#       #+ geom_text(stat='count', aes(label=after_stat(count)),hjust=-0.2,size=3,angle = 90)
#       #+ scale_x_upset(n_intersections = 90)
#       #+ scale_x_upset(order_by = "freq")
#       #+ scale_y_continuous(breaks = NULL, name = "", lim = c(0,30000))
#       + labs(y = "# of significant Cis-chromosomal interactions overlapping with SNPs",
#              x = "# of regulatory SNPs : # of non-regulatory SNPs",
#              title = "")
#)
#pdf(paste0(outprefix,"_num_inters_line_pt_plot.pdf"), width = 14, height = 4)
#p
#dev.off()
#print("# pearson correlation test between # inters and # of SNPs ratio: non-reg interactions")
#tdat <- count_SNPs_inters %>% filter(SNP == "numNonRegInters")
#cor.test(tdat$numInters,tdat$ratio,method="pearson")
#print("# pearson correlation test between # inters and # of SNPs ratio: reg interactions")
#tdat <- count_SNPs_inters %>% filter(SNP == "numRegInters")
#cor.test(tdat$numInters,tdat$ratio,method="pearson")
#
#
##combine ratio with mean zscore of interactions
#colnames(SNP_inters_bin) <- c("chr_b38","bin","cell","mPosZscore","mNegZscore", "SNP")
#count_SNPs_inters_zscore <- merge(SNPs_count,SNP_inters_bin, by = c("chr_b38","bin"))
#count_SNPs_inters_zscore <- count_SNPs_inters_zscore %>% gather(key = "Zsign", value = "zscore", 7:8)
##re-order levels for pos and neg zscores
#count_SNPs_inters_zscore <- count_SNPs_inters_zscore %>%
#  mutate(Zsign = fct_relevel(Zsign,"mPosZscore", "mNegZscore"))
#head(count_SNPs_inters_zscore)
#summary(count_SNPs_inters_zscore)
#
##line plot of data
#p <- (ggplot(count_SNPs_inters_zscore, aes(y=zscore, x=ratio, fill = SNP, colour = SNP))
#      + geom_point()
#      + geom_smooth(method = 'loess', formula = y ~ x)
#      + facet_grid(Zsign ~ ., scales = "free")
#      #+ geom_text(stat='count', aes(label=after_stat(count)), vjust=-1,hjust=-1,size=3,angle = 90)
#      #+ geom_text(stat='count', aes(label=after_stat(count)),hjust=-0.2,size=3,angle = 90)
#      #+ scale_x_upset(n_intersections = 90)
#      #+ scale_x_upset(order_by = "freq")
#      #+ scale_y_continuous(breaks = NULL, name = "", lim = c(0,30000))
#      #+ labs(x = "# of significant Cis-chromosomal interactions overlapping with SNPs",
#      #       y = "# of regulatory SNPs : # of non-regulatory SNPs",
#      #       title = "")
#)
#pdf(paste0(outprefix,"_zscore_line_pt_plot.pdf"), width = 14, height = 4)
#p
#dev.off()
#print("# pearson correlation test between # inters and # of SNPs ratio: non-reg interactions POS zscores")
#tdat <- count_SNPs_inters_zscore %>% filter(SNP == "nonReg" & Zsign == "mPosZscore")
#cor.test(tdat$zscore,tdat$ratio,method="pearson")
#print("# pearson correlation test between # inters and # of SNPs ratio: non-reg interactions NEG zscores")
#tdat <- count_SNPs_inters_zscore %>% filter(SNP == "nonReg" & Zsign == "mNegZscore")
#cor.test(tdat$zscore,tdat$ratio,method="pearson")
#print("# pearson correlation test between # inters and # of SNPs ratio: reg interactions POS zscores")
#tdat <- count_SNPs_inters_zscore %>% filter(SNP == "reg" & Zsign == "mPosZscore")
#cor.test(tdat$zscore,tdat$ratio,method="pearson")
#print("# pearson correlation test between # inters and # of SNPs ratio: reg interactions NEG zscores")
#tdat <- count_SNPs_inters_zscore %>% filter(SNP == "reg" & Zsign == "mNegZscore")
#cor.test(tdat$zscore,tdat$ratio,method="pearson")

print("DONE")
