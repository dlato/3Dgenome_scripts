########################################
# all genes present in common trans-chromosomal interactions
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: CM SNP file (tsv)
#            VSMC SNP file (tsv)
#            output file prefix
#            bin size (numeric, bp)
########################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
CM_SNP_file <- args[1]
VSMC_SNP_file <- args[2]
outprefix <- args[3]
bin_size <- as.numeric(as.character(args[4]))

##########
.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.4.1.2")
library(tidyr)
library(dplyr)
library(plyr)
library(ggplot2)
library(karyoploteR)#for karyotype plot
library(GenomicRanges)
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
####interaction data
#options(scipen = 999)
#bin_size <- 1000000
#CM_SNP_file <- "CM_test_SNPs.txt"
#outfile <- "test_num_snps"
#VSMC_SNP_file <- "VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt"
#library(tidyr)#for PCA
#library(dplyr) #for colours
#library(karyoploteR)#for PCA
#library(GenomicRanges)#for plotting

#read in SNP info
CM_SNP <- read.table(CM_SNP_file, header = TRUE,sep = "\t")
print("#bins that SNPs are in")
CM_SNP$bin <- plyr::round_any(CM_SNP$start_b38, as.numeric(as.character(bin_size)), f = floor)
CM_num_df <- CM_SNP %>% dplyr::select(chr_b38, bin) %>% group_by(chr_b38,bin) %>% dplyr::summarise(CM_num = n(), .groups = "keep") 
colnames(CM_num_df) <- c("chrom", "start","CM_num")
CM_num_df$end <- CM_num_df$start + bin_size
CM_num_df <- CM_num_df %>% dplyr::select(chrom, start, end, CM_num)
VSMC_SNP <- read.table(VSMC_SNP_file, header = TRUE,sep = "\t")
print("#bins that SNPs are in")
VSMC_SNP$bin <- plyr::round_any(VSMC_SNP$start_b38, as.numeric(as.character(bin_size)), f = floor)
VSMC_num_df <- VSMC_SNP %>% dplyr::select(chr_b38, bin) %>% group_by(chr_b38,bin) %>% dplyr::summarise(VSMC_num = n(), .groups = "keep") 
colnames(VSMC_num_df) <- c("chrom", "start","VSMC_num")
VSMC_num_df$end <- VSMC_num_df$start + bin_size
VSMC_num_df <- VSMC_num_df %>% dplyr::select(chrom, start, end, VSMC_num)
#merge dfs
dat <- merge(VSMC_num_df, CM_num_df, by = c("chrom", "start", "end"))
write.table(dat, file = as.character(paste0(outprefix,"_num_SNPs_per_",bin_size,"bp_bin.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#convert to GRanges
gdat <- toGRanges(as.data.frame(dat))
head(gdat)

###############
# karyotype of number of SNPs per bin
##############
pdf(paste0(outprefix,"_num_SNPs_per_",bin_size,"bp_bin_karyotype.pdf"), width = 14, height = 8)
kp <- plotKaryotype(genome = "hg38", plot.type = 1)
kpAddMainTitle(kp, paste0("Number of variants per ",bin_size,"bp bin"), cex=0.8)
kpAddBaseNumbers(kp)
at <- autotrack(current.track = 1, total.tracks = 2, margin = 0.2)
kpAddLabels(kp, labels = "VSMC", r0=at$r0, r1=at$r1, cex = 0.7)
kpHeatmap(kp, data=gdat,y=gdat$VSMC_num, colors = c("#FFD199","#FF8C00"),r0=0.355, r1=0.66)
at <- autotrack(current.track = 2, total.tracks = 2, margin = 0.2)
kpAddLabels(kp, labels = "CM", r0=at$r0, r1=at$r1, cex = 0.7)
kpHeatmap(kp, data=gdat,y=gdat$CM_num, colors = c("#C3BACA","#800080"),r0=0, r1=0.305)
dev.off()

for (c in unique(dat$chrom)){
  #c = "chr1"
  pdf(paste0(outprefix,"_num_SNPs_per_",bin_size,"bp_bin_",c,"_karyotype.pdf"), width = 14, height = 4)
  kp <- plotKaryotype(genome = "hg38", plot.type = 1, chromosomes = c(c))
  kpAddMainTitle(kp, paste0("Number of variants per ",bin_size,"bp bin"), cex=0.8)
  kpAddBaseNumbers(kp)
  at <- autotrack(current.track = 1, total.tracks = 2, margin = 0.2)
  kpAddLabels(kp, labels = "VSMC", r0=at$r0, r1=at$r1, cex = 0.7)
  kpHeatmap(kp, data=gdat,y=gdat$VSMC_num, colors = c("#C3BACA","#800080"),r0=0, r1=0.45)
  at <- autotrack(current.track = 2, total.tracks = 2, margin = 0.2)
  kpAddLabels(kp, labels = "CM", r0=at$r0, r1=at$r1, cex = 0.7)
  kpHeatmap(kp, data=gdat,y=gdat$CM_num, colors = c("#FFD199","#FF8C00"),r0=0.5, r1=1)
  dev.off()
}


# chrom info
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
p_chr_ord_gsub <- gsub("chr","",p_chr_ord)

#regular tickplot of number of SNPs per bin
head(dat)
for (c in unique(dat$chrom)){
  #c = "chr1"
  c_name <- gsub("chr","",c)
  xmax <- chrInf %>% filter(chrom == c) %>% dplyr::select(size)
  #VSMC
  tdat <- dat %>% filter(chrom == c)
  tp <- (ggplot(tdat, aes(start/1000000, y=1, fill = VSMC_num))
       + geom_tile(aes(fill = VSMC_num), width = 1, height = 1)
#       + scale_fill_gradient(low = "#FFD199", high = "#FF8C00", name = "# of vairants per bin")
       + scale_fill_gradient(low = "#B8B8B8", high = "#000000", name = "# of vairants per bin")
       #+ scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "OR genes per bin", na.value = 0)
       + labs(x = paste0("Chromosome ", c_name, " position [Mb]"),
              y = "",
              title = paste0("Number of variants per ",bin_size,"bp bin"))
       + expand_limits(x = c(0,xmax[1,1]/1000000))
       + scale_y_continuous(expand = c(0, 0))
       + scale_x_continuous(expand = c(0, 0))
       + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
               legend.text = element_text(size = 10),
               strip.background = element_rect(fill = "white"),
               panel.spacing = unit(0, "lines"),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
       + theme(panel.background = element_rect(fill = "white", colour = NA))
)
filename <- paste0(outprefix,"_VSMC_num_SNPs_per_",bin_size,"bp_bin_",c,"_tickplot.pdf")
pdf(filename, width = 14, height = 2)
print(tp)
dev.off()
#CM
tp <- (ggplot(tdat, aes(start/1000000, y=1, fill = CM_num))
       + geom_tile(aes(fill = CM_num), width = 1, height = 1)
#       + scale_fill_gradient(low = "#C3BACA", high = "#800080", name = "# of vairants per bin")
       + scale_fill_gradient(low = "#B8B8B8", high = "#000000", name = "# of vairants per bin")
       #+ scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "OR genes per bin", na.value = 0)
       + labs(x = paste0("Chromosome ", c_name, " position [Mb]"),
              y = "",
              title = paste0("Number of variants per ",bin_size,"bp bin"))
       + expand_limits(x = c(0,xmax[1,1]/1000000))
       + scale_y_continuous(expand = c(0, 0))
       + scale_x_continuous(expand = c(0, 0))
       + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
               legend.text = element_text(size = 10),
               strip.background = element_rect(fill = "white"),
               panel.spacing = unit(0, "lines"),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
       + theme(panel.background = element_rect(fill = "white", colour = NA))
)
filename <- paste0(outprefix,"_CM_num_SNPs_per_",bin_size,"bp_bin_",c,"_tickplot.pdf")
pdf(filename, width = 14, height = 2)
print(tp)
dev.off()
}
print("DONE")
#
