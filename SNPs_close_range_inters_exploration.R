########################################
# looking at close range (<20Kb) interactions to see where the are located (blood pressure project)
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: z-scrore interactions (from 3Dflow output, ALL inters (sig and non-sig))
#            p-values of interactions (from 3Dflow output, ALL inters (sig and non-sig))
#            cells to look for SNPs in (text file with each cell name on a different row)
#            output file prefix
#            SNP file (tsv)
#	     bin size (bp)
########################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
zdat_file <- args[1]
pdat_file <- args[2]
cells_file <- args[3]
outprefix <- args[4]
SNP_file <- args[5]
bin_size <- as.numeric(as.character(args[6]))

##########
library(tidyr)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(ggplot2)
library(ggforce)#for ridgeline
library(ggridges)#for ridgeline
#library(ggbiplot)#for PCA
#library(devtools)#for PCA
.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2")
library(circlize,lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2/") # for circos
library(regioneR,lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2/")#for permutation
library(factoextra, lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3/")#for PCA
library(harrypotter) #for colours
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
####interaction data
##options(scipen = 999)
#Atype <- "1_vs_All"
##tissue_file <- "tissue_system_info.txt"
#zdat_file <- "test_1vsAll_dat.txt"
#pdat_file <- "test_1vsAll_pvalues.txt"
#cells_file <- "cell_subset.txt"
##germlayer_file <- "germlayer_info.txt"
#bin_size <- 1000000
#anno_file <- "hg38_p13_v32_annotation.txt"
#outfile <- "test_cell_list"
#SNP_file <- "VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt"
#library(factoextra)#for PCA
#library(harrypotter) #for colours

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
#read in cells to filter
cells_sub <- read.table(cells_file)
cells_sub
#read in SNP info
SNP <- read.table(SNP_file, header = TRUE,sep = "\t")
#accounting for merged SNP files
SNP_df <- SNP
if ("logFC_comp" %in% colnames(SNP)){
  SNP_df <- SNP %>% dplyr::select(chr_b38, start_b38,end_b38, logFC_comp)
} else {
  SNP_df$logFC_comp.x <- as.numeric(as.character(SNP_df$logFC_comp.x))
  SNP_df$logFC_comp.y <- as.numeric(as.character(SNP_df$logFC_comp.y))
  SNP_df$logFC_comp <- mean(c(abs(SNP_df$logFC_comp.x),abs(SNP_df$logFC_comp.y)))
}#if else
head(SNP_df)
#z-scores
zdat <- read.table(zdat_file, header = TRUE)
print("#split ID col")
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
zdat$ID <- sub("B", "\\.B", as.character(zdat$ID))
zdat <- zdat %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
print("#remove A and B from chrom names")
zdat$chrA <- gsub("A", "", zdat$chrA)
zdat$chrB <- gsub("B", "", zdat$chrB)
print("#wide to long format")
lzdat <- zdat %>% gather(cell, zscore, 8:ncol(zdat)) %>% na.omit() %>% filter(cell %in% cells_sub$V1)
head(lzdat)
#p-values
pdat <- read.table(pdat_file, header = TRUE)
print("#split ID col")
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
pdat$ID <- sub("B", "\\.B", as.character(pdat$ID))
pdat <- pdat %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
print("#remove A and B from chrom names")
pdat$chrA <- gsub("A", "", pdat$chrA)
pdat$chrB <- gsub("B", "", pdat$chrB)
print("#wide to long format")
lpdat <- pdat %>% gather(cell, pvalue, 8:ncol(pdat))
head(lpdat) %>% na.omit() %>% filter(cell %in% cells_sub$V1)
print("merge zscore and pvals")
ldat <- merge(lzdat,lpdat) %>% na.omit() %>% filter(cell %in% cells_sub$V1)
head(ldat)
print("summary of ALL zscores and pvalues")
summary(ldat)


zdat <- data.frame()
pdat <- data.frame()
print("#filter interactions for ones that contain SNPs (for circos plot)")
print("#bins that SNPs are in")
summary(SNP_df)
SNP_df$start_b38 <- as.numeric(as.character(SNP_df$start_b38))
summary(SNP_df)
SNP_df$bin <- plyr::round_any(SNP_df$start_b38, bin_size), f = floor)
colnames(SNP_df) <- c("AllChr","SNPstart","SNPend","logFC","AllSt")
summary(SNP_df)
SNP_ID_df <- SNP_df
SNP_ID_df$ID <- format(SNP_ID_df$ID, scientific = FALSE)
SNP_ID_df$ID <- paste0(SNP_ID_df$AllChr,".",format(SNP_ID_df$AllSt, scientific = FALSE),".",format(SNP_ID_df$AllSt + as.numeric(as.character(bin_size)),scientific = FALSE))
#unique chroms containing SNPs
uSNP_chrs<- unique(SNP_ID_df$AllChr)
#filter interactions for just unique chroms in SNPs (to decrease mem)
ldat <- ldat %>% filter(chrA %in% uSNP_chrs | chrB %in% uSNP_chrs)
head(ldat)
SNP_ID_df$ID <- gsub(" ","", SNP_ID_df$ID)
head(SNP_ID_df)
#decrease SNP list to just unique IDs, find mean of |logFC|
SNP_uniq <- SNP_ID_df %>%
  group_by(ID) %>%
  dplyr::summarise(mlogFC = mean(abs(logFC)))
head(SNP_uniq)
print("done formatting ID and bin col")
print("finding SNPs that overlap with inters for circos plot")
selection <- as.logical( # 7
  rowSums( # 6
    matrix( # 5
      unlist( # 4
        lapply(SNP_uniq$ID, function(x) { #3
          as.numeric( # 2
            str_detect(ldat$ID,x) # 1
          )
        })
      ), 
      ncol = length(SNP_uniq$ID), byrow = F)
  )
)

odat <- unique(ldat[selection,] %>% dplyr::select(-ID)) %>% na.omit()
head(odat)
print("#write to df")
write.table(odat, file = as.character(paste0(outprefix,"_overlapping_intearctions_zscore_pvalue.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#calculate distance between each interaction
odat$st1 <- as.numeric(as.character(odat$st1))
odat$st2 <- as.numeric(as.character(odat$st2))
odat$dist <- abs(odat$st1 - odat$st2)
#col for sig or non-sig
odat <- odat %>%
  mutate(Sig = case_when(
    pvalue <= 0.05 ~ "sig",
    pvalue >0.05 ~ "nonSig"
  ))
#col for short or long range
odat <- odat %>%
  mutate(range = case_when(
    dist <= 20000 ~ "short",
    dist > 20000 ~ "long"
  ))
#col for short or long range
odat <- odat %>%
  mutate(sign = case_when(
    zscore < 0 ~ "pos",
    zscore > 0 ~ "neg"
  ))
head(odat)

ninters_sig <- odat %>%
  group_by(cell,Sig, range) %>%
  dplyr::summarise(n = n(), .groups = "keep")
ninters_sig
ninters_sign <- odat %>% filter(Sig == "sig")%>%
  group_by(cell,sign, range) %>%
  dplyr::summarise(n = n(), .groups = "keep")
ninters_sign

print("# boxplot/violin of sig and non-sig overlapping inters")
p <- (ggplot(ninters_sig, aes(x=range, y=n) )
      + geom_violin(alpha = 0.6, fill = "grey")
      + geom_boxplot(alpha =0.6, fill = "grey")
      + labs(title = "Distance between interacting regions",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = "Interaction range",
             y = "Number of Interactions")
      + facet_grid( .~ Sig)
      + scale_x_discrete(expand = c(0, 0))
      + scale_y_continuous(expand = c(0, 0))
)
f_name <- gsub(" ","",paste(outprefix,"_cis_short_long_range_interactions_sig_boxplot.pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()

print("# boxplot/violin of sig pos and neg overlapping inters")
p <- (ggplot(ninters_sign, aes(x=range, y=n) )
      + geom_violin(alpha = 0.6, fill = "grey")
      + geom_boxplot(alpha =0.6, fill = "grey")
      + labs(title = "Distance between interacting regions",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = "Interaction range",
             y = "Number of Interactions")
      + facet_grid( .~ sign)
      + scale_x_discrete(expand = c(0, 0))
      + scale_y_continuous(expand = c(0, 0))
)
f_name <- gsub(" ","",paste(outprefix,"_cis_short_long_range_interactions_zscore_sign_boxplot.pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()

print("#test between num of inters btwn short and long range interactions SIG INTERS")
print("#Perform the Mann-Whitney test (when dists are not normal)")
print("sig = mean is diff btwn short and long range")
wilcox.test(n~range, data=ninters_sig %>% filter(Sig == "sig"))
print("#test between num of inters btwn short and long range interactions NON-SIG INTERS")
print("#Perform the Mann-Whitney test (when dists are not normal)")
print("sig = mean is diff btwn short and long range")
wilcox.test(n~range, data=ninters_sig %>% filter(Sig == "nonSig"))

print("#test between num of inters btwn short and long range interactions SIG POS INTERS")
print("#Perform the Mann-Whitney test (when dists are not normal)")
print("sig = mean is diff btwn short and long range")
wilcox.test(n~range, data=ninters_sign %>% filter(sign == "pos"))
print("#test between num of inters btwn short and long range interactions SIG NEG INTERS")
print("#Perform the Mann-Whitney test (when dists are not normal)")
print("sig = mean is diff btwn short and long range")
wilcox.test(n~range, data=ninters_sign %>% filter( sign == "neg"))

print("DONE")
#
