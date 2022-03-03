########################################
# all genes present in common trans-chromosomal interactions
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: significant z-scrore interactions (from 3Dflow output)
#            cells to look for SNPs in (text file with each cell name on a different row)
#            annotation file (from python script) (tsv)
#            output file prefix
#            SNP file (tsv)
#            bin size (numeric, bp)
########################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
cells_file <- args[2]
anno_file <- args[3]
outfile <- args[4]
SNP_file <- args[5]
bin_size <- args[6]

##########
library(tidyr)
library(dplyr)
library(stringr)
library(GenomicRanges)
library(ggplot2)
library(ggforce)#for ridgeline
library(ggridges)#for ridgeline
library(ggbiplot)#for PCA
library(devtools)#for PCA
#library(multcomp) # for anova and tukey test
.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2")
library(circlize,lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2/") # for circos
library(regioneR,lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2/")#for permutation
library(factoextra, lib = "/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3/")#for PCA
library(harrypotter, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.3") #for colours
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
#dat_file <- "test_1vsAll_dat.txt"
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
dat <- read.table(dat_file, header = TRUE)
print("summary of ALL sig zscores per cell type")
summary(dat)
print("#split ID col")
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
dat$ID <- sub("B", "\\.B", as.character(dat$ID))
dat <- dat %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
print("#remove A and B from chrom names")
dat$chrA <- gsub("A", "", dat$chrA)
dat$chrB <- gsub("B", "", dat$chrB)
print("#wide to long format")
ldat <- dat %>% gather(cell, zscore, 8:ncol(dat))

dat <- data.frame()
print("#filter interactions for ones that contain SNPs (for circos plot)")
print("#bins that SNPs are in")
SNP_df$bin <- plyr::round_any(SNP_df$start_b38, as.numeric(as.character(bin_size)), f = floor)
colnames(SNP_df) <- c("AllChr","SNPstart","SNPend","logFC","AllSt")
#print("#df for number of SNPs per bin")
#num_SNPs <- SNP_df %>% dplyr::select(AllChr,AllSt) %>% group_by(AllChr,AllSt) %>% dplyr::summarise(numSNPs = n(), .groups = "keep") 
#head(num_SNPs)
SNP_ID_df <- SNP_df
SNP_ID_df$ID <- format(SNP_ID_df$ID, scientific = FALSE)
SNP_ID_df$ID <- paste0(SNP_ID_df$AllChr,".",format(SNP_ID_df$AllSt, scientific = FALSE),".",format(SNP_ID_df$AllSt + as.numeric(as.character(bin_size)),scientific = FALSE))
#unique chroms containing SNPs
uSNP_chrs<- unique(SNP_ID_df$AllChr)
#filter interactions for just unique chroms in SNPs (to decrease mem)
ldat <- ldat %>% filter(chrA %in% uSNP_chrs | chrB %in% uSNP_chrs)
head(ldat)
SNP_ID_df$ID <- gsub(" ","", SNP_ID_df$ID)
#decrease SNP list to just unique IDs, find mean of |logFC|
SNP_uniq <- SNP_ID_df %>%
  group_by(ID) %>%
  dplyr::summarise(mlogFC = mean(abs(logFC)))
head(SNP_uniq)
print("done formatting ID and bin col")
print("finding inters that have SNPs in each inter region")
ldat <- ldat %>% mutate(IDa = paste0(chrA,".",st1,".",end1))
ldat <- ldat %>% mutate(IDb = paste0(chrB,".",st2,".",end2))
head(ldat)

odat <- ldat %>% filter(IDa %in% unique(SNP_uniq$ID) & IDb %in% unique(SNP_uniq$ID)) %>%
        dplyr::select(-ID,-IDa,-IDb)
head(odat)
print("#write to df")
write.table(odat, file = as.character(paste0(outfile,"_overlapping_intearctions_for_circos.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


print("#counting each interaction twice (once for each chrom in interaction)")
anchD <- odat
anchD$AllChr <- anchD$chrA
anchD$AllSt <- anchD$st1
anchD$AllEnd <- anchD$end1
tarD <- odat
tarD$AllChr <- tarD$chrB
tarD$AllSt <- tarD$st2
tarD$AllEnd <- tarD$end2
r_dat2 <- rbind(anchD,tarD)
#remove NA interactions
noNAdat2 <- r_dat2 %>% select(cell, zscore, AllChr, AllSt, AllEnd) %>% na.omit() %>% filter(cell %in% cells_sub$V1)

print("check filtering by cell type")
unique(noNAdat2$cell)

#read in annotation file
anno_df <- read.table(anno_file, header = TRUE)
summary(anno_df)
head(anno_df)

#filter annotation for common interactions
colnames(noNAdat2) <- c("cell", "zscore","chr","st","end")
#positive z-scores
u_inters <- distinct(noNAdat2 %>% filter(zscore >0) %>% select(chr, st, end))
summary(u_inters)
common_genes <- c()
common_genes_metascape <- c()
for(i in 1:nrow(u_inters)) {
  #i=1
  td <- u_inters[i,]
  tgenes_df <- anno_df %>% filter(seqname == td$chr & td$st >= bin_start & td$st <= bin_end) %>%
    filter(broad_class =="prot")
  common_genes <- append(common_genes,gsub("\\..*", "",tgenes_df$gene_id, perl=TRUE))
  common_genes_metascape <- append(common_genes_metascape,gsub("\\..*", "",tgenes_df$gene_name, perl=TRUE))
} #for
common_genes <- unique(common_genes)
#common_genes
print("#### number of genes in common interactions")
length(common_genes)
write.table(common_genes, file = as.character(paste0(outfile,"_pos_zscore_common_inters_GO_analysis_gene_list.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(common_genes_metascape, file = as.character(paste0(outfile,"_pos_zscore_common_inters_metascape_analysis_gene_list.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#negative z-scores
u_inters <- distinct(noNAdat2 %>% filter(zscore <0) %>% select(chr, st, end))
summary(u_inters)
common_genes <- c()
common_genes_metascape <- c()
for(i in 1:nrow(u_inters)) {
  #i=1
  td <- u_inters[i,]
  tgenes_df <- anno_df %>% filter(seqname == td$chr & td$st >= bin_start & td$st <= bin_end) %>%
    filter(broad_class =="prot")
  common_genes <- append(common_genes,gsub("\\..*", "",tgenes_df$gene_id, perl=TRUE))
  common_genes_metascape <- append(common_genes_metascape,gsub("\\..*", "",tgenes_df$gene_name, perl=TRUE))
} #for
common_genes <- unique(common_genes)
#common_genes
print("#### number of genes in common interactions")
length(common_genes)
write.table(common_genes, file = as.character(paste0(outfile,"_neg_zscore_common_inters_GO_analysis_gene_list.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(common_genes_metascape, file = as.character(paste0(outfile,"_neg_zscore_common_inters_metascape_analysis_gene_list.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


print("# cell specific interactions")
for (c in cells_sub$V1){
  #filter per cell type
  tdat <- noNAdat2 %>% filter(cell == c) %>% na.omit()
  #positive z-scores
  u_inters <- distinct(tdat %>% filter(zscore >0) %>% select(chr, st, end))
  summary(u_inters)
  common_genes <- c()
  common_genes_metascape <- c()
  for(i in 1:nrow(u_inters)) {
    #i=1
    td <- u_inters[i,]
    tgenes_df <- anno_df %>% filter(seqname == td$chr & bin_start == td$st | bin_end == td$st) %>% 
      filter(broad_class =="prot")
    common_genes <- append(common_genes,gsub("\\..*", "",tgenes_df$gene_id, perl=TRUE))
    common_genes_metascape <- append(common_genes_metascape,gsub("\\..*", "",tgenes_df$gene_name, perl=TRUE))
  } #for
  common_genes <- unique(common_genes)
  #common_genes
  print("#### number of genes in common interactions")
  length(common_genes)
  write.table(common_genes, file = as.character(paste0(outfile,"_",c,"_pos_zscore_common_inters_GO_analysis_gene_list.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(common_genes_metascape, file = as.character(paste0(outfile,"_",c,"_pos_zscore_common_inters_metascape_analysis_gene_list.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  #negative z-scores
  u_inters <- distinct(tdat %>% filter(zscore <0) %>% select(chr, st, end))
  summary(u_inters)
  common_genes <- c()
  common_genes_metascape <- c()
  for(i in 1:nrow(u_inters)) {
    #i=1
    td <- u_inters[i,]
    tgenes_df <- anno_df %>% filter(seqname == td$chr & bin_start == td$st | bin_end == td$st) %>% 
      filter(broad_class =="prot")
    common_genes <- append(common_genes,gsub("\\..*", "",tgenes_df$gene_id, perl=TRUE))
    common_genes_metascape <- append(common_genes_metascape,gsub("\\..*", "",tgenes_df$gene_name, perl=TRUE))
  } #for
  common_genes <- unique(common_genes)
  #common_genes
  print("#### number of genes in common interactions")
  length(common_genes)
  write.table(common_genes, file = as.character(paste0(outfile,"_",c,"_neg_zscore_common_inters_GO_analysis_gene_list.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(common_genes_metascape, file = as.character(paste0(outfile,"_",c,"_neg_zscore_common_inters_metascape_analysis_gene_list.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
}

##########
## top 3000 genes in common interactions (roughly top 50 interactions)
## calculate the mean zscore per interaction (mean across all cells), then taking the top 79 interactions
## file comes from the trans_common_sig_interactions.R script
#topdat <- dat %>% arrange(desc(meanZscore))
#topdat
#topdat <- topdat[1:15,]
#topdat
#Adf <- topdat %>% select(chrA, st1,end1)
#colnames(Adf) <- c("chr","st","end")
#Bdf <- topdat %>% select(chrB, st2,end1)
#colnames(Bdf) <- c("chr","st","end")
#tu_inters <- distinct(rbind(Adf,Bdf))
#summary(tu_inters)
#tcommon_genes <- c()
#tcommon_genes_metascape <- c()
#for(i in 1:nrow(tu_inters)) {
#  #i=1
#  td <- tu_inters[i,]
#  tgenes_df <- anno_df %>% filter(seqname == td$chr & bin_start == td$st | bin_end == td$st) %>% 
#    filter(broad_class =="prot")
#  tcommon_genes <- append(tcommon_genes,gsub("\\..*", "",tgenes_df$gene_id, perl=TRUE))
#  tcommon_genes_metascape <- append(tcommon_genes_metascape,gsub("\\..*", "",tgenes_df$gene_name, perl=TRUE))
#} #for
#tcommon_genes <- unique(tcommon_genes)
#tcommon_genes
#print("#### number of genes in common interactions")
#length(tcommon_genes)
#write.table(tcommon_genes, file = as.character(paste0(outfile,"_top3000_common_inters_GO_analysis_gene_list.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#write.table(tcommon_genes_metascape, file = as.character(paste0(outfile,"_top3000_common_inters_metascape_analysis_gene_list.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#
#
#
## dealing with genes that are found in two bins: counting twice
#tdup_anno <- anno_df %>% filter(bin_start != bin_end)
#tdup_anno <- tdup_anno %>% select("seqname","source", "feature","start","end","score","strand","frame","gene_id","gene_type","gene_name","level","hgnc_id","havana_gene","transcript_id","transcript_type","transcript_name","transcript_support_level","tag","havana_transcript","exon_number","exon_id","ont","protein_id","ccdsid","broad_class","bin_end","bin_start")
#colnames(tdup_anno) <- c("seqname","source", "feature","start","end","score","strand","frame","gene_id","gene_type","gene_name","level","hgnc_id","havana_gene","transcript_id","transcript_type","transcript_name","transcript_support_level","tag","havana_transcript","exon_number","exon_id","ont","protein_id","ccdsid","broad_class","bin_start","bin_end")
#dup_anno <- rbind(anno_df,tdup_anno)
#
##filter by common interactions
#dup_anno$ID <- paste0(dup_anno$seqname,".",dup_anno$bin_start)
#head(dup_anno)
#u_inters$ID <- paste0(u_inters$chr,".",u_inters$st)
#head(u_inters)
#c_dup_anno <- dup_anno %>% filter(ID %in% u_inters$ID)
#head(c_dup_anno)
#
##sum broad category per chrom per bin
#cat_sum_bin <- c_dup_anno %>%
#  select(seqname, broad_class,bin_start) %>%
#  group_by(seqname, broad_class, bin_start) %>%
#  dplyr::summarise(n = n())
#head(cat_sum_bin)
##sum broad category per chrom
#cat_sum_chrom <- c_dup_anno %>%
#  select(seqname, broad_class,bin_start) %>%
#  group_by(seqname, broad_class) %>%
#  dplyr::summarise(n = n())
#head(cat_sum_chrom)
#
##boxplot of per chrom summary
#bdat <- cat_sum_bin
#bdat$seqname <- gsub("chr","",bdat$seqname)
#bdat$seqname <- factor(bdat$seqname, levels=p_chr_ord_gsub)
#p <- (ggplot(bdat, aes(x=seqname, y=n,fill=factor(broad_class)) )
#      + geom_boxplot()
#      + labs(title = paste(gsub("_","",Atype), "Common Interactions: Number of Genes per Chromosome"),
#             #         subtitle = "Plot of length by dose",
#             #         caption = "Data source: ToothGrowth",
#             x = paste0("Chromosome"),
#             y = "Number of Genes per 1Mb Bin",
#             fill = "Gene Category")
#      + scale_x_discrete(expand = c(0, 0))
#      + scale_y_continuous(expand = c(0, 0))
#)
#f_name <- gsub(" ","",paste("common_interactions_gene_density_per_chrom_boxplot",Atype,".pdf"))
#pdf(f_name, width = 14, height = 8)
#p
#dev.off()
#
#bdat$bin_start <- bdat$bin_start / 1000000
#p <- (ggplot(bdat, aes(x=bin_start, y=n,colour=factor(broad_class)) )
#      + geom_point()
#      + geom_smooth(method = "loess", formula = y~x)
#      + labs(title = paste(gsub("_","",Atype), "Common Interactions: Number of Genes per Chromosome"),
#             #         subtitle = "Plot of length by dose",
#             #         caption = "Data source: ToothGrowth",
#             x = paste0("Genomic Position [Mb]"),
#             y = "Number of Genes per 1Mb Bin",
#             fill = "Gene Category")
#      + scale_x_continuous(expand = c(0, 0))
#      + scale_y_continuous(expand = c(0, 0))
#      + facet_grid(seqname ~ .)
#)
#f_name <- gsub(" ","",paste("common_interactions_gene_density_per_chrom_pts_line",Atype,".pdf"))
#pdf(f_name, width = 14, height = 8)
#p
#dev.off()
#
##boxplot of all chroms together
#p <- (ggplot(bdat, aes(x=broad_class, y=n,fill=factor(broad_class)) )
#      + geom_boxplot()
#      + labs(title = paste(gsub("_","",Atype), "Common Interactions"),
#             #         subtitle = "Plot of length by dose",
#             #         caption = "Data source: ToothGrowth",
#             x = paste0("Gene Category"),
#             y = "Number of Genes per 1Mb Bin",
#             fill = "Gene Category")
#      + scale_x_discrete(expand = c(0, 0))
#      + scale_y_continuous(expand = c(0, 0))
##      + facet_grid(seqname ~ .)
#)
#f_name <- gsub(" ","",paste("common_interactions_gene_density_per_gene_boxplot",Atype,".pdf"))
#pdf(f_name, width = 14, height = 8)
#p
#dev.off()
#print("#test between means of gene classes (all chroms data)")
#print("#Test each group for normality")
#print("sig = reject normality null")
#bdat$broad_class <- as.factor(bdat$broad_class)
#bdat %>%
#  group_by(broad_class) %>%
#  summarise(W = shapiro.test(n)$statistic,
#            p.value = shapiro.test(n)$p.value)
#print("#Perform the Kruskal-Wallis test")
#print("sig = mean is diff btwn groups")
#kruskal.test(n ~ broad_class, data=bdat)
#print("# check which groups have sig diff")
#print("# perform pairwise wilcoxon test with FDR (Benjamini-Hochberg) correction")
#pairwise.wilcox.test(bdat$n, bdat$broad_class,
#                 p.adjust.method = "BH")
#
#print("lncRNA summary")
##df for non-common inters
#nc_dup_anno <- dup_anno %>% filter(!ID %in% u_inters$ID)
#nc_dup_anno$inter <- rep("nonCommon",nrow(nc_dup_anno))
##add col for common inters
#c_dup_anno$inter <- rep("common",nrow(c_dup_anno))
#comb_anno <- rbind(c_dup_anno,nc_dup_anno)
##add col for gene length
#comb_anno$len <- comb_anno$end - comb_anno$start
#lnc_df <- comb_anno %>% filter(broad_class == "lncRNA")
#head(lnc_df)
#summary(lnc_df)
#lnc_df[which(lnc_df$len == min(lnc_df$len)),]
#
##boxplot of all chroms together
#lnc_df$len <- lnc_df$len / 1000000
#p <- (ggplot(lnc_df, aes(x=inter, y=len,fill=factor(strand) ))
#      + geom_boxplot()
#      + labs(title = paste(gsub("_","",Atype), "All Interactions"),
#             #         subtitle = "Plot of length by dose",
#             #         caption = "Data source: ToothGrowth",
#             x = paste0("Interaction Category"),
#             y = "lncRNA length [Mb]",
#             fill = "Strand")
#      + scale_x_discrete(expand = c(0, 0))
#      + scale_y_continuous(expand = c(0, 0))
#      #      + facet_grid(seqname ~ .)
#)
#f_name <- gsub(" ","",paste("common_interactions_lncRNA_noncommon_boxplot",Atype,".pdf"))
#pdf(f_name, width = 14, height = 8)
#p
#dev.off()
##print("#test between means of common and non-common lncRNAs length (all chroms data)")
##print("#Test each group for normality")
##print("sig = reject normality null")
##lnc_df$inter <- as.factor(lnc_df$inter)
##lnc_df %>%
##  group_by(inter, strand) %>%
##  summarise(W = shapiro.test(len)$statistic,
##            p.value = shapiro.test(length())$p.value)
##print("#Perform ANOVA")
##res.aov2 <- aov(len ~ inter + strand, data = lnc_df)
##summary(res.aov2)
##summary(glht(res.aov2, linfct = mcp(inter = "Tukey")))
###print("sig = mean is diff btwn groups")
###kruskal.test(n ~ broad_class, data=bdat)
###print("# check which groups have sig diff")
###print("# perform pairwise wilcoxon test with FDR (Benjamini-Hochberg) correction")
###pairwise.wilcox.test(bdat$n, bdat$broad_class,
###                     p.adjust.method = "BH")
#
##boxoplot of lncRNAs in common inters ONLY per chrom
#bdat <- lnc_df
#bdat$seqname <- gsub("chr","",bdat$seqname)
#bdat$seqname <- factor(bdat$seqname, levels=p_chr_ord_gsub)
#p <- (ggplot(bdat %>% filter(inter == "common"), aes(x=seqname, y=len,fill=factor(strand) ))
#      + geom_boxplot()
#      + labs(title = paste(gsub("_","",Atype), "Common Interactions"),
#             #         subtitle = "Plot of length by dose",
#             #         caption = "Data source: ToothGrowth",
#             x = paste0("Chromosome"),
#             y = "lncRNA length [Mb]",
#             fill = "Strand")
#      + scale_x_discrete(expand = c(0, 0))
#      + scale_y_continuous(expand = c(0, 0))
#      #      + facet_grid(seqname ~ .)
#)
#f_name <- gsub(" ","",paste("common_interactions_lncRNA_common_boxplot",Atype,".pdf"))
#pdf(f_name, width = 14, height = 8)
#p
#dev.off()
print("DONE")
#
