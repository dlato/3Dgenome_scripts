########################################
# all genes present in common trans-chromosomal interactions
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: common interactions output data (from R script) (tsv)
#            bin size used for Hi-C data analysis (bp)
#            annotation file (from python script) (tsv)
#            full path and name of output file
#            type of analysis with underscore between words (i.e. 1_vs_All)
#            file with all possible interactions (tsv)
########################################

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
bin_size <- args[2]
anno_file <- args[3]
outfile <- args[4]
Atype <- args[5]
allinters_file <- args[6]

##########
library(tidyr)
library(dplyr)
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
###interaction data
#Atype <- "1_vs_All"
##tissue_file <- "tissue_system_info.txt"
#dat_file <- "test_common_inters_df.txt"
#allinters_file <- "all_trans_interactions_1Mb.txt"
##germlayer_file <- "germlayer_info.txt"
#bin_size <- 1000000
#anno_file <- "hg38_p13_v32_annotation.txt"
#outfile <- "GO_analysis_common_inters_gene_list.txt"
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
#getting data for all interactions and processing
allinters <- read.table(allinters_file, header = FALSE)
colnames(allinters) <- c("chrA", "startA", "endA", "chrB", "startB", "endB")
allinters$ID <- paste0(allinters$chrA,".",allinters$startA, ".",allinters$endA,".",allinters$chrB,".",allinters$startB, ".",allinters$endB) 
head(allinters)

dat <- read.table(dat_file, header = TRUE)
print("summary of ALL sig zscores per cell type")
summary(dat)
anno_df <- read.table(anno_file, header = TRUE)
summary(anno_df)
head(anno_df)

#filter annotation for common interactions
head(dat)
Adf <- dat %>% select(chrA, st1,end1)
colnames(Adf) <- c("chr","st","end")
Bdf <- dat %>% select(chrB, st2,end1)
colnames(Bdf) <- c("chr","st","end")
u_inters <- distinct(rbind(Adf,Bdf))
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
print("#### number of genes in common interactions")
length(common_genes)
write.table(common_genes, file = as.character(paste0(outfile,"_common_inters_GO_analysis_gene_list.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(common_genes_metascape, file = as.character(paste0(outfile,"_common_inters_metascape_analysis_gene_list.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

#########
# top 3000 genes in common interactions (roughly top 50 interactions)
# calculate the mean zscore per interaction (mean across all cells), then taking the top 79 interactions
# file comes from the trans_common_sig_interactions.R script
topdat <- dat %>% arrange(desc(meanZscore))
topdat
topdat <- topdat[1:15,]
topdat
Adf <- topdat %>% select(chrA, st1,end1)
colnames(Adf) <- c("chr","st","end")
Bdf <- topdat %>% select(chrB, st2,end1)
colnames(Bdf) <- c("chr","st","end")
tu_inters <- distinct(rbind(Adf,Bdf))
summary(tu_inters)
tcommon_genes <- c()
tcommon_genes_metascape <- c()
for(i in 1:nrow(tu_inters)) {
  #i=1
  td <- tu_inters[i,]
  tgenes_df <- anno_df %>% filter(seqname == td$chr & td$st >= bin_start & td$st <= bin_end) %>%
    filter(broad_class =="prot")
  tcommon_genes <- append(tcommon_genes,gsub("\\..*", "",tgenes_df$gene_id, perl=TRUE))
  tcommon_genes_metascape <- append(tcommon_genes_metascape,gsub("\\..*", "",tgenes_df$gene_name, perl=TRUE))
} #for
tcommon_genes <- unique(tcommon_genes)
print("#### number of genes in common interactions")
length(tcommon_genes)
write.table(tcommon_genes, file = as.character(paste0(outfile,"_top3000_common_inters_GO_analysis_gene_list.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(tcommon_genes_metascape, file = as.character(paste0(outfile,"_top3000_common_inters_metascape_analysis_gene_list.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


##############
# common per chrom
##############
#add col to all possible inters identifying common/non-common ones
dat$interType <- rep("common",nrow(dat))
dat$ID <- paste0(dat$chrA,".",dat$st1, ".",dat$end1,".",dat$chrB,".",dat$st2, ".",dat$end2) 
head(dat)
allinters$interType <- dat$interType[match(allinters$ID, dat$ID)]
allinters$interType[is.na(allinters$interType)] <- "nonCommon"
head(allinters)
summary(allinters)
unique(allinters$interType)
#mark each interaction twice
chrA_df <- allinters %>% select(chrA, startA, endA, interType)
colnames(chrA_df) <- c("chrom","start","end", "interType")
chrB_df <- allinters %>% select(chrB, startB, endB, interType)
colnames(chrB_df) <- c("chrom","start","end", "interType")
allprop <- rbind(chrA_df, chrB_df)
expInter <- allprop %>%
  select(chrom, start, interType) %>%
  group_by(chrom, interType) %>%
  dplyr::summarise(n = n())
totInter <- allprop %>%
  select(chrom, start) %>%
  group_by(chrom) %>%
  dplyr::summarise(n = n())
totInter
#combine above dfs
expInter$totInter <- totInter$n[match(expInter$chrom, totInter$chrom)]
expInter$prop <- expInter$n / expInter$totInter
expInter

#plot common inter proportion per chrom
expInter$chrom <- gsub("chr","",expInter$chrom)
expInter$chrom <- factor(expInter$chrom, levels=p_chr_ord_gsub)
bp <- (ggplot(expInter %>% filter(interType == "common"), aes(x=chrom, y=prop))
       +geom_bar(position = "dodge", stat="identity", fill = "black")
         + labs(x = "Chromosome",
                y = "Proportion of Common Interactions [%]",
                title = "Proportion of Common Interactions Per Chromosome")
      + scale_x_discrete(expand = c(0, 0))
      + scale_y_continuous(expand = c(0, 0))
)
pdf("proportion_common_inters_per_chrom.pdf", width = 14, height = 8)
bp
dev.off()

print("########")
print("# test to see if proportion of common inters is higher on gene dense chroms")
print("# gene dense chroms = top 5 highest #genes/bp")
print("########")
#df with number of genes (proportional) per chrom
gene_dense_df <- anno_df %>% filter(broad_class == "prot") %>%
                 group_by(seqname, broad_class) %>%
                 dplyr::summarise(num_genes = n(), .groups = "keep")
gene_dense_df$chrlen <- chrInf$size[match(gene_dense_df$seqname, chrInf$chrom)]
gene_dense_df <- gene_dense_df %>% mutate(gene_prop = num_genes / chrlen) %>%
                 #remove chrM
                 filter(seqname != "chrM") %>%
                 #descending order for proportion
                 arrange(desc(gene_prop))
#top 5 gene dense chroms listed as "dense", others listed as non-dense
topgd <- gene_dense_df$seqname[1:5]
topgd
gene_dense_df <- gene_dense_df %>% mutate(class = ifelse(seqname %in% topgd, 'dense', 'nonDense'))
head(gene_dense_df)
summary(gene_dense_df)
#boxplot of proportion of common inters in gene dense and non-dense chroms
p <- (ggplot(gene_dense_df, aes(x=class, y=gene_prop,fill=factor(class)) )
      + geom_boxplot()
      + labs(title = paste(gsub("_","",Atype), "Common Interactions Among Gene-dense and non-dense Chromosomes"),
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome Class"),
             y = "Proportion of Common Interactions [%]",
             fill = "Gene Category")
      + scale_x_discrete(expand = c(0, 0))
      + scale_y_continuous(expand = c(0, 0))
)
f_name <- gsub(" ","",paste("common_interactions_dense_nonDense_chroms_boxplot",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
print("#Test each group for normality")
print("sig = reject normality null")
gene_dense_df %>%
  summarise(W = shapiro.test(gene_prop)$statistic,
                  p.value = shapiro.test(gene_prop)$p.value)
print("#Perform Mann-Whitney test (non-normal)")
print("sig = mean of prop of common inters is diff btwn dense and non-dense chroms")
gene_dense_df %>%
  wilcox.test(gene_prop ~ class, data = .)




# dealing with genes that are found in two bins: counting twice
tdup_anno <- anno_df %>% filter(bin_start != bin_end)
tdup_anno <- tdup_anno %>% select("seqname","source", "feature","start","end","score","strand","frame","gene_id","gene_type","gene_name","level","hgnc_id","havana_gene","transcript_id","transcript_type","transcript_name","transcript_support_level","tag","havana_transcript","exon_number","exon_id","ont","protein_id","ccdsid","broad_class","bin_end","bin_start")
colnames(tdup_anno) <- c("seqname","source", "feature","start","end","score","strand","frame","gene_id","gene_type","gene_name","level","hgnc_id","havana_gene","transcript_id","transcript_type","transcript_name","transcript_support_level","tag","havana_transcript","exon_number","exon_id","ont","protein_id","ccdsid","broad_class","bin_start","bin_end")
dup_anno <- rbind(anno_df,tdup_anno)

#filter by common interactions
dup_anno$ID <- paste0(dup_anno$seqname,".",dup_anno$bin_start)
head(dup_anno)
u_inters$ID <- paste0(u_inters$chr,".",u_inters$st)
head(u_inters)
c_dup_anno <- dup_anno %>% filter(ID %in% u_inters$ID)
head(c_dup_anno)

#sum broad category per chrom per bin
cat_sum_bin <- c_dup_anno %>%
  select(seqname, broad_class,bin_start) %>%
  group_by(seqname, broad_class, bin_start) %>%
  dplyr::summarise(n = n())
head(cat_sum_bin)
#sum broad category per chrom
cat_sum_chrom <- c_dup_anno %>%
  select(seqname, broad_class,bin_start) %>%
  group_by(seqname, broad_class) %>%
  dplyr::summarise(n = n())
head(cat_sum_chrom)

#boxplot of per chrom summary
bdat <- cat_sum_bin
bdat$seqname <- gsub("chr","",bdat$seqname)
bdat$seqname <- factor(bdat$seqname, levels=p_chr_ord_gsub)
p <- (ggplot(bdat, aes(x=seqname, y=n,fill=factor(broad_class)) )
      + geom_boxplot()
      + labs(title = paste(gsub("_","",Atype), "Common Interactions: Number of Genes per Chromosome"),
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome"),
             y = "Number of Genes per 1Mb Bin",
             fill = "Gene Category")
      + scale_x_discrete(expand = c(0, 0))
      + scale_y_continuous(expand = c(0, 0))
)
f_name <- gsub(" ","",paste("common_interactions_gene_density_per_chrom_boxplot",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()

bdat$bin_start <- bdat$bin_start / 1000000
#facet line plot
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

#boxplot of all chroms together
p <- (ggplot(bdat, aes(x=broad_class, y=n,fill=factor(broad_class)) )
      + geom_boxplot()
      + labs(title = paste(gsub("_","",Atype), "Common Interactions"),
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Gene Category"),
             y = "Number of Genes per 1Mb Bin",
             fill = "Gene Category")
      + scale_x_discrete(expand = c(0, 0))
      + scale_y_continuous(expand = c(0, 0))
#      + facet_grid(seqname ~ .)
)
f_name <- gsub(" ","",paste("common_interactions_gene_density_per_gene_boxplot",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
print("#test between means of gene classes (all chroms data)")
print("#Test each group for normality")
print("sig = reject normality null")
bdat$broad_class <- as.factor(bdat$broad_class)
bdat %>%
  group_by(broad_class) %>%
  summarise(W = shapiro.test(n)$statistic,
            p.value = shapiro.test(n)$p.value)
print("#Perform the Kruskal-Wallis test")
print("sig = mean is diff btwn groups")
kruskal.test(n ~ broad_class, data=bdat)
print("# check which groups have sig diff")
print("# perform pairwise wilcoxon test with FDR (Benjamini-Hochberg) correction")
pairwise.wilcox.test(bdat$n, bdat$broad_class,
                 p.adjust.method = "BH")

#df for non-common inters
nc_dup_anno <- dup_anno %>% filter(!ID %in% u_inters$ID)
nc_dup_anno$inter <- rep("nonCommon",nrow(nc_dup_anno))
#add col for common inters
c_dup_anno$inter <- rep("common",nrow(c_dup_anno))
comb_anno <- rbind(c_dup_anno,nc_dup_anno)
#add col for gene length
comb_anno$len <- comb_anno$end - comb_anno$start
head(comb_anno)

############
# number of gene category per interaction
############
cat_num_interType <- comb_anno %>%
  select(ID, broad_class, inter) %>%
  group_by(ID, broad_class, inter) %>%
  dplyr::summarise(n = n())
  #dplyr::summarise(n = n(), .groups = "keep", .data = cat_num_interType)
tail(cat_num_interType)
summary(cat_num_interType)


print("#test between means of common and non-common inters in each gene class")
cat_num_interType$broad_class <- as.factor(cat_num_interType$broad_class)
for (i in unique(cat_num_interType$broad_class)){
  print("########")
  print(i)
  #i="lncRNA"
  print("########")
  print("#Test each group for normality")
  print("sig = reject normality null")
  print(cat_num_interType %>%
    filter(broad_class == i) %>%
    summarise(W = shapiro.test(n)$statistic,
              p.value = shapiro.test(n)$p.value))
  td <- cat_num_interType %>% filter(broad_class==i)
  if (length(unique(td$inter) == 2)){
  print("#Perform Mann-Whitney test (non-normal)")
  print("sig = mean is diff btwn common and non-common are diff")
  print(cat_num_interType %>%
    filter(broad_class == i) %>%
    wilcox.test(n ~ inter, data = .))
} else {
print("cannot perform Mann-Whitney test because improper grouping variables")
}#ifelse
}#for
#plot box plot
p <- (ggplot(cat_num_interType, aes(x=broad_class, y=n,fill=factor(inter)) )
      + geom_boxplot()
      + labs(title = paste(gsub("_","",Atype), "Common Interactions"),
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Gene Category"),
             y = "Number of Genes per 1Mb Bin",
             fill = "Interaction Type")
      + scale_x_discrete(expand = c(0, 0))
      + scale_y_continuous(expand = c(0, 0))
      #      + facet_grid(seqname ~ .)
)
f_name <- gsub(" ","",paste("common_nonCommon_interactions_gene_density_per_gene_class_boxplot",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()

print("lncRNA summary")
lnc_df <- comb_anno %>% filter(broad_class == "lncRNA")
head(lnc_df)
summary(lnc_df)
lnc_df[which(lnc_df$len == min(lnc_df$len)),]

print("#perform the Mann Whitney U test to determine diff betwen lncRNA len for common and non-common inters")
print("sig = diff in lncRNA len between common and non-common inters")
wilcox.test(len~inter, data = lnc_df)

print("#perform the Mann Whitney U test to determine diff betwen lncRNA len and strand within common inters")
print("sig = diff in lncRNA len between strand in common inters")
c_lnc_df <- lnc_df %>% filter(inter == "common")
wilcox.test(len~strand, data = c_lnc_df)
print("#perform the Mann Whitney U test to determine diff betwen lncRNA len and strand within non-common inters")
print("sig = diff in lncRNA len between strand in non-common inters")
nc_lnc_df <- lnc_df %>% filter(inter == "nonCommon")
wilcox.test(len~strand, data = nc_lnc_df)
#boxplot of all chroms together
lnc_df$len <- lnc_df$len / 1000000
bpdf_log <- lnc_df
bpdf_log$len <- log(bpdf_log$len,10)
p <- (ggplot(bpdf_log, aes(x=inter, y=len,fill=factor(strand) ))
      + geom_boxplot()
      + labs(title = paste(gsub("_","",Atype), "All Interactions"),
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Interaction Category"),
             y = "lncRNA length log10[Mb]",
             fill = "Strand")
      + scale_x_discrete(expand = c(0, 0))
      + scale_y_continuous(expand = c(0, 0))
      #      + facet_grid(seqname ~ .)
)
f_name <- gsub(" ","",paste("common_interactions_lncRNA_noncommon_boxplot",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
#print("#test between means of common and non-common lncRNAs length (all chroms data)")
#print("#Test each group for normality")
#print("sig = reject normality null")
#lnc_df$inter <- as.factor(lnc_df$inter)
#lnc_df %>%
#  group_by(inter, strand) %>%
#  summarise(W = shapiro.test(len)$statistic,
#            p.value = shapiro.test(length())$p.value)
#print("#Perform ANOVA")
#res.aov2 <- aov(len ~ inter + strand, data = lnc_df)
#summary(res.aov2)
#summary(glht(res.aov2, linfct = mcp(inter = "Tukey")))
##print("sig = mean is diff btwn groups")
##kruskal.test(n ~ broad_class, data=bdat)
##print("# check which groups have sig diff")
##print("# perform pairwise wilcoxon test with FDR (Benjamini-Hochberg) correction")
##pairwise.wilcox.test(bdat$n, bdat$broad_class,
##                     p.adjust.method = "BH")

#boxoplot of lncRNAs in common inters ONLY per chrom
bdat <- lnc_df
bdat$seqname <- gsub("chr","",bdat$seqname)
bdat$seqname <- factor(bdat$seqname, levels=p_chr_ord_gsub)
p <- (ggplot(bdat %>% filter(inter == "common"), aes(x=seqname, y=len,fill=factor(strand) ))
      + geom_boxplot()
      + labs(title = paste(gsub("_","",Atype), "Common Interactions"),
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome"),
             y = "lncRNA length [Mb]",
             fill = "Strand")
      + scale_x_discrete(expand = c(0, 0))
      + scale_y_continuous(expand = c(0, 0))
      #      + facet_grid(seqname ~ .)
)
f_name <- gsub(" ","",paste("common_interactions_lncRNA_common_boxplot",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
print("DONE")
