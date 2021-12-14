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
#            full path and name of output file
########################################

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
germlayer_file <- args[2]
allinters_file <- args[3]
tissue_file <- args[4]
bin_size <- args[5]
outfile <- args[6]
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
dat_file <- "test_common_inters_df.txt"
allinters_file <- "all_trans_interactions_1Mb.txt"
germlayer_file <- "germlayer_info.txt"
bin_size <- 1000000
anno_file <- "hg38_p13_v32_annotation.txt"
outfile <- "GO_analysis_common_inters_gene_list.txt"
library(factoextra)#for PCA
library(harrypotter) #for colours

#allinters <- read.table(allinters_file, header = FALSE)
#colnames(allinters) <- c("chrA", "startA", "endA", "chrB", "startB", "endB")
#head(allinters)
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
for(i in 1:nrow(u_inters)) {
  #i=1
  td <- u_inters[i,]
  tgenes_df <- anno_df %>% filter(seqname == td$chr & bin_start == td$st | bin_end == td$st) %>% 
    filter(broad_class =="prot")
  common_genes <- append(common_genes,gsub("\\..*", "",tgenes_df$gene_id, perl=TRUE))
} #for
common_genes
write.table(common_genes, file = as.character(outfile), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

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
cat_sum <- c_dup_anno %>%
  select(seqname, broad_class,bin_start) %>%
  group_by(seqname, broad_class, bin_start) %>%
  dplyr::summarise(n = n())
head(cat_sum)

