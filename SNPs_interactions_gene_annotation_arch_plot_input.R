########################################
# all genes present in interactions
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: interactions for arch plots (in modified bed file format)
#            cells to look for SNPs in (text file with each cell name on a different row)
#            annotation file (from python script) (tsv)
#            output file prefix
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
#dat_file <- "VSMC_cis50Kb_filtered_inters_SNP_MAP4_overlapping_intearctions_for_circos.txt"
#cells_file <- "cell_subset2.txt"
##germlayer_file <- "germlayer_info.txt"
#bin_size <- 1000000
#anno_file <- "hg38_p13_v32_annotation.txt"
#outfile <- "test_cell_list"
#SNP_file <- "VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt"
#library(factoextra)#for PCA
#library(harrypotter) #for colours

#read in cells file
cells_sub <- read.table(cells_file)
head(cells_sub)

dat <- read.table(dat_file, header = FALSE) 
colnames(dat) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
dat <- dat %>% filter(cell %in% cells_sub$V1)
print("summary of ALL sig zscores per cell type")
summary(dat)
unique(dat$cell)

print("# gene anno per cell")
for (c in cells_sub$V1){
  print(c)
  tdat <- dat %>% filter(cell == c) %>% na.omit()
  print("#counting each interaction twice (once for each chrom in interaction)")
  anchD <- tdat
  anchD$AllChr <- anchD$chrA
  anchD$AllSt <- anchD$st1
  anchD$AllEnd <- anchD$end1
  tarD <- tdat
  tarD$AllChr <- tarD$chrB
  tarD$AllSt <- tarD$st2
  tarD$AllEnd <- tarD$end2
  r_dat2 <- rbind(anchD,tarD)
  #remove NA interactions
  noNAdat <- r_dat2 %>% dplyr::select(cell, zscore, AllChr, AllSt, AllEnd)
  #filter interactions for ones that contain SNPs
  noNAdat2 <- as.data.frame(noNAdat)
  noNAdat2$AllSt <- as.numeric(as.character(noNAdat2$AllSt))
  #read in annotation file
  anno_df <- read.table(anno_file, header = TRUE)
  summary(anno_df)
  head(anno_df)
  #filter annotation for common interactions
  colnames(noNAdat2) <- c("cell", "zscore","chr","st","end")
  #positive zscores
  u_inters <- distinct(noNAdat2 %>% filter(zscore >0) %>% dplyr::select(chr, st, end))
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
  print("#### number of genes in pos zscore interactions")
  print(length(common_genes))
  write.table(common_genes, file = as.character(paste0(outfile,"_",c,"_pos_zscores_inters_GO_analysis_gene_list.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(common_genes_metascape, file = as.character(paste0(outfile,"_",c,"_pos_zscores_inters_metascape_analysis_gene_list.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
 if (min(noNAdat2$zscore <0)){ 
  #negative zscores
  u_inters <- distinct(noNAdat2 %>% filter(zscore <0) %>% dplyr::select(chr, st, end))
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
  print("#### number of genes in neg zscore interactions")
  print(length(common_genes))
  write.table(common_genes, file = as.character(paste0(outfile,"_",c,"_neg_zscores_inters_GO_analysis_gene_list.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(common_genes_metascape, file = as.character(paste0(outfile,"_",c,"_neg_zscores_inters_metascape_analysis_gene_list.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}#if negative zscores  
}


##########
print("DONE")
#
