########################################
# determine overlap of TARGET interactions
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: CM interactions for arch plots (in modified bed file format)
#            VSCM interactions for arch plots (in modified bed file format)
#            CM cells to look for overlap in (text file with each cell name on a different row)
#            VSMC cells to look for overlap in (text file with each cell name on a different row)
#            CM SNP file (tsv with SNP location and info)
#            VSMC SNP file (tsv with SNP location and info)
#            bin size (numeric, bp)
#            output file prefix
########################################


options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
CMSNPs_intersfile <- args[1]
VSMCSNPs_intersfile <- args[2]
CM_cells_file <- args[3]
VSMC_cells_file <- args[4]
CMSNPs_file <- args[5]
VSMCSNPs_file <- args[6]
bin_size <- as.numeric(as.character(args[7]))
outprefix <- args[8]
##########
library(tidyr)
library(purrr)
library(forcats)
library(dplyr)
library(ggplot2)
#library(nortest, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.3.6.1") #for normality test with large sample size
#library(hexbin)
library(ggforce)#for ridgeline
library(ggridges)#for ridgeline
.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.3.6.1")
library(ggupset) #for UpSet plot
library(ggVennDiagram)#for venn diagram
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
#CMSNPs_intersfile <- "CardioMyo_cis50Kb_filtered_inters_SNP_CFDP1_overlapping_intearctions_for_circos.txt"
#VSMCSNPs_intersfile <- "VSMC_cis50Kb_filtered_inters_SNP_MAP4_overlapping_intearctions_for_circos.txt"
#CM_cells_file <- "cell_subset_CM.txt"
#VSMC_cells_file <- "cell_subset2.txt"
#CMSNPs_file <- "CM_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt"
#VSMCSNPs_file <- "VSMC_diff_snps_final.ranking.withinfo.eqtl.Repeat.txt"
#bin_size <- 50000
#library(harrypotter)
#library(factoextra)
#library(hexbin)
#library(ggVennDiagram)
#library(ggupset)

print("#read in CM SNP inters")
CMSNP_inters <- read.table(CMSNPs_intersfile, header = F)
colnames(CMSNP_inters) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
CMcells <- read.table(CM_cells_file)
#filter for cells
CMSNP_inters <- CMSNP_inters %>% filter(cell %in% CMcells$V1)
print("#counting each interaction twice (once for each chrom in interaction)")
anchD <- CMSNP_inters
anchD$AllChr <- anchD$chrA
anchD$AllSt <- anchD$st1
anchD$AllEnd <- anchD$end1
tarD <- CMSNP_inters
tarD$AllChr <- tarD$chrB
tarD$AllSt <- tarD$st2
tarD$AllEnd <- tarD$end2
r_dat2 <- rbind(anchD,tarD)
summary(r_dat2)
#remove NA interactions
noNAdat <- r_dat2 %>% dplyr::select(cell, zscore, AllChr, AllSt, AllEnd)
#filter interactions for ones that contain SNPs
noNAdat2 <- as.data.frame(noNAdat)
noNAdat2$AllSt <- as.numeric(as.character(noNAdat2$AllSt))
#filter annotation for common interactions
colnames(noNAdat2) <- c("cell", "zscore","chr","st","end")
#ID col for easier filtering
CMSNP_inters_bin <- noNAdat2 %>% mutate(ID = paste0(chr,".",st, ".",end))
head(CMSNP_inters_bin)

print("#read in VSMC SNP inters")
VSMCSNPs_intersfile
VSMCSNP_inters <- read.table(VSMCSNPs_intersfile, header = FALSE,sep="\t")
colnames(VSMCSNP_inters) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
print("read in cells")
VSMCcells <- read.table(VSMC_cells_file, sep="\t")
#filter for cells
VSMCSNP_inters <- VSMCSNP_inters %>% filter(cell %in% VSMCcells$V1)
print("#counting each interaction twice (once for each chrom in interaction)")
anchD <- VSMCSNP_inters
anchD$AllChr <- anchD$chrA
anchD$AllSt <- anchD$st1
anchD$AllEnd <- anchD$end1
tarD <- VSMCSNP_inters
tarD$AllChr <- tarD$chrB
tarD$AllSt <- tarD$st2
tarD$AllEnd <- tarD$end2
r_dat2 <- rbind(anchD,tarD)
summary(r_dat2)
#remove NA interactions
noNAdat <- r_dat2 %>% dplyr::select(cell, zscore, AllChr, AllSt, AllEnd)
#filter interactions for ones that contain SNPs
noNAdat2 <- as.data.frame(noNAdat)
noNAdat2$AllSt <- as.numeric(as.character(noNAdat2$AllSt))
#filter annotation for common interactions
colnames(noNAdat2) <- c("cell", "zscore","chr","st","end")
#ID col for easier filtering
VSMCSNP_inters_bin <- noNAdat2 %>% mutate(ID = paste0(chr,".",st, ".",end))
head(VSMCSNP_inters_bin)

print("# read in CM SNP file")
CMSNP_df <- read.table(CMSNPs_file, header = T, sep = "\t")
#accounting for merged SNP files
if ("logFC_comp" %in% colnames(CMSNP_df)){
  CMSNP_df <- CMSNP_df %>% dplyr::select(chr_b38, start_b38,end_b38, logFC_comp)
} else {
  CMSNP_df$logFC_comp.x <- as.numeric(as.character(CMSNP_df$logFC_comp.x))
  CMSNP_df$logFC_comp.y <- as.numeric(as.character(CMSNP_df$logFC_comp.y))
  CMSNP_df$logFC_comp <- mean(c(abs(CMSNP_df$logFC_comp.x),abs(CMSNP_df$logFC_comp.y)))
}#if else
print("#bins that SNPs are in")
CMSNP_df$bin <- plyr::round_any(CMSNP_df$start_b38, as.numeric(as.character(bin_size)), f = floor)
colnames(CMSNP_df) <- c("AllChr","SNPstart","SNPend","logFC","AllSt")
CMSNP_df$ID <- format(CMSNP_df$ID, scientific = FALSE)
CMSNP_df$ID <- paste0(CMSNP_df$AllChr,".",format(CMSNP_df$AllSt, scientific = FALSE),".",format(CMSNP_df$AllSt + as.numeric(as.character(bin_size)),scientific = FALSE))
#unique chroms containing SNPs
CMSNP_df$ID <- gsub(" ","", CMSNP_df$ID)
#decrease SNP list to just unique IDs, find mean of |logFC|
CMSNP_uniq <- CMSNP_df %>%
  group_by(ID) %>%
  dplyr::summarise(mlogFC = mean(abs(logFC)))
head(CMSNP_uniq)

print("# read in VSMC SNP file")
VSMCSNPs_file
VSMCSNP_df <- read.table(VSMCSNPs_file, header = TRUE, sep = "\t")
#accounting for merged SNP files
if ("logFC_comp" %in% colnames(VSMCSNP_df)){
  VSMCSNP_df <- VSMCSNP_df %>% dplyr::select(chr_b38, start_b38,end_b38, logFC_comp)
} else {
  VSMCSNP_df$logFC_comp.x <- as.numeric(as.character(VSMCSNP_df$logFC_comp.x))
  VSMCSNP_df$logFC_comp.y <- as.numeric(as.character(VSMCSNP_df$logFC_comp.y))
  VSMCSNP_df$logFC_comp <- mean(c(abs(VSMCSNP_df$logFC_comp.x),abs(VSMCSNP_df$logFC_comp.y)))
}#if else
print("#bins that SNPs are in")
VSMCSNP_df$bin <- plyr::round_any(VSMCSNP_df$start_b38, as.numeric(as.character(bin_size)), f = floor)
colnames(VSMCSNP_df) <- c("AllChr","SNPstart","SNPend","logFC","AllSt")
VSMCSNP_df$ID <- format(VSMCSNP_df$ID, scientific = FALSE)
VSMCSNP_df$ID <- paste0(VSMCSNP_df$AllChr,".",format(VSMCSNP_df$AllSt, scientific = FALSE),".",format(VSMCSNP_df$AllSt + as.numeric(as.character(bin_size)),scientific = FALSE))
#unique chroms containing SNPs
VSMCSNP_df$ID <- gsub(" ","", VSMCSNP_df$ID)
#decrease SNP list to just unique IDs, find mean of |logFC|
VSMCSNP_uniq <- VSMCSNP_df %>%
  group_by(ID) %>%
  dplyr::summarise(mlogFC = mean(abs(logFC)))
head(VSMCSNP_uniq)

#filter for just target inters
CM_target <- CMSNP_inters_bin %>% filter(!ID %in% CMSNP_uniq$ID)
head(CM_target)
VSMC_target <- VSMCSNP_inters_bin %>% filter(!ID %in% VSMCSNP_uniq$ID)
head(VSMC_target)

#format table for UpSet plot
merged_target <- rbind(CM_target,VSMC_target) %>% mutate(cell=factor(cell))
head(merged_target)
allcells <- c(CMcells$V1,VSMCcells$V1)
#allcells <- c("H9hESC_day00_Zhang","Ventricular_cardiomyocyte_day80_Zhang","OmniC_pooled_V_21d")
allcells
inter_list <- vector("list", length(allcells))
for (i in 1:length(allcells)) {
  tcell <- allcells[i]
  tdat <- merged_target %>% filter(cell == tcell) %>% dplyr::select(ID)
  inter_list[[i]] =  tdat$ID
}
inter_list
#plot venn diagram
p1 <- (ggVennDiagram(inter_list,
                    category.names = allcells,
                    label_alpha = 0)
       + scale_fill_distiller(palette = "Reds", direction = 1)
       + scale_color_manual(values = rep("black", length(allcells)))
       + labs(title = "Cis-chromosomal interactions overlapping with SNPs (significant)",
              subtitle = paste("N=",length(unique(merged_target$ID))),
              fill = "# of Interactions")
       #below so you can read the cell names
       + scale_x_continuous(expand = expansion(mult = .2))
)
pdf(paste0(outprefix,"_vennDiagram_TARGET.pdf"), width = 14, height = 4)
p1
dev.off()


###########
# UpSet plot (instead of venn diagram)
###########
#make tf dataframe with interaction and cell info
tf_dat <- merged_target %>%
  group_by(ID,cell) %>%
  summarise(mZscore = mean(zscore), .groups = "keep") %>%
  mutate(tf = (!is.na(mZscore))) %>%
  dplyr::select(ID,cell, tf) %>%
  spread(ID, tf)%>%
  replace(is.na(.), FALSE)
tidy_zdat <- tf_dat %>%
  as_tibble(rownames = "cell") %>%
  gather(ID, Member, -cell) %>%
  filter(Member) %>%
  select(- Member) %>%
  group_by(ID) %>%
  summarize(cell = list(cell))
head(tidy_zdat)
up <- (ggplot(tidy_zdat, aes(x=cell))
       + geom_bar()
       #+ geom_text(stat='count', aes(label=after_stat(count)), vjust=-1,hjust=-1,size=3,angle = 90)
       + geom_text(stat='count', aes(label=after_stat(count)),hjust=-0.2,size=3,angle = 90)
       #+ scale_x_upset(n_intersections = 90)
       + scale_x_upset(order_by = "freq")
       + scale_y_continuous(breaks = NULL, name = "", lim = c(0,150))
       + labs(x = "",
              y = "",
              title = "# of Cis-chromosomal interaction targets")
)
pdf(paste0(outprefix,"_TARGET_UpSet_plot.pdf"), width = 14, height = 4)
up
dev.off()
##########
print("DONE")
#
