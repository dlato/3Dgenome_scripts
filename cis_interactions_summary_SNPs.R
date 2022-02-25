########################################
#validating known trans-chromosomal interactions
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: filtered interactions file (tab separated, first 6 columns are a bed-like format, last few columns are other info)
#            output file prefix (character)
#            cell order for graphs (file with name of cells as it matches the data frame, one cell per line)
########################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
zdat_file <- args[1]
outprefix <- args[2]
cellsfile <- args[3]

##########
library(tidyr)
library(dplyr)
library(ggplot2)
library(harrypotter, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.3.6.1")
#library(nortest, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.3.6.1") #for normality test with large sample size
#library(hexbin)
library(ggforce)#for ridgeline
library(ggridges)#for ridgeline
.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.3.6.1")
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
#zdat_file <- "cis_arch_plot_test_dat.txt"
#outprefix <- "test_cis_SNP_blood_pressure"
#cellsfile <- "cell_subset.txt"
#library(harrypotter)
#library(factoextra)
#library(hexbin)
#library(ggVennDiagram)

#chrom info
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
chrInf$chrom <- factor(chrInf$chrom, levels=p_chr_ord)

#read in cells order
cells <- read.table(cellsfile, header=FALSE)

zdat <- read.table(zdat_file, header = FALSE)
colnames(zdat) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
head(zdat)
#new column with distance between starts
zdat$dist <- abs(zdat$st1 - zdat$st2)
head(zdat)
zdat <- zdat %>% filter(cell %in% cells$V1) %>% mutate(cell=factor(cell, levels=cells$V1))
summary(zdat)

#########
# overlap of interactions
#########
#create ID column
zdat <- zdat %>% mutate(ID =paste0(chrA,".",st1,".",end1, ".",chrB,".",st2,".","end2"))
summary(zdat)
inter_list <- vector("list", length(cells$V1))
for (i in 1:length(cells$V1)) {
  tcell <- cells$V1[i]
  tdat <- zdat %>% filter(cell == tcell) %>% dplyr::select(ID)
  inter_list[[i]] =  tdat$ID
}
inter_list
#plot venn diagram
p1 <- (ggVennDiagram(inter_list,
                    category.names = cells$V1,
                    label_alpha = 0)
       + scale_fill_distiller(palette = "Reds", direction = 1)
       + scale_color_manual(values = rep("black", length(cells$V1)))
       + labs(title = "Cis-chromosomal interactions overlapping with SNPs (significant)",
              subtitle = paste("N=",length(unique(zdat$ID))),
              fill = "# of Interactions")
       #below so you can read the cell names
       + scale_x_continuous(expand = expansion(mult = .2))
)
pdf(paste0(outprefix,"_vennDiagram_sig_Interactions.pdf"), width = 14, height = 4)
p1
dev.off()

print("# overlapping interactions between select cells")
#long to wide format
wdf <- zdat %>% spread(key=cell, value=zscore)
#if we are dealing with VSMC SNPs
if ( "OmniC_pooled_M_0d" %in% cells$V1){
  #remove NAs, therefore only inters that are present between these two cells
  tdat <- wdf %>% na.omit()
  write.table(tdat %>% select(-ID, -dist), file = as.character(paste0(outprefix,"_VSMC_SNPs_overlapping_interactions.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
} else { # CM SNPs
  #find overlap only between day 80 and day 0 cells
  tdat <- wdf %>% dplyr::select(chrA, st1, end1, chrB, st2, end2, H9hESC_day00_Zhang, Ventricular_cardiomyocyte_day80_Zhang) %>% na.omit()
  write.table(tdat, file = as.character(paste0(outprefix,"_CM_SNPs_overlapping_interactions.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

print("# cell type specific interactions")
for (i in cells$V1){
  #remove current cell from cells list
#  i = "Lung"
  tcells <- cells$V1[! cells$V1 %in% i]
  tdat <- wdf
  tdat
  # NAs in all other cells
  for (x in tcells){
  #  x= "Bladder_rep1"
    tdat <- tdat %>% filter(is.na(!!sym(x)))
  }
  write.table(tdat %>% dplyr::select(-ID, -dist), file = as.character(paste0(outprefix,"_",i,"_cell_specific_interactions.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

##############

print("# number of sig inters per bin")
#count  each inter twice
tp_A <- zdat %>% dplyr::select(chrA,st1,end1,cell,zscore) 
colnames(tp_A) <- c("chr","st","end","cell","zscore")
tp_B <- zdat %>% dplyr::select(chrB,st2,end2,cell,zscore) 
colnames(tp_B) <- c("chr","st","end","cell","zscore")
tp_dat <- rbind(tp_A,tp_B)
tp_dat_sum = tp_dat %>% group_by(chr,st,cell) %>% dplyr::summarize(numInters=n(),.groups = "keep")
tp_dat_sum$chr <- gsub("chr", "", tp_dat_sum$chr)
tp_dat_sum
#tp_dat_sum <- tp_dat_sum %>% mutate(chr=factor(chr, levels=p_chr_ord))
#scale pts by 1Mb
#tp_dat_sum$st <- tp_dat_sum$st / 1000000
summary(tp_dat_sum)
head(tp_dat_sum)
for(i in unique(tp_dat_sum$chr)) {
  #i="10"
  print(i)
  xmax <- chrInf$size[match(paste0("chr",i),chrInf$chrom)] / 1000000
  tpp <- tp_dat_sum %>% filter(chr == i)
  tp <- (ggplot(tpp, aes(st/1000000, y=1, fill = numInters))
         + geom_tile(aes(fill = numInters), width = 1, height = 1)
         + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Total number of significant interactions", na.value = "grey")
         #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
         #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
         + labs(x = paste0("Chromosome ", i, " position [Mb]"),
                y = "",
                title = "Cis-chromosomal interactions overlapping with SNPs (significant)")
         + facet_grid(cell ~ .)
         #       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
         #         c("nonsig" = "Non-significant", "sig" = "Significant"))))
         #       + theme(axis.text.x = element_text(angle = 90))
         + expand_limits(x = c(0,xmax))
         + scale_y_continuous(expand = c(0, 0))
         + scale_x_continuous(expand = c(0, 0))
         + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
                 strip.background = element_rect(fill = "white"),
                 panel.spacing = unit(0, "lines"),
                 axis.text.y = element_blank(),
                 axis.ticks.y = element_blank())
         + theme(panel.background = element_rect(fill = "grey85", colour = NA))
  )
  filename <- paste0(outprefix,"_numSig_inters_chrom",i,"_tickplot_sig_Interactions.pdf")
  pdf(filename, width = 14, height = 4)
  print(tp)
  dev.off()
#  # save df from bubble graph for quantification test
#  write.table(tpp, file = as.character(paste0("chrom",i,"_bubble_data.txt")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
#  # line plot of mean zscore of all cell types per bin
#  lp = tpp %>% group_by(chr,st) %>% dplyr::summarize(mnumSig=mean(numSig, na.rm = TRUE))
#  tp <- (ggplot(lp, aes(st, y= mnumSig))
#         + geom_point(alpha = 0.5, size= 3)
#         +geom_smooth(method = 'loess',formula ='y ~ x')
#         #+ scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score per bin", na.value = "grey")
#         #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
#         #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
#         + labs(x = paste0("Chromosome ", i, " position [Mb]"),
#                y = "Mean number of significant interactions per bin across cells",
#                title = "Trans-chromosomal significant interactions (all cells)")
#         #+ facet_grid(cell ~ .)
#         + geom_vline(xintercept = interpos, colour = "red")
#         #       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
#         #         c("nonsig" = "Non-significant", "sig" = "Significant"))))
#         #       + theme(axis.text.x = element_text(angle = 90))
#         + expand_limits(x = c(0,xmax))
#         + scale_y_continuous(expand = c(0, 0))
#         + scale_x_continuous(expand = c(0, 0))
#         + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
#                 strip.background = element_rect(fill = "white"),
#                 panel.spacing = unit(0, "lines"))
#         #axis.text.y = element_blank(),
#         #axis.ticks.y = element_blank())
#         #+ theme(panel.background = element_rect(fill = "grey85", colour = NA))
#  )
#  filename <- paste0("numSig_inters_",i,"_mean_line_sig_Interactions.pdf")
#  pdf(filename, width = 14, height = 8)
#  print(tp)
#  dev.off()
}#for

######
# number of inters per cells
######
ncell_df = zdat %>% group_by(cell) %>% dplyr::summarize(numInters=n(),.groups = "keep")
ncell_df
summary(ncell_df)
# bar plot of num inters per cell
p <- (ggplot(ncell_df, aes(y=cell, x=numInters))
      + geom_bar(stat = "identity")#stack = based on counts of data, height proportional to total
      + labs(title = "",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = "Number of Interactions", y = "")
        #         tag = "A")
      + scale_x_continuous(expand = c(0, 0))
      + scale_y_discrete(expand = c(0, 0))
)
f_name <- gsub(" ","",paste(outprefix,"_bar_plot_num_inters.pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
###############
# how far apart are interacting regions
##############


print("# histogram of how far apart interacting regions are")
p <- (ggplot(zdat, aes(x=dist/1000000))
      + geom_histogram(position="dodge",alpha=.6)#stack = based on counts of data, height proportional to total
      + labs(title = "Distance between interacting regions",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = "Distance [Mb]", y = "Number of Interactions")
      #         tag = "A")
      + facet_grid(cell~. )
     + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
             strip.background = element_rect(fill = "white"),
             panel.spacing = unit(0, "lines"),
             axis.text.y = element_blank(),
             axis.ticks.y = element_blank()) 
     + scale_y_continuous(expand = c(0, 0))
     + scale_x_continuous(expand = c(0, 0))
)
f_name <- gsub(" ","",paste(outprefix,"_histogram_facet_distance_between_interactions.pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()

print("# boxplot/violin of distance between interactions per cell")
p <- (ggplot(zdat, aes(x=cell, y=dist /1000000) )
      + geom_violin(alpha = 0.6, fill = "grey")
      + geom_boxplot(alpha =0.6, fill = "grey")
      + coord_flip()
      + labs(title = "Distance between interacting regions",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = "",
             y = "Distance [Mb]",
             fill = "Cell")
      + scale_x_discrete(expand = c(0, 0))
      + scale_y_continuous(expand = c(0, 0))
)
f_name <- gsub(" ","",paste(outprefix,"_boxplot_violin_distance_between_interactions.pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()

print("# ridgeline plot for distance between interactions")
p <- (ggplot(zdat, aes(x = dist/1000000, y = cell))
      + stat_density_ridges(quantile_lines = TRUE, alpha = 0.3, scale=2, quantiles = 2, rel_min_height = 0.001)
      #+ geom_density_ridges(scale = 4, alpha = 0.3) 
      + labs(x="Distance [Mb]",
             y="Cell",
             title = "Distance between interactions")
     # + scale_fill_manual(values = gl_colours)
      #                                     gsub("F", "A", my_colors)))
      + scale_y_discrete(expand = c(0, 0))     # will generally have to set the `expand` option
      + scale_x_continuous(expand = c(0, 0))   # for both axes to remove unneeded padding
      + coord_cartesian(clip = "off") # to avoid clipping of the very top of the top ridgeline
      #+ facet_grid(.~ sig, labeller = labeller(sig= as_labeller(
      #  c("nonsig" = "Non-Significant", "sig" = "Significant"))))
)
pdf(paste0(outprefix,"_ridgeline_distance_between_interactions.pdf"), width = 14, height = 8)
p
dev.off()

print("#test between means of distance between cis interactions between cell types")
#print("#Test each group for normality")
#print("sig = reject normality null")
#for (i in unique(zdat$cell)){
#  print(i)
#  subdf <- zdat %>% filter(cell == i)
#  # apply the rest of your analysis there using subdf, for instance 
#  print(shapiro.test(subdf$dist))
#}
print("#Perform the Kruskal-Wallis rank sum test (when dists are not normal)")
print("sig = mean is diff btwn groups")
kruskal.test(dist ~ cell, data = zdat)
print("#Perform multiple pairwise tests to see which groups are different (when dists are not normal)")
print("sig = mean is diff btwn groups")
pairwise.wilcox.test(zdat$dist, zdat$cell,
                     p.adjust.method = "BH")



print("DONE")
