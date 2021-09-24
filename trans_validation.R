########################################
#validating known trans-chromosomal interactions
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: 3Dflow zscore output data (ALL INTERACTIONS)
#            first interacting region bed file (tab separated)
#            second interacting region bed file (tab separated)
#            type of analysis/plot title. (words must be separated by underscore "_")
#            3Dflow pvalue output data (ALL INTERACTIONS)
# NOTE: at the moment only ONE interacting region can be searched at a time
#       i.e. chrA:1-10 interacting with chrB:40-50
########################################

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
zdat_file <- args[1]
roi1_file <- args[2]
roi2_file <- args[3]
Atype <- args[4]
pdat_file <- args[5]

##########
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(ggplot2)
library(harrypotter, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.3.6.1")
library(nortest, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.3.6.1") #for normality test with large sample size
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
#interaction data
#Atype <- "1_vs_All"
#zdat_file <- "test_1vsAll_dat.txt"
#pdat_file <- "test_1vsAll_pvalues.txt"
#roi1_file <- "FIRRE.bed"
#roi2_file <- "ATF4.bed"
#dat <- read.table("23Jul21.primary.trans.1MB.zscores.txt", header = TRUE)
#dat <- read.table("23Jul21.primary.trans.1MB.zscores.pairwise.txt", header = TRUE)
#dat <- read.table(dat_file, header = TRUE)
zdat <- read.table(zdat_file, header = TRUE)
pdat <- read.table(pdat_file, header = TRUE)
print("summary of ALL zscores per cell type")
summary(zdat)
print("summary of ALL pvalues per cell type")
summary(pdat)
#combine zscore and pvalue into one df
zdatL <- gather(zdat, key = "cell", value = "zscore", 2:length(colnames(zdat)))
pdatL <- gather(pdat, key = "cell", value = "pvalue", 2:length(colnames(pdat)))
dat <- merge(zdatL, pdatL, by=c("ID","cell"))


datsigW <- dat %>% filter(pvalue <=0.05) %>% select(-pvalue) %>% spread(key = cell, value = zscore)
print("summary of ALL sig zscores per cell type")
summary(datsigW)
#roi1
#roi1 <- read.table("CISTR.bed", header = FALSE)
roi1 <- read.table(roi1_file)
colnames(roi1) <- c("chrom","start","end")
#roi2
#roi2 <- read.table("SOX9.bed", header = FALSE)
roi2 <- read.table(roi2_file)
colnames(roi2) <- c("chrom","start","end")

print("#new cols for interaction ID")
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
dat$ID <- sub("B", "\\.B", as.character(dat$ID))
dat2 <- dat %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
#remove A and B from chrom names
dat2$chrA <- gsub("A", "", dat2$chrA)
dat2$chrB <- gsub("B", "", dat2$chrB)
hm_zscore_df <- dat2
#convert data to GRanges object
dat_bed1 <- dat2 %>% select(chrA, st1, end1)
dat_bed1$st1 <- as.numeric(dat_bed1$st1)
dat_bed1$end1 <- as.numeric(dat_bed1$end1)
dat_chr1 <- makeGRangesFromDataFrame(dat_bed1, seqnames.field=c("chrA"),
                                     start.field="st1",
                                     end.field=c("end1"))
dat_bed2 <- dat2 %>% select(chrA, st2, end2)
dat_bed2$st2 <- as.numeric(dat_bed2$st2)
dat_bed2$end2 <- as.numeric(dat_bed2$end2)
dat_chr2 <- makeGRangesFromDataFrame(dat_bed2, seqnames.field=c("chrA"),
                                     start.field="st2",
                                     end.field=c("end2"))
roi1_G <- makeGRangesFromDataFrame(roi1)
roi2_G <- makeGRangesFromDataFrame(roi2)

print("#find interactions that are within target interacting regions")
print("#resolution with roi1")
roi1_inter1 <- unique(subsetByOverlaps(dat_chr1, roi1_G))
roi1_inter2 <- unique(subsetByOverlaps(dat_chr2, roi1_G))
roi1_inter1$overl<- countOverlaps(roi1_inter1, roi1_inter2, type="equal")
roi1_inter1 <- as.data.frame(roi1_inter1)
roi1_res <- paste0(roi1_inter1$seqnames,".",roi1_inter1$start,".", roi1_inter1$end)
print("#resolution with roi2")
roi2_inter1 <- unique(subsetByOverlaps(dat_chr1, roi2_G))
roi2_inter2 <- unique(subsetByOverlaps(dat_chr2, roi2_G))
roi2_inter1$overl<- countOverlaps(roi2_inter1, roi2_inter2, type="equal")
roi2_inter1 <- as.data.frame(roi2_inter1)
roi2_res <- paste0(roi2_inter1$seqnames,".",roi2_inter1$start,".", roi2_inter1$end)
print("#data with JUST interacting region")
Vinter <- dat2[grep(roi1_res, dat2$ID), ]
Vinter <- Vinter[grep(roi2_res, Vinter$ID),]
chrs <- c(levels(roi2$chrom),levels(roi1$chrom[1]))
print("#######")
print("# Validated interaction specific region")
print("#######")
chrs_df <- Vinter[grep(chrs[1], Vinter$ID), ]
chrs_df <- chrs_df[grep(chrs[2],chrs_df$ID),]
chrs_df_allInters <- chrs_df
chrs_df <- chrs_df %>% select(-pvalue) %>% spread(key = cell, value = zscore)
summary(chrs_df)
print("cells/datasets that have interactions btwn validated chromosomes")
allmisscols <- sapply(chrs_df, function(x) all(is.na(x) | x == '' ))
all_cells <- colnames(chrs_df)
sig_cells <- all_cells[which(allmisscols == FALSE)]
sig_cells <- sig_cells[-c(1,2,3,4,5,6,7)]
sig_cells
print("cells/datasets that have interactions btwn validated chromosomes")
nonsig_cells <- all_cells[which(allmisscols == TRUE)]
nonsig_cells <- nonsig_cells[-c(1,2,3,4,5,6,7)]
nonsig_cells
#######
print("# plot presence/absence of valid interaction with all cell types")
#######
mytitle <- gsub("_", " ", Atype)
#data for presence/absence
pa_dat <- as.data.frame(cbind(all_cells, allmisscols))
pa_dat <- pa_dat[-c(1,2,3,4,5,6,7),]
#reorder based on T/F status
pa_dat <- pa_dat %>% arrange(desc(allmisscols)) %>%    
  mutate(all_cells=factor(all_cells, levels=all_cells)) 
#change name of T/F levels
levels(pa_dat$allmisscols) <- list("Present" = "FALSE",
                                   "Absent"="TRUE")
#get column for presence z-score sign info
#convert df with just interacting region to long format
#VinterL <- gather(Vinter, cell, zscore, 8:length(colnames(Vinter)), factor_key=TRUE)
VinterL <- Vinter
VinterL$zsign <- VinterL$zscore
VinterL$zsign[VinterL$zsign>0]  <- "Pos" 
VinterL$zsign[VinterL$zsign<0]  <- "Neg" 
VinterL$zsign[is.na(VinterL$zsign)] <- "Abs"
VinterL <- VinterL %>% select(cell,zsign)
colnames(VinterL) <- c("all_cells","zsign")
summary(pa_dat_m)
levels(pa_dat_m$allmisscols)
pa_dat_m <- merge(pa_dat, VinterL, by = "all_cells")
p <- (ggplot(pa_dat_m, aes(y=all_cells, x=allmisscols, shape =allmisscols, fill = zsign))
      + geom_point(size=4)
      #      + coord_flip()
      + labs(title = mytitle,
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = "Valid interaction status in all interactions", y = "")
      #         tag = "A")
      + scale_shape_manual(values=c(21,1))
      + scale_fill_manual(values =c("Pos" = "#26C485", "Neg" = "#88E0FB", "Abs" = "black"))
      + scale_color_manual(values =c("Pos" = "black", "Neg" = "black", "Abs" = "black"))
      #  + facet_grid(seqDep ~ cell )
      + labs(shape = "", linetype = "")
      + theme(legend.position = "none")
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_present_absent_allInters",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()

#sig interactions
print("cells/datasets that have interactions btwn validated chromosomes")
chrs_df2 <- chrs_df_allInters
chrs_df2$zscore[as.numeric(chrs_df2$pvalue)>=0.05]  <- NA 
chrs_df <- chrs_df2 %>% select(-pvalue) %>% spread(key = cell, value = zscore)
summary(chrs_df)
allmisscols <- sapply(chrs_df, function(x) all(is.na(x) | x == '' ))
all_cells <- colnames(chrs_df)
sig_cells <- all_cells[which(allmisscols == FALSE)]
sig_cells <- sig_cells[-c(1,2,3,4,5,6,7)]
sig_cells
print("cells/datasets that have interactions btwn validated chromosomes")
nonsig_cells <- all_cells[which(allmisscols == TRUE)]
nonsig_cells <- nonsig_cells[-c(1,2,3,4,5,6,7)]
nonsig_cells
#######
print("# plot presence/absence of valid interaction with all cell types")
#######
mytitle <- gsub("_", " ", Atype)
#data for presence/absence
pa_dat <- as.data.frame(cbind(all_cells, allmisscols))
pa_dat <- pa_dat[-c(1,2,3,4,5,6,7),]
#reorder based on T/F status
pa_dat <- pa_dat %>% arrange(desc(allmisscols)) %>%    
  mutate(all_cells=factor(all_cells, levels=all_cells)) 
#change name of T/F levels
levels(pa_dat$allmisscols) <- list("Present" = "FALSE",
                                   "Absent"="TRUE")
#get column for presence z-score sign info
VinterL <- Vinter
VinterL$zsign <- VinterL$zscore
VinterL$zsign[VinterL$zsign>0]  <- "Pos" 
VinterL$zsign[VinterL$zsign<0]  <- "Neg" 
VinterL$zsign[is.na(VinterL$zsign)] <- "Abs"
VinterL <- VinterL %>% select(cell,zsign)
colnames(VinterL) <- c("all_cells","zsign")
pa_dat_m <- merge(pa_dat, VinterL, by = "all_cells")
p <- (ggplot(pa_dat_m, aes(y=all_cells, x=allmisscols, shape =allmisscols, fill = zsign))
      + geom_point(size=4)
      #      + coord_flip()
      + labs(title = mytitle,
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = "Valid interaction status in significant interactions", y = "")
      #         tag = "A")
      + scale_shape_manual(values=c(21,1))
      + scale_fill_manual(values =c("Pos" = "#26C485", "Neg" = "#88E0FB", "Abs" = "black"))
      + scale_color_manual(values =c("Pos" = "black", "Neg" = "black", "Abs" = "black"))
      #  + facet_grid(seqDep ~ cell )
      + labs(shape = "", linetype = "")
      + theme(legend.position = "none")
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_present_absent_sigInters",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()


print("#######")
print("# Validated interaction chromosomes (not specific region)")
print("#######")
#df with ALL interactions btwn validated chroms
chrs_df <- dat2[grep(chrs[1], dat2$ID), ]
chrs_df <- chrs_df[grep(chrs[2],chrs_df$ID),]
chrs_df <- chrs_df %>% select(-pvalue) %>% spread(key = cell, value = zscore)
summary(chrs_df)

print("cells/datasets that have interactions btwn validated chromosomes")
allmisscols <- sapply(chrs_df, function(x) all(is.na(x) | x == '' ))
all_cells <- colnames(chrs_df)
sig_cells <- all_cells[which(allmisscols == FALSE)]
sig_cells <- sig_cells[-c(1,2,3,4,5,6,7)]
sig_cells
print("cells/datasets that have NO sig interactions btwn validated chromosomes")
nonsig_cells <- all_cells[which(allmisscols == TRUE)]
nonsig_cells <- nonsig_cells[-c(1,2,3,4,5,6,7)]
nonsig_cells

#######
print("# only plot our 6 test cell types")
#######
#df with just these cells
sixcells_all <- filter(dat2, cell == "Dorsolateral_prefrontal_cortex" | cell == "Small_bowel_Schmitt" |
                         cell == "Aorta" | cell == "Right_ventricle_Schmitt" | cell == "Cardiomyocites_primitive_rep1" |
                         cell == "H1hESC_Oksuz")
head(sixcells_all)
sixcells_Vinter <- filter(Vinter, cell == "Dorsolateral_prefrontal_cortex" | cell == "Small_bowel_Schmitt" |
                         cell == "Aorta" | cell == "Right_ventricle_Schmitt" | cell == "Cardiomyocites_primitive_rep1" |
                         cell == "H1hESC_Oksuz")

head(sixcells_Vinter)
sixcells_chr <- select(chrs_df,c("ID",
                                 "Dorsolateral_prefrontal_cortex", 
                                 "Small_bowel_Schmitt", 
                                 "Aorta", 
                                 "Right_ventricle_Schmitt", 
                                 "Cardiomyocites_primitive_rep1",
                                 #                                 "Spleen"))
                                 "H1hESC_Oksuz"))
head(sixcells_chr)

#wide to long
##all data (six cells)
#all_long <- gather(sixcells_all, cell, zscore, 2:length(colnames(sixcells_all)), factor_key=TRUE)
all_long <- sixcells_all %>% select(-pvalue, -chrA, -st1, -end1, -chrB, -st2, -end2)
all_long$validType <- rep("all", length(all_long$ID))

#valid inter (six cells)
#valid_inter_long <- gather(sixcells_Vinter, cell, zscore, 2:length(colnames(sixcells_Vinter)), factor_key=TRUE)
valid_inter_long <- sixcells_Vinter %>% select(-pvalue, -chrA, -st1, -end1, -chrB, -st2, -end2)
valid_inter_long$validType <- rep("validInter", length(valid_inter_long$ID))

#valid chrs (six cells)
valid_chrs_long <- gather(sixcells_chr, cell, zscore, 2:length(colnames(sixcells_chr)), factor_key=TRUE)
valid_chrs_long$validType <- rep("validChr", length(valid_chrs_long$ID))

#combine the dfs
plot_df <- rbind(all_long, valid_inter_long, valid_chrs_long)

# rename levels
plot_df$validType <- as.factor(plot_df$validType)
plot_df$cell <- as.factor(plot_df$cell)
levels(plot_df$validType) <- list("Interactions from valid inter" = "validInter",
                                  "Interactions from chrs"="validChr",
                                  "All interactions"="all")
levels(plot_df$cell) <- list("hESCDekker" = "H1hESC_Oksuz",
                             #levels(plot_df$cell) <- list("Spleen" = "Spleen",
                             "CardPrimRep1" = "Cardiomyocites_primitive_rep1",
                             "RightVentricle" = "Right_ventricle_Schmitt",
                             "Aorta"="Aorta",
                             "SmallBowel"="Small_bowel_Schmitt",
                             "DorsoPreCort" = "Dorsolateral_prefrontal_cortex")
#new col for seq depth
plot_df$seqDep <- plot_df$cell
levels(plot_df$seqDep) <- list("High" = "hESCDekker",
                               #levels(plot_df$seqDep) <- list("High" = "Spleen",
                               "High" = "CardPrimRep1",
                               "Medium" = "RightVentricle",
                               "Medium"="Aorta",
                               "Low" = "SmallBowel",
                               "Low" = "DorsoPreCort")

#function for split violin plot
# https://stackoverflow.com/a/45614547
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, draw_group = function(self, data, ..., draw_quantiles = NULL){
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1,'group']
  newdata <- plyr::arrange(transform(data, x = if(grp%%2==1) xminv else xmaxv), if(grp%%2==1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x']) 
  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 
                                              1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function (mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, show.legend = show.legend, inherit.aes = inherit.aes, params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}



print("#plot above data for valid chroms")
target <- c("Interactions from chrs", "All interactions")
plot_d <- plot_df %>% filter(validType %in% target)
interLab <- paste("Interactions between",roi1_inter1$seqnames,"and",roi2_inter1$seqnames)
#p <- (ggplot(plot_d, aes(cell, zscore, fill = validType))
#      + geom_density(alpha=.3, stat= "identity")
#      + coord_flip()
#      + labs(title = mytitle,
#             #         subtitle = "Plot of length by dose",
#             #         caption = "Data source: ToothGrowth",
#             x = "", y = "z-score")
#      #         tag = "A")
#      #  + scale_fill_manual(values =c("plum4", "cadetblue"))
#      + scale_fill_manual(values =c("plum4", "cadetblue"),labels=c(interLab, 'All significant interactions'))
#      #  + facet_grid(seqDep ~ cell )
#)
levels(plot_d$validType) <- list("All interactions" = "All interactions",
                                 "Interactions from chrs" = "Interactions from chrs")
p <- (ggplot(plot_d, aes(x=zscore, fill = validType))
      #  + geom_split_violin()
      + stat_density(alpha=.6,position="identity")#identity = based on counts of data, height proportional to total
      #  + coord_flip()
      + labs(title = mytitle,
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = "z-score", y = "Density")
      #         tag = "A")
     #   + scale_fill_manual(values =c("plum4", "cadetblue"))
      + scale_fill_manual(values =c("grey", "cadetblue"),labels=c('All interactions', interLab))
      + facet_grid(cell~. )
)
#p
f_name <- gsub(" ","",paste("6testCells_valid_interaction_chroms_density_allInters",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
#print("#test between means of all interactions and interactions btwn valid pair of chroms")
#print("#Test each group for normality")
#print("sig = reject normality null")
#plot_d$validType <- as.factor(plot_d$validType)
#plot_d %>%
#  group_by(validType) %>%
#  summarise(W = ad.test(zscore)$statistic,
#            p.value = ad.test(zscore)$p.value)
print("#Perform the Mann-Whitney test (when dists are not normal)")
print("sig = mean is diff btwn groups")
wilcox.test(zscore ~ validType, data=plot_d, na.rm=TRUE, paired=FALSE, exact=FALSE, conf.int=TRUE)


######
#sig interactions
######
#df with ALL sig interactions btwn validated chroms
chrs_df <- dat2[grep(chrs[1], dat2$ID), ]
chrs_df$zscore[as.numeric(chrs_df$pvalue)>=0.05]  <- NA 
chrs_df <- chrs_df[grep(chrs[2],chrs_df$ID),]
chrs_df <- chrs_df %>% select(-pvalue) %>% spread(key = cell, value = zscore)
summary(chrs_df)

print("cells/datasets that have sig interactions btwn validated chromosomes")
allmisscols <- sapply(chrs_df, function(x) all(is.na(x) | x == '' ))
all_cells <- colnames(chrs_df)
sig_cells <- all_cells[which(allmisscols == FALSE)]
sig_cells <- sig_cells[-c(1,2,3,4,5,6,7)]
sig_cells
print("cells/datasets that have NO sig interactions btwn validated chromosomes")
nonsig_cells <- all_cells[which(allmisscols == TRUE)]
nonsig_cells <- nonsig_cells[-c(1,2,3,4,5,6,7)]
nonsig_cells

#######
print("# only plot our 6 test cell types")
#######
#df with just these cells
sixcells_all <- filter(dat2, cell == "Dorsolateral_prefrontal_cortex" | cell == "Small_bowel_Schmitt" |
                         cell == "Aorta" | cell == "Right_ventricle_Schmitt" | cell == "Cardiomyocites_primitive_rep1" |
                         cell == "H1hESC_Oksuz")
sixcells_all$zscore[as.numeric(sixcells_all$pvalue)>=0.05]  <- NA 
head(sixcells_all)
sixcells_Vinter <- filter(Vinter, cell == "Dorsolateral_prefrontal_cortex" | cell == "Small_bowel_Schmitt" |
                            cell == "Aorta" | cell == "Right_ventricle_Schmitt" | cell == "Cardiomyocites_primitive_rep1" |
                            cell == "H1hESC_Oksuz")
sixcells_Vinter$zscore[as.numeric(sixcells_Vinter$pvalue)>=0.05]  <- NA 

head(sixcells_Vinter)
sixcells_chr <- select(chrs_df,c("ID",
                                 "Dorsolateral_prefrontal_cortex", 
                                 "Small_bowel_Schmitt", 
                                 "Aorta", 
                                 "Right_ventricle_Schmitt", 
                                 "Cardiomyocites_primitive_rep1",
                                 #                                 "Spleen"))
                                 "H1hESC_Oksuz"))
head(sixcells_chr)

#wide to long
##all data (six cells)
#all_long <- gather(sixcells_all, cell, zscore, 2:length(colnames(sixcells_all)), factor_key=TRUE)
all_long <- sixcells_all %>% select(-pvalue, -chrA, -st1, -end1, -chrB, -st2, -end2)
all_long$validType <- rep("all", length(all_long$ID))

#valid inter (six cells)
#valid_inter_long <- gather(sixcells_Vinter, cell, zscore, 2:length(colnames(sixcells_Vinter)), factor_key=TRUE)
valid_inter_long <- sixcells_Vinter %>% select(-pvalue, -chrA, -st1, -end1, -chrB, -st2, -end2)
valid_inter_long$validType <- rep("validInter", length(valid_inter_long$ID))

#valid chrs (six cells)
valid_chrs_long <- gather(sixcells_chr, cell, zscore, 2:length(colnames(sixcells_chr)), factor_key=TRUE)
valid_chrs_long$validType <- rep("validChr", length(valid_chrs_long$ID))

#combine the dfs
plot_df <- rbind(all_long, valid_inter_long, valid_chrs_long)

# rename levels
plot_df$validType <- as.factor(plot_df$validType)
plot_df$cell <- as.factor(plot_df$cell)
levels(plot_df$validType) <- list("Interactions from valid inter" = "validInter",
                                  "Interactions from chrs"="validChr",
                                  "All significant interactions"="all")
levels(plot_df$cell) <- list("hESCDekker" = "H1hESC_Oksuz",
                             #levels(plot_df$cell) <- list("Spleen" = "Spleen",
                             "CardPrimRep1" = "Cardiomyocites_primitive_rep1",
                             "RightVentricle" = "Right_ventricle_Schmitt",
                             "Aorta"="Aorta",
                             "SmallBowel"="Small_bowel_Schmitt",
                             "DorsoPreCort" = "Dorsolateral_prefrontal_cortex")
#new col for seq depth
plot_df$seqDep <- plot_df$cell
levels(plot_df$seqDep) <- list("High" = "hESCDekker",
                               #levels(plot_df$seqDep) <- list("High" = "Spleen",
                               "High" = "CardPrimRep1",
                               "Medium" = "RightVentricle",
                               "Medium"="Aorta",
                               "Low" = "SmallBowel",
                               "Low" = "DorsoPreCort")


print("#plot above data for valid chroms")
target <- c("Significant interactions from chrs", "All significant interactions")
plot_d <- plot_df %>% filter(validType %in% target)
interLab <- paste("Significant interactions between",roi1_inter1$seqnames,"and",roi2_inter1$seqnames)
#p <- (ggplot(plot_d, aes(cell, zscore, fill = validType))
#      + geom_density(alpha=.3, stat= "identity")
#      + coord_flip()
#      + labs(title = mytitle,
#             #         subtitle = "Plot of length by dose",
#             #         caption = "Data source: ToothGrowth",
#             x = "", y = "z-score")
#      #         tag = "A")
#      #  + scale_fill_manual(values =c("plum4", "cadetblue"))
#      + scale_fill_manual(values =c("plum4", "cadetblue"),labels=c(interLab, 'All significant interactions'))
#      #  + facet_grid(seqDep ~ cell )
#)
levels(plot_d$validType) <- list("All significant interactions" = "All significant interactions",
                                 "Significant interactions from chrs" = "Significant interactions from chrs")
p <- (ggplot(plot_d, aes(x=zscore, fill = validType))
      #  + geom_split_violin()
      + stat_density(alpha=.6,position="identity")#identity = based on counts of data, height proportional to total
      #  + coord_flip()
      + labs(title = mytitle,
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = "z-score", y = "Density")
      #         tag = "A")
      #  + scale_fill_manual(values =c("plum4", "cadetblue"))
      + scale_fill_manual(values =c("plum4", "cadetblue"),labels=c(interLab, 'All significant interactions'))
      + facet_grid(cell~. )
)
#p
f_name <- gsub(" ","",paste("6testCells_valid_interaction_chroms_density_sigInters",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()

##############
# zscore
##############
##############
# heatmap along chromosome positions: ALL INTERACTIONS (including non-sig)
##############
head(hm_zscore_df)
hm_zscore_df$st1 <- as.numeric(hm_zscore_df$st1)
hm_zscore_df$st2 <- as.numeric(hm_zscore_df$st2)
summary(hm_zscore_df)
hm_zscore_df$chrs <- paste0(hm_zscore_df$chrA,hm_zscore_df$chrB)
dirA <- paste0(unique(roi1$chrom),unique(roi2$chrom))
dirB <- paste0(unique(roi2$chrom),unique(roi1$chrom))
hm_df <- hm_zscore_df %>% filter(chrs == dirA | chrs == dirB)
head(hm_df)
testdf <- hm_df %>% filter(cell == "Aorta")
xchr <- gsub("chr", "", unique(roi1$chrom))
ychr <- gsub("chr", "", unique(roi2$chrom))
#hm <- (ggplot(testdf, aes(x=st1, y=st2, fill= zscore)) 
hm <- (ggplot(hm_df, aes(x=st1, y=st2, fill= zscore)) 
           + geom_tile()
      + labs(title = "Distribution of z-scores along valid interacting chromosome (all interactions)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", xchr, " Genomic Position"),
             y = paste0("Chromosome ", ychr, " Genomic Position"))
      + scale_fill_hp(discrete = FALSE, option = "Always", name = "z-score")
      + facet_grid(cell ~.)
      + theme(panel.spacing = unit(0, "lines"),
              strip.text.y.right = element_text(angle = 0), #rotate facet labels
              strip.background = element_rect(fill = "white"),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_heatmap_allInters",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
hm
dev.off()

##############
# heatmap along chromosome positions: SIG INTERACTIONS (only sig)
##############
hm_dfsig <- hm_df %>% filter(pvalue <= 0.05)
#hm <- (ggplot(testdf, aes(x=st1, y=st2, fill= zscore)) 
hm <- (ggplot(hm_dfsig, aes(x=st1, y=st2, fill= zscore)) 
       + geom_tile()
       + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions)",
              #         subtitle = "Plot of length by dose",
              #         caption = "Data source: ToothGrowth",
              x = paste0("Chromosome ", xchr, " Genomic Position"),
              y = paste0("Chromosome ", ychr, " Genomic Position"))
       + scale_fill_hp(discrete = FALSE, option = "Always", name = "z-score")
       + facet_grid(cell ~.)
       + theme(panel.spacing = unit(0, "lines"),
               strip.text.y.right = element_text(angle = 0), #rotate facet labels
               strip.background = element_rect(fill = "white"),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_heatmap_sigInters",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
hm
dev.off()

##############
# Pvalue
##############
##############
# heatmap along chromosome positions: ALL INTERACTIONS (including non-sig)
##############
#hm <- (ggplot(testdf, aes(x=st1, y=st2, fill= zscore)) 
hm <- (ggplot(hm_df, aes(x=st1, y=st2, fill= pvalue)) 
       + geom_tile()
       + labs(title = "Distribution of p-value along valid interacting chromosome (all interactions)",
              #         subtitle = "Plot of length by dose",
              #         caption = "Data source: ToothGrowth",
              x = paste0("Chromosome ", xchr, " Genomic Position"),
              y = paste0("Chromosome ", ychr, " Genomic Position"))
       + scale_fill_hp(discrete = FALSE, option = "Always", name = "p-value")
       + facet_grid(cell ~.)
       + theme(panel.spacing = unit(0, "lines"),
               strip.text.y.right = element_text(angle = 0), #rotate facet labels
               strip.background = element_rect(fill = "white"),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_pvalue_heatmap_allInters",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
hm
dev.off()

##############
# heatmap along chromosome positions: SIG INTERACTIONS (only sig)
##############
#hm <- (ggplot(testdf, aes(x=st1, y=st2, fill= zscore)) 
hm <- (ggplot(hm_dfsig, aes(x=st1, y=st2, fill= pvalue)) 
       + geom_tile()
       + labs(title = "Distribution of p-value along valid interacting chromosome (significant interactions)",
              #         subtitle = "Plot of length by dose",
              #         caption = "Data source: ToothGrowth",
              x = paste0("Chromosome ", xchr, " Genomic Position"),
              y = paste0("Chromosome ", ychr, " Genomic Position"))
       + scale_fill_hp(discrete = FALSE, option = "Always", name = "p-value")
       + facet_grid(cell ~.)
       + theme(panel.spacing = unit(0, "lines"),
               strip.text.y.right = element_text(angle = 0), #rotate facet labels
               strip.background = element_rect(fill = "white"),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_pvalue_heatmap_sigInters",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
hm
dev.off()




#dealing with cases where the specified interaction is not sig in any cell type
sixcells_Vinter <- gather(sixcells_Vinter, cell, zscore, 2:length(colnames(sixcells_all)), factor_key=TRUE)

if (dim(Vinter)[1] == 0){
  print("dataframe is empty :( your specified interaction was not significant in any of your cell types/data")
} else {
  print("sig interactions in at least ONE cell type!")
#  #######
#  # Validated interaction specific region
#  #######
#  #plot above data for valid chroms
#  plot_d <- plot_df %>% filter(validType == "All significant interactions")
#  inter_d <- plot_df %>% filter(validType == "Interactions from valid inter")
#  head(inter_d)
#  p <- (ggplot(plot_d, aes(x=zscore))
##        + geom_density(fill = "cadetblue")
#      + stat_density(fill = "cadetblue",position="identity")#identity = based on counts of data, height proportional to total
#        #        + coord_flip()
#        + geom_vline(data = sixcells_Vinter, aes(xintercept = zscore, 
#                                                 color = cell), size=1.5)
#        + labs(title = mytitle,
#               #         subtitle = "Plot of length by dose",
#               #         caption = "Data source: ToothGrowth",
#               x = "z-score", y = "density")
#        #         tag = "A")
#        #        + scale_color_manual(values =c("plum4"))
#        #  + facet_grid(seqDep ~ cell )
#  )
#    p
#  f_name <- gsub(" ","",paste("6testCells_valid_interaction_density_",Atype,".pdf"))
#  pdf(f_name, width = 14, height = 8)
#  print(p)
#  dev.off()
}

