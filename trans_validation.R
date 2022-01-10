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
library(tidyr)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(harrypotter, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.3.6.1")
library(nortest, lib="/hpf/largeprojects/pmaass/programs/Rlib/R.3.6.1") #for normality test with large sample size
library(hexbin)
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
#xchr <- "X"
#ychr <- "22"
#xstart <- 131688779 /1000000
#ystart <- 39519695 /1000000
#library(harrypotter)
#library(factoextra)
#library(hexbin)
##dat <- read.table("23Jul21.primary.trans.1MB.zscores.txt", header = TRUE)
##dat <- read.table("23Jul21.primary.trans.1MB.zscores.pairwise.txt", header = TRUE)
##dat <- read.table(dat_file, header = TRUE)
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
print(roi1)
roi1$chrom <- as.factor(roi1$chrom)
#roi2
#roi2 <- read.table("SOX9.bed", header = FALSE)
roi2 <- read.table(roi2_file)
colnames(roi2) <- c("chrom","start","end")
roi2$chrom <- as.factor(roi2$chrom)
print(roi2)

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
chrs <- paste(chrs, ".", sep = "") 
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
pa_dat$allmisscols <- factor(pa_dat$allmisscols)
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
VinterL$zsign[VinterL$zsign>=0]  <- "Pos" 
VinterL$zsign[VinterL$zsign<0]  <- "Neg" 
VinterL$zsign[is.na(VinterL$zsign)] <- "Abs"
VinterL <- VinterL %>% select(cell,zsign)
colnames(VinterL) <- c("all_cells","zsign")
pa_dat_m <- merge(pa_dat, VinterL, by = "all_cells")
summary(pa_dat_m)
#pa_dat_m$allmisscols <- ifelse(pa_dat_m$allmisscols == TRUE, "Present", "Absent")
levels(pa_dat_m$allmisscols)
pa_dat_m$allmisscols <- as.factor(pa_dat_m$allmisscols)
pa_dat_m$zsign <- as.factor(pa_dat_m$zsign)
levels(pa_dat_m$zsign)
print("# odd values for zsign")
pa_dat_m %>% filter(zsign == "0")
p <- (ggplot(pa_dat_m, aes(y=all_cells, x=allmisscols, shape =allmisscols, fill = zsign))
      + geom_point(size=4)
      #      + coord_flip()
      + labs(title = mytitle,
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = "Valid interaction status in all interactions", y = "")
      #         tag = "A")
#      + scale_shape_manual(values=c(21,1))
#      + scale_fill_manual(values =c("Pos" = "#26C485", "Neg" = "#88E0FB", "Abs" = "black"))
#      + scale_color_manual(values =c("Pos" = "black", "Neg" = "black", "Abs" = "black"))
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
chrs_df2$zscore[as.numeric(chrs_df2$pvalue)<=0.05]  <- NA 
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
pa_dat$allmisscols <- factor(pa_dat$allmisscols)
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
#      + scale_shape_manual(values=c(21,1))
#      + scale_fill_manual(values =c("Pos" = "#26C485", "Neg" = "#88E0FB", "Abs" = "black"))
#      + scale_color_manual(values =c("Pos" = "black", "Neg" = "black", "Abs" = "black"))
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
#histogram
p <- (ggplot(plot_d, aes(x=zscore, fill = validType))
      #  + geom_split_violin()
#      + stat_density(alpha=.6,position="identity")#identity = based on counts of data, height proportional to total
#      + geom_density(alpha=.6,position="stack")#stack = based on counts of data, height proportional to total
      + geom_histogram(position="identity",alpha=.6)#stack = based on counts of data, height proportional to total
      #  + coord_flip()
      + labs(title = mytitle,
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = "z-score", y = "Count")
      #         tag = "A")
     #   + scale_fill_manual(values =c("plum4", "cadetblue"))
      + scale_fill_manual(values =c("grey", "cadetblue"),labels=c('All interactions', interLab))
      + facet_grid(cell~. )
      + theme(legend.title = element_blank())
)
#p
f_name <- gsub(" ","",paste("6testCells_valid_interaction_chroms_histogram_allInters",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
#density
p <- (ggplot(plot_d, aes(x=zscore, fill = validType))
      #  + geom_split_violin()
      + stat_density(alpha=.6,position="identity")#identity = based on counts of data, height proportional to total
#      + geom_density(alpha=.6,position="stack")#stack = based on counts of data, height proportional to total
#      + geom_histogram(position="identity",alpha=.6)#stack = based on counts of data, height proportional to total
      #  + coord_flip()
      + labs(title = mytitle,
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = "z-score", y = "Count")
      #         tag = "A")
     #   + scale_fill_manual(values =c("plum4", "cadetblue"))
      + scale_fill_manual(values =c("grey", "cadetblue"),labels=c('All interactions', interLab))
      + facet_grid(cell~. )
      + theme(legend.title = element_blank())
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
print("#df with ALL sig interaction btwn validated chroms")
chrs_df <- dat2[grep(chrs[1], dat2$ID), ]
chrs_df$zscore[as.numeric(chrs_df$pvalue)<=0.05]  <- NA 
chrs_df <- chrs_df[grep(chrs[2],chrs_df$ID),]
chrs_df <- chrs_df %>% select(-pvalue) %>% spread(key = cell, value = zscore)
summary(chrs_df)
head(chrs_df)

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
sixcells_all$zscore[as.numeric(sixcells_all$pvalue)<=0.05]  <- NA 
head(sixcells_all)
sixcells_Vinter <- filter(Vinter, cell == "Dorsolateral_prefrontal_cortex" | cell == "Small_bowel_Schmitt" |
                            cell == "Aorta" | cell == "Right_ventricle_Schmitt" | cell == "Cardiomyocites_primitive_rep1" |
                            cell == "H1hESC_Oksuz")
sixcells_Vinter$zscore[as.numeric(sixcells_Vinter$pvalue)<=0.05]  <- NA 

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
target <- c("Interactions from chrs", "All significant interactions")
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
                                 "Interactions from chrs" = "Significant interactions from chrs")
head(plot_d)
summary(plot_d)
#histogram
p <- (ggplot(plot_d, aes(x=zscore, fill = validType))
      #  + geom_split_violin()
#      + stat_density(alpha=.6,position="identity")#identity = based on counts of data, height proportional to total
#      + geom_density(alpha=.6,position="stack")#stack = based on counts of data, height proportional to total
      + geom_histogram(position="identity", alpha=.5)#identiry = overlapping histograms
#  + coord_flip()
      + labs(title = mytitle,
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = "z-score", y = "Count")
      #         tag = "A")
      #  + scale_fill_manual(values =c("plum4", "cadetblue"))
      + scale_fill_manual(values =c("grey", "cadetblue"),labels=c('All significant interactions',interLab))
      + facet_grid(cell~. )
      + theme(legend.title = element_blank())
)
#p
f_name <- gsub(" ","",paste("6testCells_valid_interaction_chroms_histogram_sigInters",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
#density
p <- (ggplot(plot_d, aes(x=zscore, fill = validType))
      #  + geom_split_violin()
     + stat_density(alpha=.6,position="identity")#identity = based on counts of data, height proportional to total
#      + geom_density(alpha=.6,position="stack")#stack = based on counts of data, height proportional to total
#      + geom_histogram(position="identity", alpha=.5)#identiry = overlapping histograms
#  + coord_flip()
      + labs(title = mytitle,
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = "z-score", y = "Count")
      #         tag = "A")
      #  + scale_fill_manual(values =c("plum4", "cadetblue"))
      + scale_fill_manual(values =c("grey", "cadetblue"),labels=c('All significant interactions',interLab))
      + facet_grid(cell~. )
      + theme(legend.title = element_blank())
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
print("# heatmap along chromosome positions: ALL INTERACTIONS (including non-sig)")
##############
head(hm_zscore_df)
hm_zscore_df$st1 <- as.numeric(hm_zscore_df$st1)
hm_zscore_df$st2 <- as.numeric(hm_zscore_df$st2)
summary(hm_zscore_df)
hm_zscore_df$chrs <- paste0(hm_zscore_df$chrA,hm_zscore_df$chrB)
head(hm_zscore_df)
dirA <- paste0(unique(roi1$chrom),unique(roi2$chrom))
print(dirA)
dirB <- paste0(unique(roi2$chrom),unique(roi1$chrom))
print(dirB)
hm_df <- hm_zscore_df %>% filter(chrs == dirA | chrs == dirB)
head(hm_df)
print("TEST chrA")
CA <- unique(hm_df$chrA)
CA
print("TEST chrB")
CB <- unique(hm_df$chrB)
CB
testdf <- hm_df %>% filter(cell == "Aorta")
xchr <- gsub("chr", "", CA)
ychr <- gsub("chr", "", CB)
#hm <- (ggplot(testdf, aes(x=st1, y=st2, fill= zscore)) 

#get max of valid chromosomes for axis expansion
#chrom info: centromere (midpoint calculated from UCSC, aprox), chrom class
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
head(chrInf)
xmax <- chrInf$size[match(CA,chrInf$chrom)] / 1000000
ymax <- chrInf$size[match(CB,chrInf$chrom)] / 1000000
# getting the valid inter positions to match with the x and y chroms in the graphs
valid_df <- rbind(roi1, roi2)
xdf <- valid_df %>% filter(chrom == paste0("chr",xchr))
xstart <- xdf$start /1000000
xdf
xstart
ydf <- valid_df %>% filter(chrom == paste0("chr",ychr))
ystart <- ydf$start /1000000
ydf
ystart
#scale genomic position by 1Mb
hm_df$st1 <- hm_df$st1/1000000
hm_df$st2 <- hm_df$st2/1000000
#gtext = textGrob("xyz", y = -5, gp = gpar(col = "red"))
#gline = linesGrob(y = c(-.02, .02),  gp = gpar(col = "red", lwd = 2))
hm <- (ggplot(hm_df, aes(x=st1, y=st2, fill= zscore)) 
           + geom_tile(width = 1, height = 1)
      + labs(title = "Distribution of z-scores along valid interacting chromosome (all interactions)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", xchr, " Genomic Position [Mb]"),
             y = paste0("Chromosome ", ychr, " Genomic Position"))
     # + annotate("rect", xmin = 39, xmax = 40, ymin = -2, ymax = -1, fill = "red")
      + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
      + facet_grid(cell ~.)
      # expand axis limits so whole chrom len is accounted for
      + expand_limits(y=c(0,ymax), x = c(0,xmax))
#      + annotation_custom(gtext, xmin=30, xmax=30, ymin=-Inf, ymax=Inf)
#      + annotation_custom(gline, xmin=30, xmax=30, ymin=-Inf, ymax=Inf)
      + geom_vline(xintercept = xstart, colour = "black")
      + scale_y_continuous(expand = c(0, 0))
      + scale_x_continuous(expand = c(0, 0))
      + theme(panel.spacing = unit(0, "lines"),
              strip.text.y.right = element_text(angle = 0), #rotate facet labels
              strip.background = element_rect(fill = "white"),
             # axis.text.x = element_blank(),
             # axis.ticks.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_heatmap_allInters",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
#g = ggplotGrob(hm)
#g$layout$clip[g$layout$name=="panel"] <- "off"
#grid.draw(g)
hm
dev.off()

##############
# heatmap along chromosome positions: SIG INTERACTIONS (only sig)
##############
hm_dfsig <- hm_df %>% filter(pvalue <= 0.05)
#hm <- (ggplot(testdf, aes(x=st1, y=st2, fill= zscore)) 
hm <- (ggplot(hm_dfsig, aes(x=st1, y=st2, fill= zscore)) 
       + geom_tile(width = 1, height = 1)
       + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions)",
              #         subtitle = "Plot of length by dose",
              #         caption = "Data source: ToothGrowth",
              x = paste0("Chromosome ", xchr, " Genomic Position [Mb]"),
              y = paste0("Chromosome ", ychr, " Genomic Position"))
       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
      + geom_vline(xintercept = xstart, colour = "black")
       + facet_grid(cell ~.)
      # expand axis limits so whole chrom len is accounted for
      + expand_limits(y=c(0,ymax), x = c(0,xmax))
      + scale_y_continuous(expand = c(0, 0))
      + scale_x_continuous(expand = c(0, 0))
       + theme(panel.spacing = unit(0, "lines"),
               strip.text.y.right = element_text(angle = 0), #rotate facet labels
               strip.background = element_rect(fill = "white"),
#               axis.text.x = element_blank(),
#               axis.ticks.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_heatmap_sigInters",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
hm
dev.off()
#having the opposite chrom on the x axis
hm <- (ggplot(hm_dfsig, aes(x=st2, y=st1, fill= zscore)) 
       + geom_tile(width = 1, height = 1)
       + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions)",
              #         subtitle = "Plot of length by dose",
              #         caption = "Data source: ToothGrowth",
              x = paste0("Chromosome ", ychr, " Genomic Position [Mb]"),
              y = paste0("Chromosome ", xchr, " Genomic Position"))
      + geom_vline(xintercept = ystart, colour = "black")
       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
       + facet_grid(cell ~.)
       # expand axis limits so whole chrom len is accounted for
       + expand_limits(y=c(0,xmax), x = c(0,ymax))
      + scale_y_continuous(expand = c(0, 0))
      + scale_x_continuous(expand = c(0, 0))
       + theme(panel.spacing = unit(0, "lines"),
               strip.text.y.right = element_text(angle = 0), #rotate facet labels
               strip.background = element_rect(fill = "white"),
               #               axis.text.x = element_blank(),
               #               axis.ticks.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_heatmap_sigInters_YX",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
hm
dev.off()

print("##########")
print("# linear highly interacting regions between valid chrom pair")
#only sig inters btwn that valid chroms
hir_df <- hm_df %>% filter(pvalue <= 0.05) #filter for sig interactions
length(hir_df$zscore)
aorta_df <- hir_df %>% filter(cell == "Aorta")
#p <- (ggplot(aorta_df, aes(x=st1, y=zscore)) 
print("# pts and line")
print("# two graphs for pos and neg")
hir_df <- hir_df %>% mutate(sign = ifelse(zscore >= 0, "pos", "neg"))
head(hir_df)
p <- (ggplot(hir_df, aes(x=st1, y=zscore, fill=sign)) 
      + geom_point(alpha = 0.4, aes(color = sign))
      + geom_vline(aes(xintercept = xstart), colour = "red")
      + geom_smooth(aes(colour = sign), method = 'loess', formula = y ~ x)
      + labs(title = "Distribution of z-scores between valid interacting chromosomes (significant interactions)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", xchr, " Genomic Position [Mb]"),
             y = "z-score",
             fill = "z-score sign")
      + scale_fill_manual(values =c("pos" = "#EE9B00", "neg" = "#005F73"), name = "z-score sign",labels = c("pos" = "Positive","neg" ="Negative"))
      + scale_colour_manual(values =c("pos" = "#EE9B00", "neg" = "#005F73"), name = "z-score sign", labels = c("pos" = "Positive","neg" ="Negative"))
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,xmax))
      + scale_x_continuous(expand = c(0, 0))
#      + facet_wrap(.~ sign, scales = "free")
      #       + theme(panel.spacing = unit(0, "lines"),
      #               strip.text.y.right = element_text(angle = 0), #rotate facet labels
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_pts_posneg_sigInters_all_cells",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
p <- (ggplot(hir_df, aes(x=st2, y=zscore, fill=sign)) 
      + geom_point(alpha = 0.4, aes(color = sign))
      + geom_vline(aes(xintercept = ystart), colour = "red")
      + geom_smooth(aes(colour = sign), method = 'loess', formula = y ~ x)
      + labs(title = "Distribution of z-scores between valid interacting chromosomes (significant interactions)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", ychr, " Genomic Position [Mb]"),
             y = "z-score",
             fill = "z-score sign")
      + scale_fill_manual(values =c("pos" = "#EE9B00", "neg" = "#005F73"), name = "z-score sign",labels = c("pos" = "Positive","neg" ="Negative"))
      + scale_colour_manual(values =c("pos" = "#EE9B00", "neg" = "#005F73"), name = "z-score sign", labels = c("pos" = "Positive","neg" ="Negative"))
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,ymax))
      + scale_x_continuous(expand = c(0, 0))
      #      + facet_wrap(.~ sign, scales = "free")
      #       + theme(panel.spacing = unit(0, "lines"),
      #               strip.text.y.right = element_text(angle = 0), #rotate facet labels
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_pts_posneg_YX_sigInters_all_cells",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()

########
# facet line and pts
########
#change positive to be listed first
hir_df$sign = factor(hir_df$sign, levels=c("pos" ="pos","neg" = "neg"))
p <- (ggplot(hir_df, aes(x=st1, y=zscore, fill=sign)) 
      + geom_point(alpha = 0.4, aes(color = sign))
      + geom_vline(aes(xintercept = xstart), colour = "red")
      + geom_smooth(aes(colour = sign), method = 'loess', formula = y ~ x)
      + labs(title = "Distribution of z-scores between valid interacting chromosomes (significant interactions)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", xchr, " Genomic Position [Mb]"),
             y = "z-score",
             fill = "z-score sign")
      + scale_fill_manual(values =c("pos" = "#EE9B00", "neg" = "#005F73"), name = "z-score sign",labels = c("pos" = "Positive","neg" ="Negative"))
      + scale_colour_manual(values =c("pos" = "#EE9B00", "neg" = "#005F73"), name = "z-score sign", labels = c("pos" = "Positive","neg" ="Negative"))
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,xmax))
      + scale_x_continuous(expand = c(0, 0))
            + facet_grid(sign ~ ., scales = "free")
             + theme(panel.spacing = unit(0, "lines"),
                     strip.text.y.right = element_blank())
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_pts_posneg_facet_sigInters_all_cells",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
p <- (ggplot(hir_df, aes(x=st2, y=zscore, fill=sign)) 
      + geom_point(alpha = 0.4, aes(color = sign))
      + geom_vline(aes(xintercept = ystart), colour = "red")
      + geom_smooth(aes(colour = sign), method = 'loess', formula = y ~ x)
      + labs(title = "Distribution of z-scores between valid interacting chromosomes (significant interactions)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", ychr, " Genomic Position [Mb]"),
             y = "z-score",
             fill = "z-score sign")
      + scale_fill_manual(values =c("pos" = "#EE9B00", "neg" = "#005F73"), name = "z-score sign",labels = c("pos" = "Positive","neg" ="Negative"))
      + scale_colour_manual(values =c("pos" = "#EE9B00", "neg" = "#005F73"), name = "z-score sign", labels = c("pos" = "Positive","neg" ="Negative"))
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,ymax))
      + scale_x_continuous(expand = c(0, 0))
            + facet_grid(sign ~ ., scales = "free")
             + theme(panel.spacing = unit(0, "lines"),
                     strip.text.y.right = element_blank())
      #      + facet_wrap(.~ sign, scales = "free")
      #       + theme(panel.spacing = unit(0, "lines"),
      #               strip.text.y.right = element_text(angle = 0), #rotate facet labels
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_pts_posneg_facet_YX_sigInters_all_cells",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()


#hex map
hir_df$sign <- factor(hir_df$sign,      # Reordering group factor levels
                         levels = c("pos", "neg"))
p <- (ggplot(hir_df, aes(y = zscore, fill = sign, x = st1))
      + geom_vline(aes(xintercept = xstart), colour = "red")
      + facet_wrap(~sign, ncol = 1, scales  = 'free_y') 
      + stat_binhex(aes(alpha = ..count..), colour = 'grey90')  
      + scale_alpha(name = 'Frequency', range = c(0,1))  
      + labs(title = "Distribution of z-scores between valid interacting chromosomes (significant interactions)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", xchr, " Genomic Position [Mb]"),
             y = "z-score",
             fill = "z-score sign")
      + theme(strip.background = element_blank(), strip.text = element_blank())
      + scale_fill_manual(values =c("pos" = "#EE9B00", "neg" = "#005F73"), name = "z-score sign",labels = c("pos" = "Positive","neg" ="Negative"))
      + expand_limits(x = c(0,xmax))
      + scale_x_continuous(expand = c(0, 0))
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_hex_posneg_sigInters_all_cells",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
p <- (ggplot(hir_df, aes(y = zscore, fill = sign, x = st2))
      + geom_vline(aes(xintercept = ystart), colour = "red")
  + facet_wrap(~sign, ncol = 1, scales  = 'free_y') 
  + stat_binhex(aes(alpha = ..count..), colour = 'grey80')  
  + scale_alpha(name = 'Frequency', range = c(0,1))  
  + labs(title = "Distribution of z-scores between valid interacting chromosomes (significant interactions)",
         #         subtitle = "Plot of length by dose",
         #         caption = "Data source: ToothGrowth",
         x = paste0("Chromosome ", ychr, " Genomic Position [Mb]"),
         y = "z-score",
         fill = "z-score sign")
  + theme(strip.background = element_blank(), strip.text = element_blank())
  + scale_fill_manual(values =c("pos" = "#EE9B00", "neg" = "#005F73"), name = "z-score sign",labels = c("pos" = "Positive","neg" ="Negative"))
      + expand_limits(x = c(0,ymax))
      + scale_x_continuous(expand = c(0, 0))
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_hex_posneg_YX_sigInters_all_cells",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()


print("# boxplot of above data")
p <- (ggplot(hir_df, aes(x=st1, y=zscore, group=st1, fill=factor(sign)) )
      #p <- (ggplot(hir_df, aes(x=st1, y=zscore)) 
      #      + geom_point(alpha = 0.4)
      + geom_boxplot()
      + geom_vline(aes(xintercept = xstart), colour = "red")
      + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", xchr, " Genomic Position [Mb]"),
             y = "z-score",
             fill = "z-score sign")
      #       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
      #       + facet_grid(cell ~.)
      + scale_fill_manual(values =c("pos" = "#EE9B00", "neg" = "#005F73"), name = "z-score sign",labels = c("pos" = "Positive","neg" ="Negative"))
      + scale_colour_manual(values =c("pos" = "#EE9B00", "neg" = "#005F73"), name = "z-score sign", labels = c("pos" = "Positive","neg" ="Negative"))
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,xmax))
      + scale_x_continuous(expand = c(0, 0))
      #       + theme(panel.spacing = unit(0, "lines"),
      #               strip.text.y.right = element_text(angle = 0), #rotate facet labels
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_box_posneg_sigInters_all_cells",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
p <- (ggplot(hir_df, aes(x=st2, y=zscore, group=st2, fill=sign)) 
      #p <- (ggplot(hir_df, aes(x=st1, y=zscore)) 
      #      + geom_point(alpha = 0.4)
      + geom_boxplot()
      + geom_vline(aes(xintercept = ystart), colour = "red")
      + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", ychr, " Genomic Position [Mb]"),
             y = "z-score",
             fill = "z-score sign")
      #       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
      #       + facet_grid(cell ~.)
      + scale_fill_manual(values =c("pos" = "#EE9B00", "neg" = "#005F73"), labels = c("pos" = "Positive","neg" ="Negative"))
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,ymax))
      + scale_x_continuous(expand = c(0, 0))
      #       + theme(panel.spacing = unit(0, "lines"),
      #               strip.text.y.right = element_text(angle = 0), #rotate facet labels
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_box_posneg_YX_sigInters_all_cells",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
########
# facet box plot
########
p <- (ggplot(hir_df, aes(x=st1, y=zscore, group=st1, fill=factor(sign)) )
      #p <- (ggplot(hir_df, aes(x=st1, y=zscore)) 
      #      + geom_point(alpha = 0.4)
      + geom_boxplot()
      + geom_vline(aes(xintercept = xstart), colour = "red")
      + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", xchr, " Genomic Position [Mb]"),
             y = "z-score",
             fill = "z-score sign")
      #       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
      #       + facet_grid(cell ~.)
      + scale_fill_manual(values =c("pos" = "#EE9B00", "neg" = "#005F73"), name = "z-score sign",labels = c("pos" = "Positive","neg" ="Negative"))
      + scale_colour_manual(values =c("pos" = "#EE9B00", "neg" = "#005F73"), name = "z-score sign", labels = c("pos" = "Positive","neg" ="Negative"))
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,xmax))
      + scale_x_continuous(expand = c(0, 0))
            + facet_grid(sign ~ ., scales = "free")
             + theme(panel.spacing = unit(0, "lines"),
                     strip.text.y.right = element_blank())
      #       + theme(panel.spacing = unit(0, "lines"),
      #               strip.text.y.right = element_text(angle = 0), #rotate facet labels
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_box_posneg_facet_sigInters_all_cells",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
p <- (ggplot(hir_df, aes(x=st2, y=zscore, group=st2, fill=sign)) 
      #p <- (ggplot(hir_df, aes(x=st1, y=zscore)) 
      #      + geom_point(alpha = 0.4)
      + geom_boxplot()
      + geom_vline(aes(xintercept = ystart), colour = "red")
      + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", ychr, " Genomic Position [Mb]"),
             y = "z-score",
             fill = "z-score sign")
      #       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
      #       + facet_grid(cell ~.)
      + scale_fill_manual(values =c("pos" = "#EE9B00", "neg" = "#005F73"), labels = c("pos" = "Positive","neg" ="Negative"))
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,ymax))
      + scale_x_continuous(expand = c(0, 0))
            + facet_grid(sign ~ ., scales = "free")
             + theme(panel.spacing = unit(0, "lines"),
                     strip.text.y.right = element_blank())
      #       + theme(panel.spacing = unit(0, "lines"),
      #               strip.text.y.right = element_text(angle = 0), #rotate facet labels
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_box_posneg_facet_YX_sigInters_all_cells",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()


p <- (ggplot(hir_df, aes(x=st1, y=zscore)) 
       + geom_point(alpha = 0.4)
       + geom_vline(aes(xintercept = xstart), colour = "red")
      + geom_smooth(colour = "black", method = 'loess', formula = y ~ x)
       + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions)",
              #         subtitle = "Plot of length by dose",
              #         caption = "Data source: ToothGrowth",
              x = paste0("Chromosome ", xchr, " Genomic Position [Mb]"),
              y = "z-score")
#       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
#       + facet_grid(cell ~.)
#       # expand axis limits so whole chrom len is accounted for
       + expand_limits(x = c(0,xmax))
#       + theme(panel.spacing = unit(0, "lines"),
#               strip.text.y.right = element_text(angle = 0), #rotate facet labels
#               strip.background = element_rect(fill = "white"),
#               #               axis.text.x = element_blank(),
#               #               axis.ticks.x = element_blank(),
#               axis.text.y = element_blank(),
#               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_pts_sigInters_all_cells",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
p <- (ggplot(aorta_df, aes(x=st1, y=zscore)) 
#p <- (ggplot(hir_df, aes(x=st1, y=zscore)) 
      + geom_point(alpha = 0.4)
      + geom_smooth(colour = "black", method = 'loess', formula = y ~ x)
+ geom_vline(aes(xintercept = xstart), colour = "red")
      + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions, Aorta)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", xchr, " Genomic Position [Mb]"),
             y = "z-score")
      #       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
      #       + facet_grid(cell ~.)
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,xmax))
      #       + theme(panel.spacing = unit(0, "lines"),
      #               strip.text.y.right = element_text(angle = 0), #rotate facet labels
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_pts_sigInters_Aorat_test",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
p <- (ggplot(hir_df, aes(x=st2, y=zscore)) 
      + geom_point(alpha = 0.4)
      + geom_smooth(colour = "black", method = 'loess', formula = y ~ x)
      + geom_vline(aes(xintercept = ystart), colour = "red")
      + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", ychr, " Genomic Position [Mb]"),
             y = "z-score")
      #       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
      #       + facet_grid(cell ~.)
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,ymax))
      #       + theme(panel.spacing = unit(0, "lines"),
      #               strip.text.y.right = element_text(angle = 0), #rotate facet labels
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_pts_sigInters_YX_all_cells",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
p <- (ggplot(aorta_df, aes(x=st2, y=zscore)) 
      #p <- (ggplot(hir_df, aes(x=st1, y=zscore)) 
      + geom_point(alpha = 0.4)
      + geom_smooth(colour = "black", method = 'loess', formula = y ~ x)
      + geom_vline(aes(xintercept = ystart), colour = "red")
      + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions, Aorta)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", ychr, " Genomic Position [Mb]"),
             y = "z-score")
      #       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
      #       + facet_grid(cell ~.)
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,ymax))
      #       + theme(panel.spacing = unit(0, "lines"),
      #               strip.text.y.right = element_text(angle = 0), #rotate facet labels
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_pts_sigInters_YX_Aorat_test",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
print("# just the line, no pts")
p <- (ggplot(hir_df, aes(x=st1, y=zscore)) 
#      + geom_point(alpha = 0.4)
      + geom_smooth(colour = "black", method = 'loess', formula = y ~ x)
+ geom_vline(aes(xintercept = xstart), colour = "red")
      + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", xchr, " Genomic Position [Mb]"),
             y = "z-score")
      #       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
      #       + facet_grid(cell ~.)
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,xmax))
      #       + theme(panel.spacing = unit(0, "lines"),
      #               strip.text.y.right = element_text(angle = 0), #rotate facet labels
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_line_sigInters_all_cells",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
p <- (ggplot(aorta_df, aes(x=st1, y=zscore)) 
      #p <- (ggplot(hir_df, aes(x=st1, y=zscore)) 
#      + geom_point(alpha = 0.4)
      + geom_smooth(colour = "black", method = 'loess', formula = y ~ x)
+ geom_vline(aes(xintercept = xstart), colour = "red")
      + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions, Aorta)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", xchr, " Genomic Position [Mb]"),
             y = "z-score")
      #       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
      #       + facet_grid(cell ~.)
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,xmax))
      #       + theme(panel.spacing = unit(0, "lines"),
      #               strip.text.y.right = element_text(angle = 0), #rotate facet labels
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_line_sigInters_Aorat_test",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
p <- (ggplot(hir_df, aes(x=st2, y=zscore)) 
#      + geom_point(alpha = 0.4)
      + geom_smooth(colour = "black", method = 'loess', formula = y ~ x)
+ geom_vline(aes(xintercept = ystart), colour = "red")
      + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", ychr, " Genomic Position [Mb]"),
             y = "z-score")
      #       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
      #       + facet_grid(cell ~.)
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,ymax))
      #       + theme(panel.spacing = unit(0, "lines"),
      #               strip.text.y.right = element_text(angle = 0), #rotate facet labels
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_line_sigInters_YX_all_cells",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
p <- (ggplot(aorta_df, aes(x=st2, y=zscore)) 
      #p <- (ggplot(hir_df, aes(x=st1, y=zscore)) 
#      + geom_point(alpha = 0.4)
      + geom_smooth(colour = "black", method = 'loess', formula = y ~ x)
+ geom_vline(aes(xintercept = ystart), colour = "red")
      + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions, Aorta)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", ychr, " Genomic Position [Mb]"),
             y = "z-score")
      #       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
      #       + facet_grid(cell ~.)
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,ymax))
      #       + theme(panel.spacing = unit(0, "lines"),
      #               strip.text.y.right = element_text(angle = 0), #rotate facet labels
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_line_sigInters_YX_Aorat_test",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()

print("# boxplot of above data")
p <- (ggplot(hir_df, aes(x=st1, y=zscore, group=st1)) 
      #p <- (ggplot(hir_df, aes(x=st1, y=zscore)) 
      #      + geom_point(alpha = 0.4)
      + geom_boxplot()
      + geom_vline(aes(xintercept = xstart), colour = "red")
      + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", xchr, " Genomic Position [Mb]"),
             y = "z-score")
      #       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
      #       + facet_grid(cell ~.)
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,xmax))
      #       + theme(panel.spacing = unit(0, "lines"),
      #               strip.text.y.right = element_text(angle = 0), #rotate facet labels
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_box_sigInters_all_cells",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
p <- (ggplot(aorta_df, aes(x=st1, y=zscore, group=st1)) 
      #p <- (ggplot(hir_df, aes(x=st1, y=zscore)) 
      #      + geom_point(alpha = 0.4)
      + geom_boxplot()
      + geom_vline(aes(xintercept = xstart), colour = "red")
      + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions, Aorta)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", xchr, " Genomic Position [Mb]"),
             y = "z-score")
      #       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
      #       + facet_grid(cell ~.)
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,xmax))
      #       + theme(panel.spacing = unit(0, "lines"),
      #               strip.text.y.right = element_text(angle = 0), #rotate facet labels
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_box_sigInters_Aorat_test",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
p <- (ggplot(hir_df, aes(x=st2, y=zscore, group=st2)) 
      #p <- (ggplot(hir_df, aes(x=st1, y=zscore)) 
      #      + geom_point(alpha = 0.4)
      + geom_boxplot()
      + geom_vline(aes(xintercept = ystart), colour = "red")
      + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", ychr, " Genomic Position [Mb]"),
             y = "z-score")
      #       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
      #       + facet_grid(cell ~.)
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,ymax))
      #       + theme(panel.spacing = unit(0, "lines"),
      #               strip.text.y.right = element_text(angle = 0), #rotate facet labels
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_box_sigInters_YX_all_cells",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
p <- (ggplot(aorta_df, aes(x=st2, y=zscore, group=st2)) 
      #p <- (ggplot(hir_df, aes(x=st1, y=zscore)) 
      #      + geom_point(alpha = 0.4)
      + geom_boxplot()
      + geom_vline(aes(xintercept = ystart), colour = "red")
      + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions, Aorta)",
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", ychr, " Genomic Position [Mb]"),
             y = "z-score")
      #       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
      #       + facet_grid(cell ~.)
      #       # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,ymax))
      #       + theme(panel.spacing = unit(0, "lines"),
      #               strip.text.y.right = element_text(angle = 0), #rotate facet labels
      #               strip.background = element_rect(fill = "white"),
      #               #               axis.text.x = element_blank(),
      #               #               axis.ticks.x = element_blank(),
      #               axis.text.y = element_blank(),
      #               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_box_sigInters_YX_Aorat_test",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()



print("##########")
print("# tickplot per chrom, each bin represented fact with cells")
print("# mean zscore per bin")
head(hm_df)
#count  each inter twice
tp_A <- hm_df %>% select(chrA,st1,end1,cell,zscore,pvalue) 
colnames(tp_A) <- c("chr","st","end","cell","zscore","pvalue")
tp_B <- hm_df %>% select(chrB,st2,end2,cell,zscore,pvalue) 
colnames(tp_B) <- c("chr","st","end","cell","zscore","pvalue")
#select only sig inters
tp_dat <- rbind(tp_A,tp_B) %>% filter(pvalue <=0.05)
tp_dat_sum = tp_dat %>% group_by(chr,st,cell) %>% dplyr::summarize(mzscore=mean(zscore, na.rm = TRUE))
tp_dat_sum$chr <- gsub("chr", "", tp_dat_sum$chr)
#tp_dat_sum <- tp_dat_sum %>% mutate(chr=factor(chr, levels=p_chr_ord))
#scale pts by 1Mb
#tp_dat_sum$st <- tp_dat_sum$st / 1000000
summary(tp_dat_sum)
head(tp_dat_sum)
for(i in unique(tp_dat_sum$chr)) {
  #i="22"
  interpos <- xstart
  if (xchr == i) {
    interpos = xstart
  } else {
    interpos = ystart
  }
  print("######")
  print(i)
  xmax <- chrInf$size[match(paste0("chr",i),chrInf$chrom)] / 1000000
  tpp <- tp_dat_sum %>% filter(chr == i)
  tp <- (ggplot(tpp, aes(st, y=1, fill = mzscore))
         + geom_tile(aes(fill = mzscore), width = 1, height = 1)
         + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score per bin", na.value = "grey")
         #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
         #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
         + labs(x = paste0("Chromosome ", i, " position [Mb]"),
                y = "",
                title = "Trans-chromosomal interactions (significant) z-scores")
         + facet_grid(cell ~ .)
         + geom_vline(xintercept = interpos, colour = "black")
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
  filename <- paste0("zscore_chrom",i,"_mean_tickplot_sig_Interactions.pdf")
  pdf(filename, width = 14, height = 8)
  print(tp)
  dev.off()
  # line plot of mean zscore of all cell types per bin
  lp = tpp %>% group_by(chr,st) %>% dplyr::summarize(mmzscore=mean(mzscore, na.rm = TRUE))
  tp <- (ggplot(lp, aes(st, y= mmzscore))
         + geom_point(alpha = 0.5, size= 3)
         +geom_smooth(method = 'loess',formula ='y ~ x')
         #+ scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score per bin", na.value = "grey")
         #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
         #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
         + labs(x = paste0("Chromosome ", i, " position [Mb]"),
                y = "Mean z-score per bin across cells",
                title = "Trans-chromosomal interactions (significant) z-scores (all cells)")
         #+ facet_grid(cell ~ .)
         + geom_vline(xintercept = interpos, colour = "red")
         #       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
         #         c("nonsig" = "Non-significant", "sig" = "Significant"))))
         #       + theme(axis.text.x = element_text(angle = 90))
         + expand_limits(x = c(0,xmax))
         + scale_y_continuous(expand = c(0, 0))
         + scale_x_continuous(expand = c(0, 0))
         + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
                 strip.background = element_rect(fill = "white"),
                 panel.spacing = unit(0, "lines"))
                 #axis.text.y = element_blank(),
                 #axis.ticks.y = element_blank())
         #+ theme(panel.background = element_rect(fill = "grey85", colour = NA))
  )
  filename <- paste0("zscore_chrom",i,"_mean_line_sig_Interactions.pdf")
  pdf(filename, width = 14, height = 8)
  print(tp)
  dev.off()
}#for
#df of mean zscore for bubble graph
bgz_dat <- tp_dat_sum %>% group_by(chr,st) %>% dplyr::summarize(mmzscore=mean(mzscore, na.rm = TRUE))
print("##########")
print("##########")
print("# tickplot per chrom, each bin represented fact with cells")
print("# number of sig inters per bin")
head(hm_df)
#count  each inter twice
tp_A <- hm_df %>% select(chrA,st1,end1,cell,zscore,pvalue) 
colnames(tp_A) <- c("chr","st","end","cell","zscore","pvalue")
tp_B <- hm_df %>% select(chrB,st2,end2,cell,zscore,pvalue) 
colnames(tp_B) <- c("chr","st","end","cell","zscore","pvalue")
#select only sig inters
tp_dat <- rbind(tp_A,tp_B) %>% filter(pvalue <=0.05)
tp_dat_sum = tp_dat %>% group_by(chr,st,cell) %>% dplyr::summarize(numSig=n())
tp_dat_sum$chr <- gsub("chr", "", tp_dat_sum$chr)
tp_dat_sum
#tp_dat_sum <- tp_dat_sum %>% mutate(chr=factor(chr, levels=p_chr_ord))
#scale pts by 1Mb
#tp_dat_sum$st <- tp_dat_sum$st / 1000000
summary(tp_dat_sum)
head(tp_dat_sum)
for(i in unique(tp_dat_sum$chr)) {
  #i="22"
  interpos <- xstart
  if (xchr == i) {
    interpos = xstart
  } else {
    interpos = ystart
  }
  xmax <- chrInf$size[match(paste0("chr",i),chrInf$chrom)] / 1000000
  tpp <- tp_dat_sum %>% filter(chr == i)
  tp <- (ggplot(tpp, aes(st, y=1, fill = numSig))
         + geom_tile(aes(fill = numSig), width = 1, height = 1)
         + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Total number of significant interactions", na.value = "grey")
         #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
         #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
         + labs(x = paste0("Chromosome ", i, " position [Mb]"),
                y = "",
                title = "Trans-chromosomal interactions (significant)")
         + facet_grid(cell ~ .)
         + geom_vline(xintercept = interpos, colour = "black")
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
  filename <- paste0("numSig_inters_chrom",i,"_tickplot_sig_Interactions.pdf")
  pdf(filename, width = 14, height = 8)
  print(tp)
  dev.off()
  # line plot of mean zscore of all cell types per bin
  lp = tpp %>% group_by(chr,st) %>% dplyr::summarize(mnumSig=mean(numSig, na.rm = TRUE))
  tp <- (ggplot(lp, aes(st, y= mnumSig))
         + geom_point(alpha = 0.5, size= 3)
         +geom_smooth(method = 'loess',formula ='y ~ x')
         #+ scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score per bin", na.value = "grey")
         #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
         #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
         + labs(x = paste0("Chromosome ", i, " position [Mb]"),
                y = "Mean number of significant interactions per bin across cells",
                title = "Trans-chromosomal significant interactions (all cells)")
         #+ facet_grid(cell ~ .)
         + geom_vline(xintercept = interpos, colour = "red")
         #       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
         #         c("nonsig" = "Non-significant", "sig" = "Significant"))))
         #       + theme(axis.text.x = element_text(angle = 90))
         + expand_limits(x = c(0,xmax))
         + scale_y_continuous(expand = c(0, 0))
         + scale_x_continuous(expand = c(0, 0))
         + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
                 strip.background = element_rect(fill = "white"),
                 panel.spacing = unit(0, "lines"))
         #axis.text.y = element_blank(),
         #axis.ticks.y = element_blank())
         #+ theme(panel.background = element_rect(fill = "grey85", colour = NA))
  )
  filename <- paste0("numSig_inters_",i,"_mean_line_sig_Interactions.pdf")
  pdf(filename, width = 14, height = 8)
  print(tp)
  dev.off()
}#for
#df of mean number of inters for bubble graph
bgi_dat <- tp_dat_sum %>% group_by(chr,st) %>% dplyr::summarize(mnumSig=mean(numSig, na.rm = TRUE))
head(bgi_dat)

print("########")
print("# bubble graph for mean zscore and mean number of sig inters per chrom per bin across cells")
print("########")
#combine dfs
bg_dat <- merge(bgz_dat, bgi_dat, by = c('chr','st'))
head(bg_dat)
summary(bg_dat)
for(i in unique(bg_dat$chr)) {
  #i="22"
  interpos <- xstart
  if (xchr == i) {
    interpos = xstart
  } else {
    interpos = ystart
  }
  xmax <- chrInf$size[match(paste0("chr",i),chrInf$chrom)] / 1000000
  tpp <- bg_dat %>% filter(chr == i)
bg <- (ggplot(tpp, aes(st, y= mmzscore, size=mnumSig))
       + geom_point(alpha = 0.5)
       #+geom_smooth(method = 'loess',formula ='y ~ x')
       #+ scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "Mean z-score per bin", na.value = "grey")
       #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
       #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
       + labs(x = paste0("Chromosome ", i, " position [Mb]"),
              y = "Mean z-score per bin across cells",
              title = "Trans-chromosomal significant interactions (all cells)")
       #+ facet_grid(cell ~ .)
       + geom_vline(xintercept = interpos, colour = "red")
       #       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
       #         c("nonsig" = "Non-significant", "sig" = "Significant"))))
       #       + theme(axis.text.x = element_text(angle = 90))
       + expand_limits(x = c(0,xmax))
       + scale_y_continuous(expand = c(0, 0))
       + scale_x_continuous(expand = c(0, 0))
       + scale_size("Mean number of significant interactions")
       + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
               strip.background = element_rect(fill = "white"),
               panel.spacing = unit(0, "lines"))
       #axis.text.y = element_blank(),
       #axis.ticks.y = element_blank())
       #+ theme(panel.background = element_rect(fill = "grey85", colour = NA))
)
filename <- paste0("numSig_inters_zscore_",i,"_mean_bubble_sig_Interactions.pdf")
pdf(filename, width = 14, height = 8)
print(bg)
dev.off()
}#for

print("##########")


##################
print("# point graph of all zscores on each valid inter chrom and genomic position")
##################
sig_pts <- hm_dfsig
summary(sig_pts)
p <- (ggplot(sig_pts, aes(x=st1, y=zscore)) 
       + geom_point(colour = "#55828B", alpha = 0.6)
       + geom_smooth(colour = "#013040", method = "loess", formula = y ~ x)
       + labs(title = paste0("Significant trans-chromosomal interactions between chromosomes ",xchr, " and ", ychr),
              #         subtitle = "Plot of length by dose",
              #         caption = "Data source: ToothGrowth",
              x = paste0("Chromosome ", xchr, " Genomic Position [Mb]"),
              y = "z-score")
#       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "p-value")
       + facet_grid(cell ~.)
       + scale_x_continuous(expand = c(0, 0))
       + scale_y_continuous(expand = c(0, 0))
       # expand axis limits so whole chrom len is accounted for
       + expand_limits(x = c(0,xmax))
       + theme(panel.spacing = unit(0, "lines"),
               strip.text.y.right = element_text(angle = 0), #rotate facet labels
               strip.background = element_rect(fill = "white"),
               #               axis.text.x = element_blank(),
               #               axis.ticks.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("zscore_along_",xchr,"_facet_",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()


p <- (ggplot(sig_pts, aes(x=st2, y=zscore)) 
      + geom_point(colour = "#55828B", alpha = 0.6)
      + geom_smooth(colour = "#013040", method = "loess", formula = y ~ x)
      + labs(title = paste0("Significant trans-chromosomal interactions between chromosomes ",xchr, " and ", ychr),
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = paste0("Chromosome ", ychr, " Genomic Position [Mb]"),
             y = "z-score")
      #       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "p-value")
      + facet_grid(cell ~.)
      + scale_x_continuous(expand = c(0, 0))
      + scale_y_continuous(expand = c(0, 0))
      # expand axis limits so whole chrom len is accounted for
      + expand_limits(x = c(0,ymax))
      + theme(panel.spacing = unit(0, "lines"),
              strip.text.y.right = element_text(angle = 0), #rotate facet labels
              strip.background = element_rect(fill = "white"),
              #               axis.text.x = element_blank(),
              #               axis.ticks.x = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("zscore_along_",ychr,"_facet_",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()
print("#zscore point plot of sig inters per cell")
ucells <- unique(sig_pts$cell)
xchr
ychr
for(i in ucells) {
  #i="Aorta"
  print(i)
  c_dat <- sig_pts %>% filter(cell == i)
  p <- (ggplot(c_dat, aes(x=st2, y=zscore)) 
        + geom_point(colour = "#55828B", alpha = 0.6)
       + geom_vline(aes(xintercept = ystart), colour = "red")
        + geom_smooth(colour = "#013040", method = "loess", formula = y ~ x)
        + labs(title = paste0("Significant trans-chromosomal interactions in ",i," between chromosomes ",xchr, " and ", ychr),
               #         subtitle = "Plot of length by dose",
               #         caption = "Data source: ToothGrowth",
               x = paste0("Chromosome ", ychr, " Genomic Position [Mb]"),
               y = "z-score")
        #       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "p-value")
        + scale_x_continuous(expand = c(0, 0))
        + scale_y_continuous(expand = c(0, 0))
        # expand axis limits so whole chrom len is accounted for
        + expand_limits(x = c(0,ymax))
  
  )
  f_name <- gsub(" ","",paste("zscore_along_",ychr,"_",i,"_",Atype,".pdf"))
  pdf(f_name, width = 14, height = 8)
  print(p)
  dev.off()
  summary(c_dat)
  p <- (ggplot(c_dat, aes(x=st1, y=zscore)) 
        + geom_point(colour = "#55828B", alpha = 0.6)
       + geom_vline(aes(xintercept = xstart), colour = "red")
        + geom_smooth(colour = "#013040", method = "loess", formula = y ~ x)
        + labs(title = paste0("Significant trans-chromosomal interactions in ",i," between chromosomes ",xchr, " and ", ychr),
               #         subtitle = "Plot of length by dose",
               #         caption = "Data source: ToothGrowth",
               x = paste0("Chromosome ", xchr, " Genomic Position [Mb]"),
               y = "z-score")
        #       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "p-value")
        + scale_x_continuous(expand = c(0, 0))
        + scale_y_continuous(expand = c(0, 0))
        # expand axis limits so whole chrom len is accounted for
        + expand_limits(x = c(0,xmax))
  )
  f_name <- gsub(" ","",paste("zscore_along_",xchr,"_",i,"_",Atype,".pdf"))
  pdf(f_name, width = 14, height = 8)
  print(p)
  dev.off()
}#for

# heatmap with grey in place of all white space
##############
# heatmap along chromosome positions: SIG INTERACTIONS (only sig)
##############
hm_dfsig <- hm_df %>% filter(pvalue <= 0.05)
#hm <- (ggplot(testdf, aes(x=st1, y=st2, fill= zscore)) 
hm <- (ggplot(hm_dfsig, aes(x=st1, y=st2, fill= zscore)) 
       + geom_tile(width = 1, height = 1)
       + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions)",
              #         subtitle = "Plot of length by dose",
              #         caption = "Data source: ToothGrowth",
              x = paste0("Chromosome ", xchr, " Genomic Position [Mb]"),
              y = paste0("Chromosome ", ychr, " Genomic Position"))
       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
       + geom_vline(xintercept = xstart, colour = "black")
       + facet_grid(cell ~.)
       # expand axis limits so whole chrom len is accounted for
       + expand_limits(y=c(0,ymax), x = c(0,xmax))
      + scale_y_continuous(expand = c(0, 0))
      + scale_x_continuous(expand = c(0, 0))
#       +theme_grey(base_size=10)
       + theme(panel.spacing = unit(0, "lines"),
               strip.text.y.right = element_text(angle = 0), #rotate facet labels
               strip.background = element_rect(fill = "white"),
              panel.background = element_rect(fill = "grey85", colour = NA),
               #               axis.text.x = element_blank(),
               #               axis.ticks.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_heatmap_sigInters_GREYBG",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
hm
dev.off()
#having the opposite chrom on the x axis
hm <- (ggplot(hm_dfsig, aes(x=st2, y=st1, fill= zscore)) 
       + geom_tile(width = 1, height = 1)
       + labs(title = "Distribution of z-scores along valid interacting chromosome (significant interactions)",
              #         subtitle = "Plot of length by dose",
              #         caption = "Data source: ToothGrowth",
              x = paste0("Chromosome ", ychr, " Genomic Position [Mb]"),
              y = paste0("Chromosome ", xchr, " Genomic Position"))
       + geom_vline(xintercept = ystart, colour = "black")
       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "z-score")
       + facet_grid(cell ~.)
       # expand axis limits so whole chrom len is accounted for
       + expand_limits(y=c(0,xmax), x = c(0,ymax))
      + scale_y_continuous(expand = c(0, 0))
      + scale_x_continuous(expand = c(0, 0))
       + theme(panel.spacing = unit(0, "lines"),
               strip.text.y.right = element_text(angle = 0), #rotate facet labels
               strip.background = element_rect(fill = "white"),
              panel.background = element_rect(fill = "grey85", colour = NA),
               #               axis.text.x = element_blank(),
               #               axis.ticks.x = element_blank(),
               axis.text.y = element_blank(),
               axis.ticks.y = element_blank())
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_chroms_zscore_heatmap_sigInters_YX_GREYBG",Atype,".pdf"))
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
       + geom_tile(width = 1, height = 1)
       + labs(title = "Distribution of p-value along valid interacting chromosome (all interactions)",
              #         subtitle = "Plot of length by dose",
              #         caption = "Data source: ToothGrowth",
              x = paste0("Chromosome ", xchr, " Genomic Position [Mb]"),
              y = paste0("Chromosome ", ychr, " Genomic Position"))
       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "p-value")
       + facet_grid(cell ~.)
       + geom_vline(xintercept = xstart, colour = "black")
      # expand axis limits so whole chrom len is accounted for
      + expand_limits(y=c(0,ymax), x = c(0,xmax))
      + scale_y_continuous(expand = c(0, 0))
      + scale_x_continuous(expand = c(0, 0))
       + theme(panel.spacing = unit(0, "lines"),
               strip.text.y.right = element_text(angle = 0), #rotate facet labels
               strip.background = element_rect(fill = "white"),
#               axis.text.x = element_blank(),
#               axis.ticks.x = element_blank(),
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
       + geom_tile(width = 1, height = 1)
       + labs(title = "Distribution of p-value along valid interacting chromosome (significant interactions)",
              #         subtitle = "Plot of length by dose",
              #         caption = "Data source: ToothGrowth",
              x = paste0("Chromosome ", xchr, " Genomic Position [Mb]"),
              y = paste0("Chromosome ", ychr, " Genomic Position"))
       + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "p-value")
       + geom_vline(xintercept = xstart, colour = "black")
       + facet_grid(cell ~.)
      # expand axis limits so whole chrom len is accounted for
      + expand_limits(y=c(0,ymax), x = c(0,xmax))
      + scale_y_continuous(expand = c(0, 0))
      + scale_x_continuous(expand = c(0, 0))
       + theme(panel.spacing = unit(0, "lines"),
               strip.text.y.right = element_text(angle = 0), #rotate facet labels
               strip.background = element_rect(fill = "white"),
#               axis.text.x = element_blank(),
#               axis.ticks.x = element_blank(),
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

