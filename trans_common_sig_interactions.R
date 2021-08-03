########################################
#validating known trans-chromosomal interactions
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: 3Dflow output data
########################################

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
Atype <- args[4]

##########
library(dplyr)
library(tidyr)
#library(GenomicRanges)
library(ggplot2)
library(ggforce)
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
                  legend.title = element_blank(),
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
dat <- read.table("test_pairwise_dat.txt", header = TRUE)
dat <- read.table(dat_file, header = TRUE)
dat <- as.data.frame(dat)
print("summary of ALL sig zscores per cell type")
summary(dat)
dat$ID <- as.character(dat$ID)

#remove interactions involving x and y chrs
#dat <- dat[grep("chrY", df$ID, invert=TRUE), ]
#dat2 <- dat[grep("chrY", df$ID), ]
#dat <- dat[grep("chrX", df$ID, invert=TRUE), ]

#select only rows with NO NAs in any cell type
df <- na.omit(dat)
print("percent of ALL common interactions across genome")
(length(df$ID)/length(dat$ID)) *100


#prep data for parallel sets plot
#split ID col
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
df$ID <- sub("B", "\\.B", as.character(df$ID))
dat2 <- df %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
#remove A and B from chrom names
dat2$chrA <- gsub("A", "", dat2$chrA)
dat2$chrB <- gsub("B", "", dat2$chrB)
ps_df <- dat2 %>% select(chrA, chrB)
ps_df <- unique(ps_df)
ps_df$fake <- rep(1,length(ps_df$chrA))
ps_df<- ps_df %>%
  gather_set_data(1:3)
ps_df
#plot
ps <- (ggplot(data =ps_df, aes(chrA, id=id, split = chrB, value = 1))
       #  + geom_parallel_sets(aes(fill = U00096000))
       + geom_parallel_sets()
#       + scale_fill_manual(values = c("#2E294E","#BEBEBE"))
       #  + geom_parallel_sets(aes(fill = U00096 ))
       + xlab("test") 
       + ylab("Genomic Position")
       + coord_flip()
#       + scale_x_discrete(expand = c(0,0))
#       + theme(legend.title=element_blank())
)
ps


data <- reshape2::melt(Titanic)
head(data)
data <- gather_set_data(data, 1:4)
tail(data)






#roi1
#roi1 <- read.table("CISTR.bed", header = FALSE)
#roi1 <- read.table("FIRRE.bed", header = FALSE)
roi1 <- read.table(roi1_file)
colnames(roi1) <- c("chrom","start","end")
#roi2
#roi2 <- read.table("SOX9.bed", header = FALSE)
#roi2 <- read.table("ATF4.bed", header = FALSE)
roi2 <- read.table(roi2_file)
colnames(roi2) <- c("chrom","start","end")

print("#new cols for interaction ID")
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
dat$ID <- sub("B", "\\.B", as.character(dat$ID))
head(dat)
dat2 <- dat %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
#remove A and B from chrom names
dat2$chrA <- gsub("A", "", dat2$chrA)
dat2$chrB <- gsub("B", "", dat2$chrB)
head(dat2)
#convert data to GRanges object
dat_bed1 <- dat2 %>% select(chrA, st1, end1)
dat_bed1$st1 <- as.numeric(dat_bed1$st1)
dat_bed1$end1 <- as.numeric(dat_bed1$end1)
summary(dat_bed1)
dat_chr1 <- makeGRangesFromDataFrame(dat_bed1, seqnames.field=c("chrA"),
                                     start.field="st1",
                                     end.field=c("end1"))
dat_bed2 <- dat2 %>% select(chrA, st2, end2)
dat_bed2$st2 <- as.numeric(dat_bed2$st2)
dat_bed2$end2 <- as.numeric(dat_bed2$end2)
summary(dat_bed2)
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
roi1_res
print("#resolution with roi2")
roi2_inter1 <- unique(subsetByOverlaps(dat_chr1, roi2_G))
roi2_inter2 <- unique(subsetByOverlaps(dat_chr2, roi2_G))
roi2_inter1$overl<- countOverlaps(roi2_inter1, roi2_inter2, type="equal")
roi2_inter1 <- as.data.frame(roi2_inter1)
roi2_inter1
roi2_res <- paste0(roi2_inter1$seqnames,".",roi2_inter1$start,".", roi2_inter1$end)
roi2_res
print("#data with JUST interacting region")
Vinter <- dat2[grep(roi1_res, dat2$ID), ]
Vinter <- Vinter[grep(roi2_res, Vinter$ID),]
#Vinter[grep(roi2_res, Vinter$ID),]
#test_D<- dat2[grep(roi2_res, dat2$ID),]
#test_D[grep(roi1_res, test_D$ID),]
#Vinter
#dim(Vinter)
chrs <- c(levels(roi2$chrom),levels(roi1$chrom[1]))
print("#######")
print("# Validated interaction specific region")
print("#######")
chrs_df <- Vinter[grep(chrs[1], Vinter$ID), ]
chrs_df <- chrs_df[grep(chrs[2],chrs_df$ID),]
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
p <- (ggplot(pa_dat, aes(y=all_cells, x=allmisscols, shape =allmisscols))
      + geom_point(size=4)
#      + coord_flip()
      + labs(title = mytitle,
             #         subtitle = "Plot of length by dose",
             #         caption = "Data source: ToothGrowth",
             x = "Valid interaction status", y = "")
      #         tag = "A")
      + scale_shape_manual(values=c(16,1))
      + scale_fill_manual(values =c("plum4", "cadetblue"))
      #  + facet_grid(seqDep ~ cell )
      + theme(legend.position = "none")
)
f_name <- gsub(" ","",paste("allCells_valid_interaction_present_absent_",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()


print("#######")
print("# Validated interaction chromosomes (not specific region)")
print("#######")
#df with ALL sig interactions btwn validated chroms
chrs_df <- dat2[grep(chrs[1], dat2$ID), ]
chrs_df <- chrs_df[grep(chrs[2],chrs_df$ID),]
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
sixcells_all <- select(dat2,c("ID",
                            "Dorsolateral_Prefrontal_cortex", 
                            "Small_bowell", 
                            "Aorta", 
                            "Right_Ventricle", 
                            "Cardiomyocites_primitive_rep1",
                            "Spleen"))
head(sixcells_all)
sixcells_Vinter <- select(Vinter,c("ID",
                                 "Dorsolateral_Prefrontal_cortex", 
                                 "Small_bowell", 
                                 "Aorta", 
                                 "Right_Ventricle", 
                                 "Cardiomyocites_primitive_rep1",
                                 "Spleen"))

head(sixcells_Vinter)
sixcells_chr <- select(chrs_df,c("ID",
                                 "Dorsolateral_Prefrontal_cortex", 
                                 "Small_bowell", 
                                 "Aorta", 
                                 "Right_Ventricle", 
                                 "Cardiomyocites_primitive_rep1",
                                 "Spleen"))
head(sixcells_chr)

#wide to long
#all data (six cells)
all_long <- gather(sixcells_all, cell, zscore, 2:length(colnames(sixcells_all)), factor_key=TRUE)
all_long$validType <- rep("all", length(all_long$ID))
summary(all_long)

#valid inter (six cells)
valid_inter_long <- gather(sixcells_Vinter, cell, zscore, 2:length(colnames(sixcells_Vinter)), factor_key=TRUE)
valid_inter_long$validType <- rep("validInter", length(valid_inter_long$ID))
summary(valid_inter_long)

#valid chrs (six cells)
valid_chrs_long <- gather(sixcells_chr, cell, zscore, 2:length(colnames(sixcells_chr)), factor_key=TRUE)
valid_chrs_long$validType <- rep("validChr", length(valid_chrs_long$ID))
summary(valid_chrs_long)

#combine the dfs
plot_df <- rbind(all_long, valid_inter_long, valid_chrs_long)
summary(plot_df)

# rename levels
plot_df$validType <- as.factor(plot_df$validType)
plot_df$cell <- as.factor(plot_df$cell)
levels(plot_df$validType) <- list("Interactions from valid inter" = "validInter",
                                  "Interactions from chrs"="validChr",
                                  "All significant interactions"="all")
levels(plot_df$cell) <- list("Spleen" = "Spleen",
                                  "CardPrimRep1" = "Cardiomyocites_primitive_rep1",
                                  "RightVentricle" = "Right_Ventricle",
                                  "Aorta"="Aorta",
                                  "SmallBowel"="Small_bowell",
                                  "DorsoPreCort" = "Dorsolateral_Prefrontal_cortex")
#new col for seq depth
plot_df$seqDep <- plot_df$cell
levels(plot_df$seqDep) <- list("High" = "Spleen",
                                  "High" = "CardPrimRep1",
                                  "Medium" = "RightVentricle",
                                  "Medium"="Aorta",
                                  "Low" = "SmallBowel",
                                  "Low" = "DorsoPreCort")
head(plot_df)

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
target <- c("Interactions from chrs", "All significant interactions")
plot_d <- plot_df %>% filter(validType %in% target)
p <- (ggplot(plot_d, aes(cell, zscore, fill = validType))
  + geom_split_violin()
  + coord_flip()
  + labs(title = mytitle,
#         subtitle = "Plot of length by dose",
#         caption = "Data source: ToothGrowth",
         x = "", y = "z-score")
#         tag = "A")
  + scale_fill_manual(values =c("plum4", "cadetblue"))
#  + facet_grid(seqDep ~ cell )
)
#p
f_name <- gsub(" ","",paste("6testCells_valid_interaction_chroms_density_",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
p
dev.off()

#dealing with cases where the specified interaction is not sig in any cell type
if (dim(Vinter)[1] == 0){
  print("dataframe is empty :( your specified interaction was not significant in any of your cell types/data")
} else {
  print("sig interactions in at least ONE cell type!")
  #######
  # Validated interaction specific region
  #######
  #plot above data for valid chroms
plot_d <- plot_df %>% filter(validType == "All significant interactions")
inter_d <- plot_df %>% filter(validType == "Interactions from valid inter")
summary(inter_d)
plot_d[which(plot_d$validType == "Interactions from valid inter"),]
  p <- (ggplot(plot_d, aes(x=zscore))
        + geom_density(fill = "cadetblue")
#        + coord_flip()
        + geom_vline(data = inter_d, aes(xintercept = zscore, 
                                               color = cell), size=1.5)
        + labs(title = mytitle,
               #         subtitle = "Plot of length by dose",
               #         caption = "Data source: ToothGrowth",
               x = "z-score", y = "density")
        #         tag = "A")
#        + scale_color_manual(values =c("plum4"))
        #  + facet_grid(seqDep ~ cell )
  )
#  p
f_name <- gsub(" ","",paste("6testCells_valid_interaction_density_",Atype,".pdf"))
pdf(f_name, width = 14, height = 8)
print(p)
dev.off()
}


