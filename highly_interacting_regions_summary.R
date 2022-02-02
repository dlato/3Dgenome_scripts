########################################
#summary plots/stats for the automatically generated highly interacting regions
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: bed-like file for highly interacing regions on CHROM PAIRS (tsv, (chr, start, end, chrpair)))
#            bed-like file for highly interacting regions on CHROM VS ALL OTHER INTERACTIONS (tsv, (chr, start, end, chrpair)))
########################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
pairs_file <- args[1]
onevall_file <- args[2]

##########
library(tidyr)
library(dplyr)
library(ggplot2)
.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2")
library(karyoploteR)#for karyotype plot
#library(karyoploteR)#for karyotype plot
#library(BRGenomics)#for karyotype plot
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
                  legend.text = element_text(size = 16),
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


#########################################################################
print("#read in files")
####interaction data
##library(circlize) # for circos
#options(scipen = 999)
#pairs_file <- "highly_interacting_regions_chrPairs_0.90_percentile.bed"
#onevall_file <- "highly_interacting_regions_all_chroms_0.90_percentile.bed"

pairs_df <- read.table(pairs_file, header = FALSE)
colnames(pairs_df) <- c("chr","start","end","chrPair")
head(pairs_df)

onevall_df <- read.table(onevall_file, header = FALSE)
colnames(onevall_df) <- c("chr","start","end","chrPair")
head(onevall_df)

#combine above files into one df
pairs_df$class <- rep("pairs",nrow(pairs_df))
onevall_df$class <- rep("onevall",nrow(onevall_df))
dat <- rbind(pairs_df,onevall_df)
#get length of highly interacting region
dat$len <- dat$end - dat$start
head(dat)
summary(dat)

print("# number of highly interacting regions for each analysis (pairs and 1vsAll)")
dat %>% group_by(class) %>% dplyr::summarize(nSig=n())
dat$len <- dat$len / 1000000
dat$chr <- sub("chr","", dat$chr)


print("#############")
print("#check if means are different between pairwise and chrA vs ALL groups")
print("#Test each group for normality")
print("sig = reject normality null")
#chrClass_dat$chrClass <- as.factor(chrClass_dat$chrClass)
dat %>%
  group_by(class) %>%
  summarise(W = shapiro.test(len)$statistic,
            p.value = shapiro.test(len)$p.value)
print("#Perform the Mann-Whitney (non-parametric) test because not normal")
print("sig = mean is diff btwn groups")
wilcox.test(len ~ class, data=dat)



#change order of class levels
dat <- dat %>%
  mutate(class = factor(class, levels=c("pairs", "onevall")))  # This trick update the factor levels

#histogram of highly interacting region sizes
bp <- (ggplot(dat, aes(x=len,fill=class))
       +geom_histogram(position = "dodge", color = "black")
       + scale_fill_manual(values = c("#FAC9A1", "#013040"), labels= c("Pairwise", "ChrA vs All"))
       #  + scale_fill_hp(discrete = FALSE, option = "Always", name = "Mean z-score per chromosomal pair", na.value = "grey")
       #  #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
       #  #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
       + labs(y = "Number of Highly Interacting Regions",
              x = "Length of Highly Interacting Region [Mb]",
              title = "",
              fill = "Type of Analysis")
       + scale_y_continuous(expand = c(0,0))
       + scale_x_continuous(expand = c(0,0))
       + theme(strip.text.y.right = element_text(angle = 0))
)
pdf("hightly_interacting_regions_histogram_length.pdf", width = 14, height = 8)
bp
dev.off()


#boxplot/violin of region sizes
set.seed(369)
levels(dat$class) <- c("Pairwise",  # Relevel factor labels
                     "ChrA vs All")
p <- (ggplot(dat, aes(y=len,x=class, fill = class, color = class))
      +geom_violin(alpha = 0.7, color = "black")
      + geom_boxplot(width = 0.25, alpha = 0.7, color = "black", outlier.colour = NA)
      + geom_jitter(alpha = 0.3, color = "black")
       + scale_fill_manual(values = c("#FAC9A1", "#013040"), labels= c("Pairwise", "ChrA vs All"))
       + scale_color_manual(values = c("#FAC9A1", "#013040"), labels= c("Pairwise", "ChrA vs All"))
      + labs(y="Length of Highly Interacting Region [Mb]",
             x="Type of Analysis",
             title = "",
             fill = "Type of Analysis")
)
pdf("hightly_interacting_regions_box_violin_length.pdf", width = 14, height = 8)
p
dev.off()

#boxplot of region sizes per chrom
p <- (ggplot(dat, aes(y=len,x=chr, fill = class, color = class))
      + geom_boxplot(alpha = 0.9, color = "black", outlier.colour = NA)
#      + geom_jitter(width = 0.1, alpha = 0.4, color = "black")
      + scale_fill_manual(values = c("#FAC9A1", "#013040"), labels= c("Pairwise", "ChrA vs All"))
      + scale_color_manual(values = c("#FAC9A1", "#013040"), labels= c("Pairwise", "ChrA vs All"))
      + labs(y="Length of Highly Interacting Region [Mb]",
             x="Type of Analysis",
             title = "",
             fill = "Chromosome")
)
pdf("hightly_interacting_regions_boxplot_per_chromosome_length.pdf", width = 14, height = 8)
p
dev.off()

###############
# karyotype of regions
##############
dat$chr <- sub("^","chr", dat$chr)
head(dat)
kdat <- toGRanges(dat %>% select(chr,start,end, len))
head(kdat)
pdf("highly_interacting_regions_karyotype.pdf", width = 14, height = 8)
kp <- plotKaryotype(genome = "hg38")
kpAddBaseNumbers(kp)
kpHeatmap(kp, data=kdat,y=kdat$len, colors = c("#FFAA00","#5F0B32"),r0=0.05, r1=0.5)
dev.off()

print("DONE")