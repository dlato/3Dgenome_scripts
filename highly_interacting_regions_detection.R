########################################
# determining consecutive genomic regions/bins that have high mean z-score and high # of interactions along a chromosome
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: 3Dflow zscore output data (ONLY SIGNIFICANT INTERACTIONS)
#            percentile for determining highly interacting regions. If you want the top 10% of z-scores and number of interactions to be considered high, choose 0.9 (numeric)
#            bin size used for Hi-C data processing (numeric, in bp)
#
########################################

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
perc <- args[2]
bin_size <- args[3]

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


# for testing only!
library(harrypotter)
dat_file <- ""
perc <- 0.4
bin_size <- 1000000
chrom_in <- "12"
options(scipen = 999) 
# create fake df for now
bin_start <- seq(0,10000000,1000000)
zscore <- c(1.2,1.5,2,-0.2,-1,0,1.6,1.4,1,2,1.5)
numInters <- c(17,13,20,4,5,2,19,20,5,17,18)
tpp <- as.data.frame(cbind(bin_start,zscore,numInters))


print("#read in files")
tpp <- read.table(dat_file, header = TRUE)
head(tpp)
bg <- (ggplot(tpp, aes(bin_start, y= zscore, size=numInters, colour=numInters))
       + geom_point(alpha = 0.5)
       #+geom_smooth(method = 'loess',formula ='y ~ x')
       + scale_color_hp(discrete = FALSE, option = "ronweasley2", guide="legend")
       #       + scale_fill_hp_d(option = "Always", name = "Mean z-score") 
       #+ scale_fill_gradient(low = "white", high = "steelblue", name = "Mean z-score")
       + labs(x = paste0("Chromosome 12 position [Mb]"),
              y = "Mean z-score per bin across cells",
              title = "Trans-chromosomal significant interactions (all cells)",
              size = "Mean number of significant interactions per bin",
              color = "Mean number of significant interactions per bin")
       #+ facet_grid(cell ~ .)
       #+ geom_vline(xintercept = interpos, colour = "red")
       #       + facet_wrap(.~sig, labeller = labeller(sig= as_labeller(
       #         c("nonsig" = "Non-significant", "sig" = "Significant"))))
       #       + theme(axis.text.x = element_text(angle = 90))
       #+ expand_limits(x = c(0,xmax))
       + scale_y_continuous(expand = c(0, 0))
       + scale_x_continuous(expand = c(0, 0))
       #+ scale_size("Mean number of significant interactions")
       + theme(strip.text.y.right = element_text(angle = 0), #rotate facet labels
               strip.background = element_rect(fill = "white"),
               panel.spacing = unit(0, "lines"))
       #axis.text.y = element_blank(),
       #axis.ticks.y = element_blank())
       #+ theme(panel.background = element_rect(fill = "grey85", colour = NA))
)
bg

#get userspecified nth percentile
toppercZ <- quantile(tpp$zscore, probs = perc)
toppercZ
toppercN <- quantile(tpp$numInters, probs = perc)
toppercN
#find regions that are above percent cutoff
highinter <- tpp %>% arrange(bin_start) %>%
  filter(zscore >= toppercZ) %>%
  filter(numInters >= toppercN)
highinter
# go through df and create bed file with the regions
# skipping up to 2 consecutive regions is ok
highinter$tdiff <- c(NA,diff(highinter$bin_start))
highinter
bedstart <- c()
bedend <- c()
tmps <- NA
tmps <- NA
for (i in 1:nrow(highinter)) {
  if (is.na(highinter$tdiff[i])){
    tmps <- highinter$bin_start[i]
  } else {
  # end of the df
  if (i == nrow(highinter)){
      tmpe <- highinter$bin_start[i] + bin_size
      bedstart <- c(bedstart,tmps)
      bedend <- c(bedend,tmpe)
  }#if
    # if difference is less than 2 bins, keep going until we get a diff that is >2 bins
    if (highinter$tdiff[i] <= bin_size *2) {
      next
    # difference is more than 2 bins
    } else {
      tmpe <- highinter$bin_start[i-1] + bin_size
      bedstart <- c(bedstart,tmps)
      bedend <- c(bedend,tmpe)
      tmps <- highinter$bin_start[i]
    }# if else
  }#if else
}#for
bedstart
bedend
chrom <- rep(paste0("chr",chrom_in),length(bedstart))
bed_df <- as.data.frame(cbind(chrom,bedstart,bedend))
bed_df
#save bed file to table
write.table(bed_df, file = as.character(paste0("highly_interacting_regions_",perc,"_percentile.bed")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
highinter