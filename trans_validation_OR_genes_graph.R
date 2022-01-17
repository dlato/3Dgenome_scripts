########################################
# number of OR genes per bin per chromosome
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: gene annotation for genes of interest (subset of gtf file with only genes of interest on one chromosome)
#            bin size (bp)
########################################

options(echo=F)
args <- commandArgs(trailingOnly = TRUE)
gtf_file <- args[1]
bin_size<- args[2]

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
##interaction data
#gtf_file <- "OR_gene_chr14.gtf"
#bin_size <- 1000000
#library(harrypotter)
gtf <- read.table(gtf_file, header = FALSE, sep = "\t")
colnames(gtf) <- c("chrom", "db", "geneType","st","end","score","strand","frame","otherInfo")
head(gtf)
summary(gtf)
# get bin start for each gene
gtf$bin_start <- plyr::round_any(gtf$st, bin_size, f = floor)
# sum up the number of genes per bin
num_genes <- gtf %>% group_by(bin_start) %>% dplyr::summarize(n=n())
# adjust bin start for graph
num_genes$scalestart <- num_genes$bin_start / 1000000
head(num_genes)
# chrom info
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

print("##########")
print("# tickplot with number of OR genes per bin")
xchr <- gsub("chr","", gtf$chrom[1])
xmax <- chrInf$size[match(paste0("chr",xchr),chrInf$chrom)] / 1000000
# adding values of 0 for all other bins
allbins = seq(0, xmax*1000000, 1000000)
allbins
# Filter out bins that are already present
bins0 = allbins[!(allbins %in% num_genes$bin_start)]
bin0 = data.frame(bin_start = bins0, n = 0, scalestart = bins0/1000000)

# Append this `data.frame` and resort in time:
num_genes2 = rbind(num_genes, bin0)
num_genes2 = num_genes2[order(num_genes2$scalestart),]

tp <- (ggplot(num_genes2, aes(scalestart, y=1, fill = n))
         + geom_tile(aes(fill = n), width = 1, height = 1)
         + scale_fill_hp(discrete = FALSE, option = "ronweasley2", name = "OR genes per bin", na.value = 0)
         + labs(x = paste0("Chromosome ", xchr, " position [Mb]"),
                y = "",
                title = "Number of Olfactory Receptor (OR) genes")
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
  filename <- paste0("num_OR_genes_chrom",xchr,"_tickplot.pdf")
  pdf(filename, width = 14, height = 2)
  print(tp)
  dev.off()
  