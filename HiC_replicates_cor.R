########################################
# correlation between replicates for HiC (raw normalized data from cooler)
# used as justification for pooling replicates
######
# Developer: Daniella F. Lato
#            email:  daniellalato@gmail.com
#            github: https://github.com/dlato
######
# arguments: cooler output with all replicates. tsv similar to 3Dflow output
#            outfile prefix
########################################

options(echo=F)
options(scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
outprefix <- args[2]

##########
library(tidyr)
library(dplyr)
library(ggplot2)
.libPaths("/hpf/largeprojects/pmaass/programs/Rlib/R.4.0.2")
library(ggcorrplot)
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
###interaction data
#options(scipen = 999)
#dat_file <- "HiC_replicates_subset.txt"

dat <- read.table(dat_file, header = TRUE)
print("summary of cooler normalized inter freq per cell type")
summary(dat)

#compute correlation matrix
corr <- round(cor(na.omit(dat %>% dplyr::select(-ID))))
corr
write.table(corr, file = as.character(paste0(outprefix,"_correlation_coef.txt")), sep = "\t", quote = FALSE, row.names = T, col.names = T)
#compute matrix of cor pvals
p.mat <- cor_pmat(dat %>% dplyr::select(-ID))
p.mat
write.table(p.mat, file = as.character(paste0(outprefix,"_correlation_pvalues.txt")), sep = "\t", quote = FALSE, row.names = T, col.names = T)

#plot correlation matrix
p <- ggcorrplot(corr,
           hc.order = TRUE, type = "upper",
           outline.color = "white",
           lab = TRUE,
           p.mat = p.mat,
           insig = "blank",
           colors = c("#264653", "white", "#e76f51")
)

pdf(paste0(outprefix,"_correlation_between_replicates.pdf"), width = 10, height = 10)
print(p)
dev.off()
print("DONE")
#
