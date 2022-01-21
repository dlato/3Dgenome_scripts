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
options(scipen = 999) 
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
perc <- args[2]
bin_size <- args[3]

##########
library(tidyr)
library(dplyr)
##########
##########################################################################

# for testing only!
#library(harrypotter)
#dat_file <- "test_1vsAll_dat.txt"
#perc <- 0.4
#bin_size <- 1000000
#chrom_in <- "12"
### create fake df for now
##bin_start <- seq(0,10000000,1000000)
##zscore <- c(1.2,1.5,2,-0.2,-1,0,1.6,1.4,1,2,1.5)
##numInters <- c(17,13,20,4,5,2,19,20,5,17,18)
##tpp <- as.data.frame(cbind(bin_start,zscore,numInters))
## reading in chr 12 chr 17 sig zscore data from 3Dflow run
#dat_file <- "chr12chr17_zscores.txt"


print("#read in files")
dat <- read.table(dat_file, header = TRUE)
#split up ID col
colnm <- c("chrA", "st1", "end1","chrB","st2","end2")
dat$ID <- sub("B", "\\.B", as.character(dat$ID))
dat <- dat %>% separate(ID, sep = "\\.", into = colnm, remove = FALSE)
#remove A and B from chrom names
dat$chrA <- gsub("A", "", dat$chrA)
dat$chrB <- gsub("B", "", dat$chrB)
#long format for the cell type
dat <- gather(dat, cell, zscore, 8:ncol(dat), factor_key=TRUE)
#col for chrom pair
dat$chrPair <- paste0(dat$chrA,dat$chrB)
dat$st1 <- as.numeric(dat$st1)
dat$st2 <- as.numeric(dat$st2)
dat$zscore <- as.numeric(dat$zscore)
# remove NA values (missing or non-sig) so the num sig inters calculation does not mess up
dat <- dat %>% drop_na(zscore)
head(dat)
summary(dat)
uchrPair <- unique(dat$chrPair)

#loop through each chrom pair and determine highly interacting regions
fbed_df <- data.frame()
for (p in uchrPair){
#  p = "chr12chr17"
print(p)
  tdat <- dat %>% filter(chrPair == p)
  uchrs <- c(tdat$chrA[1],tdat$chrB[1])
  #######
  #######
  mZdat <- tdat %>% group_by(chrA,st1) %>% dplyr::summarize(mzscore=mean(zscore, na.rm = TRUE))
  tmNdat <- tdat %>% group_by(chrA,st1, cell) %>% dplyr::summarize(nSig=n())
  head(tmNdat)
  mNdat <- tmNdat %>% group_by(chrA,st1) %>% dplyr::summarize(mnSig=mean(nSig, na.rm = TRUE))
  #merge the two dfs
  mdat <- merge(mZdat,mNdat, c("chrA","st1"))
print(mdat)
  #get userc specified nth percentile
  toppercZ <- quantile(mdat$mzscore, probs = as.numeric(perc))
  toppercN <- quantile(mdat$mnSig, probs = as.numeric(perc))
print("top perc z-score")
print(toppercZ)
print("top perc num inters")
print(toppercN)
  #find regions that are above percent cutoff
  highinter <- mdat %>% arrange(st1) %>%
    filter(mzscore >= toppercZ) %>%
    filter(mnSig >= toppercN)
# if there are no interactions that meet the cutoff
if (dim(highinter)[1] == 0) {
 next
}
  print(highinter)
  # go through df and create bed file with the regions
  highinter$tdiff <- c(NA,diff(highinter$st1))
  bedstart <- c()
  bedend <- c()
  tmps <- NA
  tmps <- NA
if (nrow(highinter) <=2){
  for (z in 1:nrow(highinter)){
    tmpe <- as.numeric(highinter$st1[z]) + as.numeric(bin_size)
    tmps <- highinter$st1[z]
    bedstart <- c(bedstart,tmps)
    bedend <- c(bedend,tmpe)
  }#for
 } else {
  for (i in 1:nrow(highinter)) {
    if (is.na(highinter$tdiff[i])){
      tmps <- highinter$st1[i]
    } else {
      # end of the df
      if (i == nrow(highinter)){
        tmpe <- highinter$st1[i] + as.numeric(bin_size)
        bedstart <- c(bedstart,tmps)
        bedend <- c(bedend,tmpe)
      }#if
      # if difference is less than 2 bins, keep going until we get a diff that is >2 bins
      if (highinter$tdiff[i] <= as.numeric(bin_size) *3) {
        next
        # difference is more than 2 bins
      } else {
          tmpe <- highinter$st1[i-1] + as.numeric(bin_size)
          bedstart <- c(bedstart,tmps)
          bedend <- c(bedend,tmpe)
          tmps <- highinter$st1[i]
      }# if else
    }#if else
  }#for
}# if else df length
  chrom <- rep(tdat$chrA[1],length(bedstart))
  chrompair <- rep(tdat$chrPair[1],length(bedstart))
  bed_df <- as.data.frame(cbind(chrom,bedstart,bedend,chrompair))
  fbed_df <- rbind(fbed_df,bed_df)
  #######
  #for chrB
  #######
  mZdat <- tdat %>% group_by(chrB,st2) %>% dplyr::summarize(mzscore=mean(zscore, na.rm = TRUE))
  tmNdat <- tdat %>% group_by(chrB,st2, cell) %>% dplyr::summarize(nSig=n())
  mNdat <- tmNdat %>% group_by(chrB,st2) %>% dplyr::summarize(mnSig=mean(nSig, na.rm = TRUE))
  #merge the two dfs
  mdat <- merge(mZdat,mNdat, c("chrB","st2"))
  #get userc specified nth percentile
  toppercZ <- quantile(mdat$mzscore, probs = as.numeric(perc))
  toppercN <- quantile(mdat$mnSig, probs = as.numeric(perc))
print("top perc z-score")
print(toppercZ)
print("top perc num inters")
print(toppercN)
  #find regions that are above percent cutoff
  highinter <- mdat %>% arrange(st2) %>%
    filter(mzscore >= toppercZ) %>%
    filter(mnSig >= toppercN)
# if there are no interactions that meet the cutoff
if (dim(highinter)[1] == 0) {
 next
}
  print(highinter)
  # go through df and create bed file with the regions
  highinter$tdiff <- c(NA,diff(highinter$st2))
  bedstart <- c()
  bedend <- c()
  tmps <- NA
  tmps <- NA
  if (nrow(highinter) <=2){
    for (z in 1:nrow(highinter)){
      tmpe <- as.numeric(highinter$st2[z]) + as.numeric(bin_size)
      tmps <- highinter$st2[z]
      bedstart <- c(bedstart,tmps)
      bedend <- c(bedend,tmpe)
    }#for
  } else {
    for (i in 1:nrow(highinter)) {
      if (is.na(highinter$tdiff[i])){
        tmps <- highinter$st2[i]
      } else {
        # end of the df
        if (i == nrow(highinter)){
          tmpe <- as.numeric(highinter$st2[i]) + as.numeric(bin_size)
          bedstart <- c(bedstart,tmps)
          bedend <- c(bedend,tmpe)
        }#if
        # if difference is less than 2 bins, keep going until we get a diff that is >2 bins
        if (highinter$tdiff[i] <= as.numeric(bin_size) *2) {
          next
          # difference is more than 2 bins
        } else {
          tmpe <- as.numeric(highinter$st2[i-1]) + as.numeric(bin_size)
          bedstart <- c(bedstart,tmps)
          bedend <- c(bedend,tmpe)
          tmps <- highinter$st2[i]
        }# if else
      }#if else
    }#for
  }# if else df length
  chrom <- rep(tdat$chrB[1],length(bedstart))
  chrompair <- rep(tdat$chrPair[1],length(bedstart))
  bed_df <- as.data.frame(cbind(chrom,bedstart,bedend,chrompair))
  fbed_df <- rbind(fbed_df,bed_df)
} #for each chrom pair

#save bed file to table
write.table(fbed_df, file = as.character(paste0("highly_interacting_regions_chrPairs_",perc,"_percentile.bed")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

##################
# highly interacting regions for one chrom vs all other chroms
##################
#count each interaction as two rows, one for each chrom
datA <- dat %>% select(chrA, st1, cell, zscore, chrPair)
colnames(datA) <- c("chr", "st", "cell", "zscore", "chrPair")
datB <- dat %>% select(chrB, st2, cell, zscore, chrPair)
colnames(datB) <- c("chr", "st", "cell", "zscore", "chrPair")
ddat <- rbind(datA,datB)
uchrom <- unique(ddat$chr)
uchrom
#loop through each chrom and determine highly interacting regions
fbed_df <- data.frame()
for (p in uchrom){
  #  p = "chr12chr17"
  print(p)
  tdat <- dat %>% filter(chr == p)
  #######
  mZdat <- tdat %>% group_by(chr,st) %>% dplyr::summarize(mzscore=mean(zscore, na.rm = TRUE))
  tmNdat <- tdat %>% group_by(chr,st, cell) %>% dplyr::summarize(nSig=n())
  head(tmNdat)
  mNdat <- tmNdat %>% group_by(chr,st) %>% dplyr::summarize(mnSig=mean(nSig, na.rm = TRUE))
  #merge the two dfs
  mdat <- merge(mZdat,mNdat, c("chr","st"))
  print(mdat)
  #get userc specified nth percentile
  toppercZ <- quantile(mdat$mzscore, probs = as.numeric(perc))
  toppercN <- quantile(mdat$mnSig, probs = as.numeric(perc))
  print("top perc z-score")
  print(toppercZ)
  print("top perc num inters")
  print(toppercN)
  #find regions that are above percent cutoff
  highinter <- mdat %>% arrange(st1) %>%
    filter(mzscore >= toppercZ) %>%
    filter(mnSig >= toppercN)
  # if there are no interactions that meet the cutoff
  if (dim(highinter)[1] == 0) {
    next
  }
  print(highinter)
  # go through df and create bed file with the regions
  highinter$tdiff <- c(NA,diff(highinter$st1))
  bedstart <- c()
  bedend <- c()
  tmps <- NA
  tmps <- NA
  if (nrow(highinter) <=2){
    for (z in 1:nrow(highinter)){
      tmpe <- as.numeric(highinter$st1[z]) + as.numeric(bin_size)
      tmps <- highinter$st1[z]
      bedstart <- c(bedstart,tmps)
      bedend <- c(bedend,tmpe)
    }#for
  } else {
    for (i in 1:nrow(highinter)) {
      if (is.na(highinter$tdiff[i])){
        tmps <- highinter$st1[i]
      } else {
        # end of the df
        if (i == nrow(highinter)){
          tmpe <- highinter$st1[i] + as.numeric(bin_size)
          bedstart <- c(bedstart,tmps)
          bedend <- c(bedend,tmpe)
        }#if
        # if difference is less than 2 bins, keep going until we get a diff that is >2 bins
        if (highinter$tdiff[i] <= as.numeric(bin_size) *2) {
          next
          # difference is more than 2 bins
        } else {
          tmpe <- highinter$st1[i-1] + as.numeric(bin_size)
          bedstart <- c(bedstart,tmps)
          bedend <- c(bedend,tmpe)
          tmps <- highinter$st1[i]
        }# if else
      }#if else
    }#for
  }# if else df length
  chrom <- rep(tdat$chrA[1],length(bedstart))
  chrompair <- rep(tdat$chrPair[1],length(bedstart))
  bed_df <- as.data.frame(cbind(chrom,bedstart,bedend,chrompair))
  fbed_df <- rbind(fbed_df,bed_df)
} #for each chrom pair

#save bed file to table
write.table(fbed_df, file = as.character(paste0("highly_interacting_regions_all_chroms_",perc,"_percentile.bed")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

print("DONE")
