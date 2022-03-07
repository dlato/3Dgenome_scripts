#CM SNPs interactions
dat <- read.table("CardioMyo_cis50Kb_SNPs_MERGED_overlapping_intearctions_for_circos.txt")
colnames(dat) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
head(dat)

MAP_ULK <- dat %>% filter(chrA == "chr3") %>%
  #MAP4 filter
  filter(((st1 >= 47850000 & st1 <= 48100000) | (st2 >= 47850000 & st2 <= 48100000)) | ((end1 >= 47850000 & end1 <= 48100000) | (end2 >= 47850000 & end2 <= 48100000))) %>%
  #ULK4 filter
  filter(((st1 >= 41700000 & st1 <= 42000000) | (st2 >= 41700000 & st2 <= 42000000)) | ((end1 >= 41700000 & end1 <= 42000000) | (end2 >= 41700000 & end2 <= 42000000)))
MAP_ULK %>% select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")

MAP_SMARCC <- dat %>% filter(chrA == "chr3") %>%
  #MAP4 filter
  filter(((st1 >= 47850000 & st1 <= 48100000) | (st2 >= 47850000 & st2 <= 48100000)) | ((end1 >= 47850000 & end1 <= 48100000) | (end2 >= 47850000 & end2 <= 48100000))) %>%
  #SMARCC1 filter
  filter(((st1 >= 47700000 & st1 <= 47800000) | (st2 >= 47700000 & st2 <= 47800000)) | ((end1 >= 47700000 & end1 <= 47800000) | (end2 >= 47700000 & end2 <= 47800000)))
MAP_SMARCC %>% select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")

ULK_SMARCC <- dat %>% filter(chrA == "chr3") %>%
  #ULK4 filter
  filter(((st1 >= 41700000 & st1 <= 42000000) | (st2 >= 41700000 & st2 <= 42000000)) | ((end1 >= 41700000 & end1 <= 42000000) | (end2 >= 41700000 & end2 <= 42000000)))%>%
  #SMARCC1 filter
  filter(((st1 >= 47700000 & st1 <= 47800000) | (st2 >= 47700000 & st2 <= 47800000)) | ((end1 >= 47700000 & end1 <= 47800000) | (end2 >= 47700000 & end2 <= 47800000)))
ULK_SMARCC %>% select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")

#VSMC SNPs interactions
dat <- read.table("VSMC_cis50Kb_SNPs_MERGED_overlapping_intearctions_for_circos.txt")
colnames(dat) <- c("chrA","st1","end1","chrB","st2","end2","cell","zscore")
head(dat)

MAP_ULK <- dat %>% filter(chrA == "chr3") %>%
  #MAP4 filter
  filter(((st1 >= 47850000 & st1 <= 48100000) | (st2 >= 47850000 & st2 <= 48100000)) | ((end1 >= 47850000 & end1 <= 48100000) | (end2 >= 47850000 & end2 <= 48100000))) %>%
  #ULK4 filter
  filter(((st1 >= 41700000 & st1 <= 42000000) | (st2 >= 41700000 & st2 <= 42000000)) | ((end1 >= 41700000 & end1 <= 42000000) | (end2 >= 41700000 & end2 <= 42000000)))
MAP_ULK %>% select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")

MAP_SMARCC <- dat %>% filter(chrA == "chr3") %>%
  #MAP4 filter
  filter(((st1 >= 47850000 & st1 <= 48100000) | (st2 >= 47850000 & st2 <= 48100000)) | ((end1 >= 47850000 & end1 <= 48100000) | (end2 >= 47850000 & end2 <= 48100000))) %>%
  #SMARCC1 filter
  filter(((st1 >= 47700000 & st1 <= 47800000) | (st2 >= 47700000 & st2 <= 47800000)) | ((end1 >= 47700000 & end1 <= 47800000) | (end2 >= 47700000 & end2 <= 47800000)))
MAP_SMARCC %>% select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")

ULK_SMARCC <- dat %>% filter(chrA == "chr3") %>%
  #ULK4 filter
  filter(((st1 >= 41700000 & st1 <= 42000000) | (st2 >= 41700000 & st2 <= 42000000)) | ((end1 >= 41700000 & end1 <= 42000000) | (end2 >= 41700000 & end2 <= 42000000)))%>%
  #SMARCC1 filter
  filter(((st1 >= 47700000 & st1 <= 47800000) | (st2 >= 47700000 & st2 <= 47800000)) | ((end1 >= 47700000 & end1 <= 47800000) | (end2 >= 47700000 & end2 <= 47800000)))
ULK_SMARCC %>% select(cell) %>% group_by(cell) %>% summarise(n = n(), .groups = "keep")
