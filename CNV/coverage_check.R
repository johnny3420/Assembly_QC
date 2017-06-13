### Looking at read depth across genome and genes

library(dplyr)
library(ggplot2)


## Read and format function
ReadnBin <- function(file){
  tmp <- read.delim(file, header = F) %>%
    group_by(V1) %>%
    count(cut_width(V2, 100000), wt = V3)
  colnames(tmp) <- c("Chromosome","Position","Count")

  tmp$log2 <- log2(tmp$Count)
  tmp$Scaled <- scale(tmp$Count, center = T)
  return(tmp)
}

## CDS regions
genes <- read.delim("Brassica_napus.annotation_v5.bed.gz", header = F) %>%
  filter(V8 == "CDS") %>%
  select(V1:V4) %>%
  rename(Chromosome = V1, Start = V2, End = V3, Gene = V4)

##BWA mapped reads

BWA <- list.files("BWA/",full.names = T)

for(i in BWA){
  temp <- ReadnBin(i)
  all.BWA <- rbind(all.BWA,temp)
}

chrA01 <- read.delim("BWA/BWA.chrA01.coverage.gz", header = F) %>%
  rename(Chromosome = V1, Position = V2, Depth = V3) %>%
  count(cut_width(Position, 100000), wt = Depth)
colnames(chrA01) <- c("Position","Count")

chrA01$log2 <- log2(chrA01$Count)
ggplot(chrA01, aes(Position,log2)) +
  geom_point()

chrA01$Scaled <- scale(chrA01$Count, center = T)

ggplot(chrA01, aes(Position,Scaled)) +
  geom_point()
###
test <- read.delim("BWA.CDS.coverage.gz", header = F) %>%
  rename(Chromosome = V1, Start = V2, End = V3, Depth = V4)
test2 <- cbind(genes, test[,4])
test2 <- test2 %>% rename(Depth = `test[, 4]`)
test2$CDS_Num <- seq_len(nrow(test2))
test2$Depth <- log2(test2$Depth)
plot(test2$CDS_Num,test2$Depth)
