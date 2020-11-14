## ---------------------------
##
## Script name: Atha_NLRome_FindNeighbors.R
##
## Purpose of script: find neighbor NLRs using GenomicRanges toolkit
##
## Author: Daniil Prigozhin
##
## Date Created: 2020-11-03
##
## Copyright (c) Daniil Prigozhin, 2020
## Email: daniilprigozhin@lbl.gov
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------
## load packages

require(tidyverse)
BiocManager::install("GenomicRanges")
library(GenomicRanges)
## set working directory

setwd("~/BioInf/ArabidopsisRENSEQ/Proximity/")

data <- read_delim("ImmuneGenes.bed", delim = "\t", col_names = c('chr','start','end','strand','id'))
gr <- GRanges(seqnames = data$chr,ranges = IRanges(data$start,data$end),strand = data$strand, id = data$id)

?GenomicRanges
distances <- vector(length = length(gr))
for (i in 1:length(gr)){
  a <- GenomicRanges::distance(gr[i], gr, select="all",ignore.strand = T)
  b <- tibble(column = a)
  colnames(b) <- gr[i]@elementMetadata@listData$id
  distances[[i]] <- b
  }
distances

dist_matrix <- distances[[1]]
for (jj in 2:length(distances)){dist_matrix<-cbind(dist_matrix,distances[[jj]])}
dist_matrix <- as_tibble(dist_matrix)
names <- colnames(dist_matrix)
dist <- vector()
for (ii in 1:nrow(dist_matrix)){
  for (jj in ii:nrow(dist_matrix)){
    dist <- rbind(dist,c(names[[ii]],names[[jj]],dist_matrix[[ii,jj]]))
  }
}
dist <- as_tibble(dist) 
colnames(dist) <- c("Gene1", "Gene2", "Distance")
dist %>% filter(!is.na(Distance)) %>% mutate(Distance = as.integer(Distance)) %>% filter(Distance >0)->dist
dist %>% filter(Distance < 5000)

## now we have all the distances and need to convert that to a pairs dataset
## important to note that names are of genes not proteins
## this is tricky as we do not want to draw lines where none should exist
## this should motivate removal of alternate spliced proteins from trees
## for now, can avoid the issue by converting all to *.1

dist_1 <- dist %>% mutate(Gene1 = paste0("Athaliana_",Gene1,".1"),Gene2 = paste0("Athaliana_",Gene2,".1"))
dist_1 <- dist_1 %>% filter(Distance < 5000)
dist_1 %>% filter(grepl("12010",Gene1) |grepl("12010",Gene2))

sink("Connections5K.txt",append = F)
cat(
  "DATASET_CONNECTION
SEPARATOR COMMA
DATASET_LABEL,Distance5K
COLOR,#ff0ff0
#if set to 1, arrows will be drawn on the destination end of each connection line
DRAW_ARROWS,0
#when arrows are displayed, this option sets their size
ARROW_SIZE,20
#maximum width specified in the dataset will be drawn with this value. All other widths will be scaled down proportionally.
MAXIMUM_LINE_WIDTH,10
#Angle can be any value between -80 and 80. With angle set to 0 lines will be straight.
CURVE_ANGLE,0
#if CENTER_CURVES is set to 1, center of the tree (in circular display mode only) will be used as the control point for all curves, ignoring 'CURVE_ANGLE'
CENTER_CURVES,1

#if ALIGN_TO_LABELS is set to 1, connections from/to leaf nodes will start/end next to the leaf label, and not on the branch
ALIGN_TO_LABELS,1

DATA
")
for (ii in 1:length(dist_1$Gene1)){
  cat(dist_1[[ii,1]])
  cat(",")
  cat(dist_1[[ii,2]])
  cat(",5,#000000,normal,neighbor\n")
}
sink()
getwd()
