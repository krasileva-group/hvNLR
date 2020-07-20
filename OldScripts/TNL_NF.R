#Set Working Directory-----------------------------------
setwd(dir = "~/Box/NLR_binding_site_prediction/analyses/ArabidopsisRENSEQ/Phylogeny/Clades/TNL_NF/")

#Installing Packages for alignment manipulation-----------
install.packages("seqinr")
install.packages("entropy")
install.packages("tidyverse")
source("http://bioconductor.org/biocLite.R")
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("msa")
BiocManager::install("odseq")
install.packages("htmlwidgets")
install.packages("rowr")
install.packages("fractal")

#Loading libraries-----------------------------------------
library(seqinr)
library(msa)
library(entropy)
library(odseq)
library(tidyverse)
library(ggplot2)
library(stringr)
library(htmlwidgets)
library(rowr)
library(fractal)

#The works---------------------------------------------------
#Loop over sub folders in RNL/CNL clade lists, import alignments, maks gaps, output maksed alignment to file,
#Draw the entropy figure with correct title 
#getwd()
list <- list.files(path = ".", pattern = "best.fa", recursive = T)

#Initialize data frame to be appended in the loop-------------
Entropy <- as.data.frame(c(1:5000))
colnames(Entropy) <- c("Position")


#Loop that imports alignments one at a time and appends dataframe with columns of entropy values--------
for (i in list) {
  #i <- c("Athaliana_AT3G44480.1/output.best.fas")
  #Reading in sequence alignment-----------------------------
  
  maa <- readAAMultipleAlignment(paste(i))
  
  #Extracting gene name and folder name--------------------------------------
  name <- str_split(i,"/")[[1]][2]
  folder <- str_split(i, "/")[[1]][1]
  
  #masking gappy sequences and columns--------------------------------------------
#  y <- odseq(maa, threshold = 0.05, distance_metric = "affine", B = 1000)
#  rowmask(maa) <- IRanges(y)
#  ?maskGaps
  maa
  autoMasked <- maskGaps(maa, min.fraction=0.9, min.block.width=3)
  MinAli <- as(autoMasked, "AAStringSet")
  
  #Calculating Consensus Matrix---------------------------------------
  CMmask <-as.data.frame(consensusMatrix(MinAli))
  consensusString(MinAli)
  
  
  #Entropy Calculation------------------------------------------------
  ent <- apply(CMmask, 2, entropy,unit="log2")
  ent <- as.data.frame(ent)
  colnames(ent) <- name

  #head(Entropy)
  Entropy <- rowr::cbind.fill(Entropy, ent, fill = NA)
#  writeXStringSet(as(MinAli, "AAStringSet"), file=paste(folder,"minimal.pep.fa", sep = '/'))
}
#?apply
#At2g01210
Ent <- as_tibble(Entropy)

ggplot(Entropy, mapping = aes(x=Position, 
                              y = AT1G58807.1))+
  geom_line(color = "red")+theme_classic()+
  xlim(1,1200)+
  ylim(0,2.2)


for (l in 1:length(Ent)){
  if (l>1){
    d <- as.data.frame(select(Ent,1,l))
    #head(d)    
    d <- filter(d,!is.na(d[,2]))
    name <- colnames(d)[2]
    folder <- str_split(name,".best")[[1]][1]
 #   folder <- str_split(folder,"X")[[1]][2]
        plot <- 
      ggplot(d, aes(x=d[,1], y = d[,2]))+
      geom_line(color = "blue")+
      ylim(0,2.2)+
      ggtitle(paste(name)) +
      xlab("Position") + 
      ylab("Shannon's Entropy (bits)")+
      theme_classic()
    pdf(paste0(folder, "_Entropy", ".pdf"))
    print(plot)
    dev.off()
  }
}

#Evaluate the plots


F2GeneList <- read_delim("../../../Box/NLR_binding_site_prediction/figures/Figure 2/InterestingGenes.txt","\t")
for (l in F2GeneList$ID){
  l <-"AT3G44480.1"
  line <- filter(F2GeneList, ID == l)
  name <- line$Name
  (d <- as.data.frame(select(Ent,1,l)))
  #head(d)    
  d <- filter(d,!is.na(d[,2]))
  #   name <- colnames(d)[2]
  plot <- 
    ggplot(d, aes(x=d[,1], y = d[,2]))+
    geom_ribbon(ymin = 0, ymax = d[,2], color = "red", fill = "red", alpha = 0.5)+
    ylim(0,3.5)+
    ggtitle(paste(name)) +
    xlab("Position") + 
    ylab("Entropy (bits)")+
    theme_classic()
  print(plot)
  write_delim(d, " ", "\t")
}

