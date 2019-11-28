#Set Working Directory-----------------------------------
setwd(dir = "~/Box/NLR_binding_site_prediction/analyses/ArabidopsisRENSEQ/Phylogeny/HV Clades/")
# 
# #Installing Packages for alignment manipulation-----------
# install.packages("seqinr")
# install.packages("rowr")
# install.packages("tidyverse")
# source("http://bioconductor.org/biocLite.R")
# if (!requireNamespace("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install("msa")
# BiocManager::install("odseq")
# install.packages("htmlwidgets")
# install.packages("rowr")
# install.packages("fractal")

#Loading libraries-----------------------------------------
library(reshape2)
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

#Loop that imports alignments one at a time and appends dataframe with columns of entropy values--------
for (j in list) {
  #Reading in sequence alignment-----------------------------
  #j<-"AT1G31540.2/AT1G31540.1.best.fa"
  maa <- readAAMultipleAlignment(paste(j))
    
  #Extracting alignment name and folder name--------------------------------------
  name <- str_split(j,"/")[[1]][2]
  folder <- str_split(j, "/")[[1]][1]
  folder
  #Masking columns by reference gene----------------
  RefGene <- folder
  #grep(pattern = folder, x=rownames(maa))
  #maa
  #nrow(maa)
  #ncol(maa)
  RefSeq <- as(maa, "AAStringSet")[grep(pattern = folder, x=rownames(maa))]
  AliLen <- width(RefSeq)
  GapMask<- NULL
  for (i in 1:AliLen){
    c<-as.vector(RefSeq[[1]][i]) %in% c("-")
    GapMask<-append(GapMask,c,length(GapMask))
  }
  
  colmask(maa) <- IRanges(GapMask)
  
  #Masking reference genes---------------------------
  RefIDs <- grep(pattern = "pacid", x=names(unmasked(maa)))
  rowmask(maa) <- IRanges(start = RefIDs, end = RefIDs)
    
  #Retrieving the non-masked subset------------------
  RefAli <- as(maa, "AAStringSet")
  RefLen <- width(RefAli[1])

  #Calculating Consensus Matrix---------------------------------------
  CMmask <-as.data.frame(consensusMatrix(RefAli))
  NoGaps<-CMmask[c(2:20,22),]
  
  #Entropy Calculation------------------------------------------------
  #ent <- apply(CMmask, 2, entropy,unit="log2")
  #ent <- as.data.frame(ent)
  #colnames(ent) <- folder
  #ent
  #Entropy Calculation Ignoring Gaps----------------------------------
  entNG <- apply(NoGaps, 2, entropy,unit="log2")
  entNG <- as.data.frame(entNG)
  #colnames(entNG) <- folder
  
  #frame <- cbind(1:RefLen,entNG)
  #ggplot(frame,aes(x=frame[,1], y = frame[,2]))+
  #  geom_line(color = "blue")+theme_classic()
    
  #Output entropy results to file-------------------------------------
  sink(file = paste0(folder,"/",folder,".ChimeraEntropy.txt"),append = F)
  cat("attribute: shannonEntropy\n")
  cat("match mode: 1-to-1\n")
  cat("recipient: residues\n")
    for (j in 1:(RefLen-1)){
    cat("\t")
    cat(paste0(":",j))
    cat("\t")
    cat(sprintf("%.5f", entNG[[1]][j]))
    cat("\n")
  }
  sink()   
  #file.show(paste0(folder,"/",folder,".Entropy.txt"))
  
  ###Save the alignment filtered by reference ------------------------
  writeXStringSet(as(RefAli, "AAStringSet"), file = paste(folder,"ReferenceFiltered.pep.fa", sep = '/'))
}

##Import LRR Annotation
LRR_Annotation <- read.csv(file = "LRR_Annotation.csv", header = F, sep = ";")
colnames(LRR_Annotation) <- c("Gene", "Motiff", "Start", "Stop", "Length")
LRR_Annotation <- as_tibble(LRR_Annotation)
LRR_Annotation %>% filter(Motiff == "NB-ARC")


###For proteins with annotated LRR's, output Entropy and Hydrophobicity scores---------------
RepNum<-NULL

for (i in 1:length(LRR_Annotation$Gene)){
  RepNum <- append(RepNum, base::strsplit(as.character(LRR_Annotation$Motiff[i]),split = "LRR",fixed = T)[[1]][2], length(RepNum))
}
RepNum <-as.integer(RepNum)
LRR_Annotation <- cbind(LRR_Annotation,RepNum) %>% as_tibble()
LRR_Annotation
max(LRR_Annotation$RepNum, na.rm = T)

Entropy<-vector()
for (i in levels(LRR_Annotation$Gene)){
#  i<-"AT1G31540.1"
  alignment <- readAAMultipleAlignment(filepath = paste(i,"ReferenceFiltered.pep.fa", sep = '/'))
  CMmask <-as.data.frame(consensusMatrix(alignment))
  NoGaps<-CMmask[c(2:20,22),]
  entNG <- apply(NoGaps, 2, entropy,unit="log2")
  entNG
  annot <- LRR_Annotation %>% filter(Gene == i & !is.na(RepNum))  
  for (j in 1:length(annot$Motiff))  {
    #j<-3
    line <- annot$Motiff[j]
    values <- c(2,3,5,7,8)+annot$Start[j]-1
    Entropy <- rbind(Entropy, c(i,as.character(line),as.numeric(entNG[values])))
  }
}

Entropy <-as_tibble(Entropy)
names(Entropy)<-c("Gene","Repeat","P2","P3","P5","P7","P8")
Entropy

Entropy <- mutate(Entropy, P2 = as.numeric(P2))
Entropy <- mutate(Entropy, P3 = as.numeric(P3))
Entropy <- mutate(Entropy, P5 = as.numeric(P5))
Entropy <- mutate(Entropy, P7 = as.numeric(P7))
Entropy <- mutate(Entropy, P8 = as.numeric(P8))
Entropy <- mutate(Entropy, Gene = as.factor(Gene))
Entropy <- mutate(Entropy, Repeat = as.factor(Repeat))
Entropy
levels(Entropy$Repeat)
Entropy %>% print(n=400)


LRR_Annotation %>% filter(Gene == "AT1G31540.2")

###Cleanup required for AT1G31540.1 and AT1G61310.1-------------
#Will ignore for the time being



########Plot heatmap with geom tile
for (Name in levels(LRR_Annotation$Gene)){
#Name <- "AT1G58807.1"
plot <- Entropy %>% filter(Gene == Name) %>% select(Repeat,P2,P3,P5,P7,P8)
plot.m <- melt(plot)
plot.m <- plot.m %>% mutate(N = as.integer(str_remove(as.character(Repeat), "LRR"))) %>% arrange(N)
p <- ggplot(plot.m, aes(variable, N)) + 
  geom_tile(aes(fill = value), colour = "white") + 
  scale_fill_gradient2(low = "#1a9641", mid = "#ffffbf",high = "#d7191c", midpoint = 1)
print(p+theme_void()+ggtitle(paste(Name))+theme(legend.position = "none"))
}
write_csv(Entropy,"Entropy.txt")


HPhob<-vector()
for (i in levels(LRR_Annotation$Gene)){
  #  i<-"AT1G31540.1"
  alignment <- readAAMultipleAlignment(filepath = paste(i,"ReferenceFiltered.pep.fa", sep = '/'))
  as_tibble(consensusMatrix(alignment))
  CMmask <-as.data.frame(consensusMatrix(alignment))
  NoGaps<-CMmask[c(2:20,22),]
  Hydrophobes<-CMmask[c(6,9,11,12,14,19,20,22),]
  (SumsNG <-apply(NoGaps, 2, sum))
  (SumsHP <- apply(Hydrophobes,2,sum) ) 
  FractionHP <- SumsHP/SumsNG
    annot <- LRR_Annotation %>% filter(Gene == i & !is.na(RepNum))  
  for (j in 1:length(annot$Motiff))  {
    #j<-3
    line <- annot$Motiff[j]
    values <- c(2,3,5,7,8)+annot$Start[j]-1
    HPhob <- rbind(HPhob, c(i,as.character(line),as.numeric(FractionHP[values])))
  }
}

HPhob <-as_tibble(HPhob)
names(HPhob)<-c("Gene","Repeat","P2","P3","P5","P7","P8")
HPhob

HPhob <- mutate(HPhob, P2 = as.numeric(P2))
HPhob <- mutate(HPhob, P3 = as.numeric(P3))
HPhob <- mutate(HPhob, P5 = as.numeric(P5))
HPhob <- mutate(HPhob, P7 = as.numeric(P7))
HPhob <- mutate(HPhob, P8 = as.numeric(P8))
HPhob <- mutate(HPhob, Gene = as.factor(Gene))
HPhob <- mutate(HPhob, Repeat = as.factor(Repeat))
HPhob
levels(HPhob$Repeat)
HPhob %>% print(n=400)


########Plot heatmap with geom tile
for (Name in levels(LRR_Annotation$Gene)){
  #Name <- "AT1G58807.1"
  plot <- HPhob %>% filter(Gene == Name) %>% select(Repeat,P2,P3,P5,P7,P8)
  plot.m <- melt(plot)
  plot.m <- plot.m %>% mutate(N = as.integer(str_remove(as.character(Repeat), "LRR"))) %>% arrange(N)
  p <- ggplot(plot.m, aes(variable, N)) + 
    geom_tile(aes(fill = value), colour = "white") + 
    scale_fill_gradient(low = "white",high = "red")
  print(p+theme_void()+ggtitle(paste(Name))+theme(legend.position = "none"))
}
write_csv(HPhob,"HPhob.txt")



###Make a simple product combined measure----------
HPhob.m<-melt(HPhob)
Entropy.m<-melt(Entropy)
Combi <- mutate(HPhob.m,value = value*Entropy.m$value)
Combi

for (Name in levels(LRR_Annotation$Gene)){
  #Name <- "AT1G58807.1"
  plot.m <- Combi %>% filter(Gene == Name) 
  plot.m <- plot.m %>% mutate(N = as.integer(str_remove(as.character(Repeat), "LRR"))) %>% arrange(N)
  p <- ggplot(plot.m, aes(variable, N)) + 
    geom_tile(aes(fill = value), colour = "white") + 
    scale_fill_gradient(low = "white",high = "blue")
  print(p+theme_void()+ggtitle(paste(Name))+theme(legend.position = "none"))
}
