## ---------------------------
##
## Script name: Atha_NLRome_GeneEntropy.R
##
## Purpose of script: Integrate entropy, hydrophobicity data from alignments with LRR prediction to output "surface maps"
##
## Author: Daniil Prigozhin
##
## Date Created: 2020-04-28
##
## Copyright (c) Daniil Prigozhin, 2020
## Email: daniilprigozhin@lbl.gov
##
## ---------------------------
##
## Notes: Last tested under R version 3.6.3
##   
##
## ---------------------------

  
#Set Working Directory-----------------------------------
setwd(dir = "~/BioInf/ArabidopsisRENSEQ/Phylogeny/Autoclades_70/")
# 
# #Installing Packages for alignment manipulation-----------
# install.packages("entropy")
# install.packages("tidyverse")
# source("http://bioconductor.org/biocLite.R")
# if (!requireNamespace("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install("msa")
# BiocManager::install("odseq")

#Loading libraries-----------------------------------------
library(msa)
library(entropy)
library(odseq)
library(tidyverse)

#Get list of genes to produce Entropy plots for-------
getwd()
(genes<-read_delim("AthaHV.txt",col_names = c("Gene","CommonName"),delim = "\t" ))
genes %>% print(n=200)


#Get the associated clades and alignment files from Atha_NLRome_CladeAnalysis.R output i.e. the main tibble called Common-----
(CladesFiles<- Common %>% filter(Gene %in% genes$Gene))
left_join(genes,CladesFiles) %>% arrange(Clade) %>% select(Gene,CommonName,Clade) %>% print(n=40)
                                                                                            


#Loop that imports alignments one at a time and appends dataframe with columns of entropy values--------
Alph_21 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","-")
Alph_20 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

for (j in seq_along(genes$Gene)) {
  #Reading in sequence alignment-----------------------------
  gene <- genes$Gene[[j]]
  CN <- genes$CommonName[[j]]
  file <- (CladesFiles %>% filter(Gene == gene))$File
  maa <- readAAMultipleAlignment(file)

  #Masking columns by reference gene----------------
  RefGene <- gene
  RefSeq <- as(maa, "AAStringSet")[grep(pattern = gene, x=rownames(maa))]
  AliLen <- width(RefSeq)
  GapMask<- NULL
  for (i in 1:AliLen){
    c<-as.vector(RefSeq[[1]][i]) %in% c("-")
    GapMask<-append(GapMask,c,length(GapMask))
  }
  
  colmask(maa) <- IRanges(GapMask)
  
  #Masking reference genes---------------------------
  RefIDs <- grep(pattern = "Athaliana", x=names(unmasked(maa)))
  rowmask(maa) <- IRanges(start = RefIDs, end = RefIDs)
    
  #Retrieving the non-masked subset------------------
  RefAli <- as(maa, "AAStringSet")
  RefLen <- width(RefAli[1])

  ## Calculating Consensus Matrix
  (Tidy_CM<-as_tibble(t(consensusMatrix(RefAli, baseOnly = T))))
  ## Compensating for consensus matrix not keeping full alphabet in output
  for (a in setdiff(Alph_21,colnames(Tidy_CM))){
    vec <- as_tibble(0*(1:nrow(Tidy_CM)))
    colnames(vec) <- paste(a)
    Tidy_CM <- as_tibble(cbind(Tidy_CM,vec))
  } 
  ##Selecting relevant columns
  (Tidy_CM_NoGaps <- select(Tidy_CM,all_of(Alph_20)))
  
  ##Entropy Calculation Ignoring Gaps
  entNG <- apply(Tidy_CM_NoGaps, 1, entropy,unit="log2") %>% as_tibble()
  colnames(entNG)<-paste0("EntropyNoGaps_",gene)
  entNG
  
  OutputDirectory <- paste0("GeneAnalysisHVPlus/",gene,"_",CN,"/")
  if (!dir.exists(OutputDirectory)){dir.create(OutputDirectory)}
  
  #Output entropy results to file-------------------------------------
  sink(file = paste0(OutputDirectory,gene,"_",CN,".ChimeraEntropy.txt"),append = F)
  cat("attribute: shannonEntropy\n")
  cat("match mode: 1-to-1\n")
  cat("recipient: residues\n")
    for (ii in seq_along(entNG[[1]])){
    cat("\t")
    cat(paste0(":",ii))
    cat("\t")
    cat(sprintf("%.5f", entNG[[1]][ii]))
    cat("\n")
  }
  sink()   
  
  ### Save the alignment filtered by reference ------------------------
  writeXStringSet(as(RefAli, "AAStringSet"), file = paste0(OutputDirectory,CN,"_ReferenceFiltered.pep.fa"))

  ### Generate and print the entropy plot to the same folder-----------  
  Ent <- as_tibble(cbind(1:length(entNG[[1]]),entNG[[1]]))
  colnames(Ent)<-c("Position","Entropy")
  
  ggplot(Ent, aes(x = Position))+
    geom_line(aes(y = Entropy), color = "blue")+
    ylim(0,3)+
    ggtitle(paste0(gene,"_",CN)) +
    xlab("Residue") + 
    ylab("Shannon Entropy")+
    xlim(0,1500)+
    theme_classic()
  ggsave(filename = paste0(OutputDirectory,gene,"_",CN,"_Entropy_Masked",".pdf"),
         width = 15,
         height = 7,
         units = "cm",
         dpi = 300
         )
}

##Import LRR Annotation
LRR_Annotation <- read_delim(file = "GeneAnalysisHV/SnapGene/All_features.tsv", col_names = F,delim = "\t")
# LRR_Annotation <- read_delim(file = "~/Box/NLR_binding_site_prediction/figures/Figure_4/ZAR1.FeatureList.txt", col_names = F,delim = "\t")
colnames(LRR_Annotation) <- c("Gene", "Motiff", "Start", "Stop", "Length")
LRR_Annotation %>% filter(Motiff == "NB-ARC")
LRR_Annotation %>% print(n=600)
getwd()
###For proteins with annotated LRR's, output Entropy and Hydrophobicity scores---------------
RepNum<-NULL

for (i in 1:length(LRR_Annotation$Gene)){
  RepNum <- append(RepNum, base::strsplit(as.character(LRR_Annotation$Motiff[i]),split = "LRR",fixed = T)[[1]][2], length(RepNum))
}
RepNum <-as.integer(RepNum)
LRR_Annotation <- mutate(LRR_Annotation,RepNum = RepNum)
LRR_Annotation %>% mutate(Gene = paste0("Athaliana_",Gene)) -> LRR_Annotation
LRR_Annotation
max(LRR_Annotation$RepNum, na.rm = T) # the greatest number of repeats is 26

#Get the associated clades and alignment files from Atha_NLRome_CladeAnalysis.R output i.e. the main tibble called Common-----
(genes<-read_delim("AthaHV.txt",col_names = c("Gene","CommonName"),delim = "\t" ))
(genes %>% unique() -> genes)
(missing <- filter(genes, Gene %ni% LRR_Annotation$Gene))
(CladesFiles<- Common %>% filter(Gene %in% genes$Gene))
(CladesFiles<- Common %>% filter(Gene %in% LRR_Annotation$Gene))

#Loop that imports alignments one at a time and appends dataframe with columns of entropy values--------
Alph_21 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","-")
Alph_20 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
Hydrophob <- c("A","I","L","M","F","P","W","Y","V")
# genes <- tibble(Gene = "Athaliana_AT3G50950.1", CommonName = "ZAR1")

All_LRR_Entropy <- vector("list", length = nrow(genes))
All_LRR_Hydrophob <-vector("list", length = nrow(genes))
All_LRR_AA <-vector("list", length = nrow(genes))
for (j in seq_along(genes$Gene)) {
  #Reading in sequence alignment-----------------------------
  gene <- genes$Gene[[j]]
  CN <- genes$CommonName[[j]]
  file <- (CladesFiles %>% filter(Gene == gene))$File
  maa <- readAAMultipleAlignment(file)
  
  #Masking columns by reference gene----------------
  RefGene <- gene
  RefSeq <- as(maa, "AAStringSet")[grep(pattern = gene, x=rownames(maa))]
  AliLen <- width(RefSeq)
  GapMask<- NULL
  for (i in 1:AliLen){
    c<-as.vector(RefSeq[[1]][i]) %in% c("-")
    GapMask<-append(GapMask,c,length(GapMask))
  }
  
  colmask(maa) <- IRanges(GapMask)
  (RefSeqNG <- as(maa, "AAStringSet")[grep(pattern = gene, x=rownames(maa))])
  #Masking reference genes---------------------------
  RefIDs <- grep(pattern = "Athaliana", x=names(unmasked(maa)))
  rowmask(maa) <- IRanges(start = RefIDs, end = RefIDs)
  
  #Retrieving the non-masked subset------------------
  RefAli <- as(maa, "AAStringSet")
  RefLen <- width(RefAli[1])
  
  
  ## Calculating Consensus Matrix-----------------
  (Tidy_CM<-as_tibble(t(consensusMatrix(RefAli, baseOnly = T))))
  ## Compensating for consensus matrix not keeping full alphabet in output
  for (a in setdiff(Alph_21,colnames(Tidy_CM))){
    vec <- as_tibble(0*(1:nrow(Tidy_CM)))
    colnames(vec) <- paste(a)
    Tidy_CM <- as_tibble(cbind(Tidy_CM,vec))
  } 
  ##Selecting relevant columns
  (Tidy_CM_NoGaps <- select(Tidy_CM,Alph_20))
  (Tidy_CM_Hydrophob <- select(Tidy_CM_NoGaps, Hydrophob))
  ##Entropy Calculation Ignoring Gaps-----------
  entNG <- apply(Tidy_CM_NoGaps, 1, entropy,unit="log2") %>% as_tibble()
  colnames(entNG)<-paste0("EntropyNoGaps_",gene)
  entNG
  ## Hydrophobicity Claculation--------
  Hphob <- apply(Tidy_CM_Hydrophob, 1, sum) %>% as_tibble()
  Total <- apply(Tidy_CM_NoGaps, 1, sum) %>% as_tibble()
  Hphob/Total %>% as_tibble() -> FracHphob
  colnames(FracHphob)<-paste0("FrHphob_",gene)
  head(FracHphob)
  
  
  ### Get all positions in sink using start coordinates of LRR---------
  (annot <- LRR_Annotation %>% filter(Gene == gene & !is.na(RepNum)))
  LRR_Entropy <-vector()
  LRR_Hphob <- vector()
  LRR_AA <-vector()
  for (ii in 1:length(annot$Motiff))  {
    line <- annot$Motiff[ii]
    (values <- c(-4:13)+annot$Start[ii]-1)
    (LRR_Entropy <- rbind(LRR_Entropy, c(gene,as.character(line),as.numeric(annot$Start[ii]),as.numeric(entNG[[1]][values]))))
    (LRR_Hphob <- rbind(LRR_Hphob, c(gene,as.character(line),as.numeric(annot$Start[ii]),as.numeric(FracHphob[[1]][values]))))  
    (LRR_AA <- rbind(LRR_AA, c(gene,as.character(line),as.numeric(annot$Start[ii]),strsplit(as.character(RefSeqNG[[1]][values]),"")[[1]])))
    }
  All_LRR_Entropy[[j]] <- as_tibble(LRR_Entropy)
  All_LRR_Hydrophob[[j]] <- as_tibble(LRR_Hphob)
  All_LRR_AA[[j]] <- as_tibble(LRR_AA)
  ### Collect statistics on where the high entropy residues are

}
AAlist<-vector()
for(jj in seq_along(All_LRR_Hydrophob)){AAlist<-rbind(AAlist,All_LRR_AA[[jj]])}
colnames(AAlist) <- c("Gene","LRR","Start","Pm5","Pm4","Pm3","Pm2","Pm1","P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12","P13")
AAlist

Hphob<-vector()
for(jj in seq_along(All_LRR_Hydrophob)){Hphob<-rbind(Hphob,All_LRR_Hydrophob[[jj]])}
colnames(Hphob) <- c("Gene","LRR","Start","Pm5","Pm4","Pm3","Pm2","Pm1","P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12","P13")
Hphob

Entropy<-vector()
for(jj in seq_along(All_LRR_Entropy)){Entropy<-rbind(Entropy,All_LRR_Entropy[[jj]])}
colnames(Entropy) <- c("Gene","LRR","Start","Pm5","Pm4","Pm3","Pm2","Pm1","P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","P12","P13")
Entropy

#######################################
## Export resulting tables------------
# write_delim(Entropy,"Entropy_LRR.txt",delim = "\t")
# write_delim(Hphob,"Hphob_LRR.txt",delim = "\t")
# write_delim(AAlist,"ZAR1.AAlist.txt",delim = "\t")