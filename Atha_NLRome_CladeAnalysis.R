## ---------------------------
##
## Script name: Atha_NLRome_CladeAnalysis.R
##
## Purpose of script: Synthesize info on clade assignment, ecotypes, alleles, etc. in one "Common" table.
##
## Author: Daniil Prigozhin
##
## Date Created: 2019-08-20
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

setwd("~/BioInf/Zenodo/hvNLR/Athaliana_NLR_Phylogeny/Autoclades_70/")

library(tidyverse)
########################################################################
### Build a Common table of clade assignments for every NLR ------------
########################################################################

'%ni%' <- Negate('%in%')
(zeroth <-read_delim("./partition_0.list",col_names = c("Clade_0","Gene"), delim = "\t"))
(first <-read_delim("../Autoclades_70_Refinement_1/partition_1.list",col_names = c("Clade_1","Gene"), delim = "\t"))
(second <-read_delim("../Autoclades_70_Refinement_2/partition_2.list",col_names = c("Clade_2","Gene"), delim = "\t"))
(third <-read_delim("../Autoclades_70_Refinement_3/partition_3.list",col_names = c("Clade_3","Gene"), delim = "\t"))

####Remove duplicate assignment
(first <- first %>% filter(Clade_1 != "Int9836_17_15_L_16"))

###JOIN DATA TOGETHER############
Common<-left_join(zeroth,first,by = "Gene")
Common<-left_join(Common,second,by = "Gene")
Common<-left_join(Common,third,by = "Gene")
Common
Common<-Common[c(2,1,3,4,5)] ###Rearrange Columns

###FIND Multiple assignments for individual hits###########################
Common %>% select(Gene) %>% unique()
Common  %>% group_by(Gene) %>% filter(n()>1)  %>% summarise(N=n())%>% arrange(Gene) %>% print(n=300)
Common %>% group_by(Gene) %>% filter(n()>1)  %>% arrange(Gene) %>% print(n=300)

######MAKE a Column that uniquely establishes final assigned clade######
CladeA<-vector()
for (l in c(1:nrow(Common))){
  b<-Common[l,]
  a<-paste(b[2], sep ="")
  if (!is.na(b[3])){a<-paste(b[3], sep ="")}
  if (!is.na(b[4])){a<-paste(b[4], sep ="")}
  if (!is.na(b[5])){a<-paste(b[5], sep ="")}
  CladeA<-append(CladeA, a, after = length(CladeA))
}


Common<-mutate(Common, Clade = CladeA)
Common
Common %>% select(Gene, Clade) %>% unique() %>% group_by(Gene) %>% filter(n()>1) %>% print(n=39) 
Common %>% select(Gene, Clade) %>% unique() 
Common %>% select(Clade_0, Clade) %>% filter(grepl("Int9878", Clade_0)) %>% unique() %>% arrange(Clade)

####Find problem Clades#########
Common %>% select(Gene, Clade) %>% unique() %>% group_by(Gene) %>% filter(n()>1) %>% ungroup() %>% select(Clade) %>% unique()

#####Make a column that contains Ecotype ID only######
Ecotype<-vector()
for (l in c(1:nrow(Common))){
  b<-Common[l,]
  a<-strsplit(b$Gene,"_")[[1]][1]
  Ecotype<-append(Ecotype, a, after = length(Ecotype))
}
(Common<-mutate(Common,Ecotype = as.factor(Ecotype)))

###GET allele number##############
Allele <-vector()
for (a in Common$Gene){
  if (grepl("AT.G",a)){b<-strsplit(a,".",fixed=T)[[1]][2]}else if(grepl("-R",a)){b<-strsplit(a,"-R",fixed=T)[[1]][2]}else{b<-vector()}
  Allele <- append(Allele,b,after = length(Allele))
}
(Common <-mutate(Common, Allele = as.integer(Allele)))
Common %>% filter(Allele !=1) %>%print(n=200)
Common %>% filter(Allele ==1, Ecotype != "Athaliana") ###Number of unique genes estimated ~7566

####Get a column of Gene names---------------
LookUp <-vector()
for (Tr in Common$Gene){
  #Tr <- "6909_T215-R1"
  if(grepl(x = Tr,"_T",fixed = T)){
  a <- str_replace(str_remove(Tr, "-R."), "T","G")
  LookUp <- rbind(LookUp, c(Tr,a))  
  } else {LookUp <- rbind(LookUp, c(Tr,Tr)) }
}
colnames(LookUp) <- c("Gene", "Name")
LookUp <- as_tibble(LookUp)
Common <- mutate(Common, Name = LookUp$Name)
Common %>% select(Name) %>% unique() ##7716 unique genes

#####Find non-unique gene assignments#########
FixList <- Common %>% select(Gene,Clade) %>% distinct() %>% group_by(Gene) %>% filter(n()>1) %>% select(Gene) %>% distinct()
FixList
Common %>% filter(Gene %in% FixList$Gene) %>% arrange(Gene)
####Need to think about how to fix these. Just 1 at the moment########

###FIND Singletons#########
SingletonClades <- Common %>% group_by(Clade) %>% summarise(n=n()) %>% filter(n<2) %>% select(Clade)
Common %>% filter(Gene %in% (Common %>% filter(Clade %in% SingletonClades$Clade))$Gene) %>% arrange(Gene)
##Currently just 12 singletons. Can make the fusion gene a singleton and rerun the two clades it appears in to simlplify their alignments

#### Create a Pointer table to link individual genes to assemblies. Add column to Common.
dirs <- dir("..", pattern = "^Autoclades_70.*", full.names = T, recursive = F)
dirs2 <- dir(path=dirs, pattern = "Int.*",full.names = T)
files <- list.files(path = dirs2, pattern = "best.fa$",full.names = T)
files
(clades <- unlist(strsplit(files, "/"))[4*1:length(files)-1])
Pointer <- tibble(Clade = clades, File = files)
Pointer
Common <- left_join(Common,Pointer, by = "Clade")
Common %>% filter(is.na(File))%>% select(Clade) %>% unique() #singleton clades

###############################################
### Look at clade properties ------------------
###############################################

###CLADE Sizes####
Frequency_table <- Common %>% filter(Ecotype != "Athaliana", Allele ==1)%>% group_by(Clade) %>% summarise(n=n()) %>% arrange(n)
Frequency_table <- Common %>% group_by(Clade) %>% summarise(n=n()) %>% arrange(n)
ggplot(Frequency_table, aes (x=n))+geom_density()+theme_classic()
ggplot(Frequency_table, aes (x=n))+geom_density()+theme_classic()+xlim(0,85)
ggplot(Frequency_table, aes (x=n))+geom_bar()+theme_classic()#+xlim(0,85)
Frequency_table %>% select(n) %>% sum() ###sum of 7610 due to filtering for Allele == 1, and excluding reference genes
###Currently 237 clades based on the 7818-tip tree
Frequency_table <- Common %>% group_by(Clade) %>% summarise(n=n()) %>% arrange(n)
Frequency_table %>% filter(n<20) %>% select(n) %>% sum()

###Clade Characteristics###
AllClades<-levels(as.factor(Common$Clade))
Common$Clade
AllClades
AllEcotypes<-levels(Common$Ecotype)
N_Ecotype<-vector()
NDupl_Ecotype<-vector()
for (l in AllClades){
  #l <- "2_AT4G16920.1"
  (a <- Common %>% filter(Ecotype != "Athaliana" & Allele == 1)%>%filter(Clade == l) %>% select (Gene,Clade,Ecotype) %>% distinct() %>% group_by(Ecotype) %>% filter(n()>1) %>% distinct(Ecotype) %>% nrow())
  (b <- Common %>% filter(Ecotype != "Athaliana" & Allele == 1)%>%filter(Clade == l) %>% select(Ecotype) %>% distinct()%>% nrow())
  N_Ecotype<-rbind(N_Ecotype, c(l,b))
  NDupl_Ecotype<-rbind(NDupl_Ecotype, c(l,a))
}
N_Ecotype
colnames(N_Ecotype)<-c("Clade","N_Ecotype")
colnames(NDupl_Ecotype)<-c("Clade","NDupl_Ecotype")
(N_Ecotype<-as_tibble(N_Ecotype))
(NDupl_Ecotype<-as_tibble(NDupl_Ecotype))

EcoTibble<-left_join(N_Ecotype,NDupl_Ecotype,by = "Clade")
EcoTibble$N_Ecotype <- as.numeric(EcoTibble$N_Ecotype)
EcoTibble$NDupl_Ecotype <- as.numeric(EcoTibble$NDupl_Ecotype)
EcoTibble

EcoTibble %>% arrange(-NDupl_Ecotype)%>% print(n=200)

pd <- position_dodge(width = 1)
ggplot(EcoTibble, aes(x=N_Ecotype, y=NDupl_Ecotype))+geom_point(position = pd)

EcoTibble %>% filter(NDupl_Ecotype >0) %>%arrange(N_Ecotype) %>% print(n=200)


##################################################
### Calculate entropy of the clade alignments  ---
##################################################

### Now it is possible to create a list of genes with entropy strings 
### with positions corresponding to the residues of the gene (i.e. remove gaps at the level of the gene)
### Can do this while cycling over Clades to save on reading
library(msa)
library(entropy)
hvSiteEntCutoff <- 1.5
MinGapFraction <- 1 
MinGapBlockWidth <- 1
Alph_21 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","-")
Alph_20 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")
files<-levels(as.factor(Common$File))
EntropyNG <- vector("list",length(files))
Entropy <- vector("list",length(files))
stats<-vector()
## Entropy calculation
for (i in seq_along(files)) {
  ## Read alignment
  
  maa <- readAAMultipleAlignment(files[[i]])
  
  ## Extracting folder name
  str <- str_split(files[[i]], "/")[[1]]
  folder <- str[[length(str)-1]]
  print(paste("Looking at clade ",folder,"..."))
  ## Filter out gappy columns
  
  if ("-" %in% rownames(consensusMatrix(maa))){
    autoMasked <- maskGaps(maa, min.fraction = MinGapFraction, min.block.width = MinGapBlockWidth) ##KEY FILTERING PARAMETERS
    MinAli <- as(autoMasked, "AAStringSet")
  }else{MinAli<-as(maa, "AAStringSet")}
  MinAli
  ## Calculating Consensus Matrix
  (Tidy_CM<-as_tibble(t(consensusMatrix(MinAli, baseOnly = T))))
  ## Compensating for consensus matrix not keeping full alphabet in output
  for (a in setdiff(Alph_21,colnames(Tidy_CM))){
    vec <- as_tibble(0*(1:nrow(Tidy_CM)))
    colnames(vec) <- paste(a)
    Tidy_CM <- as_tibble(cbind(Tidy_CM,vec))
  } 
  ##Selecting relevant columns
  (Tidy_CM_Gaps <- select(Tidy_CM,Alph_21))
  (Tidy_CM_NoGaps <- select(Tidy_CM,Alph_20))
  
  ##Entropy Calculation
  ent <- apply(Tidy_CM_Gaps, 1, entropy,unit="log2") %>% as_tibble()
  colnames(ent)<-paste0("Entropy_",folder)
  
  
  ##Entropy Calculation Ignoring Gaps
  entNG <- apply(Tidy_CM_NoGaps, 1, entropy,unit="log2") %>% as_tibble()
  colnames(entNG)<-paste0("EntropyNoGaps_",folder)
  
  
  ##Save fraction of invariant positions with/without gaps and number of highly variable positions (without gaps)
  frbl <- length(which(ent == 0))/nrow(ent)
  frblng <- length(which(entNG == 0))/nrow(entNG)
  nHVsites <- length(which(entNG > hvSiteEntCutoff))                          ####KEY CUTOFF PARAMETER
  stats <- rbind(stats,c(folder,frbl,frblng,nrow(ent),nHVsites))  
  
  ##Save results of entropy calculation
  Entropy[[i]] <- ent
  EntropyNG[[i]] <- entNG
}

##Format stats tibble-------
colnames(stats)<-c("Clade","FractionZero","FractionZeroNG","Ali_Length","N_HV_Sites")
stats<-as_tibble(stats)
stats
stats<-mutate(stats, Clade = Clade,
              FractionZero = as.numeric(FractionZero),
              FractionZeroNG = as.numeric(FractionZeroNG),
              N_HV_Sites = as.numeric(N_HV_Sites),
              Ali_Length = as.numeric(Ali_Length)
)
stats %>% arrange(-N_HV_Sites) %>% print(n=500)

N_Tips <- Common %>% group_by(Clade)%>%select(Gene)%>%summarise(N=n())
EcoTibble <- left_join(EcoTibble,N_Tips)
EcoTibble <- left_join(EcoTibble,stats)
HV_Clades <- stats %>% filter(N_HV_Sites>9)
HV_Clades
Common <- mutate(Common, HV = ifelse(Common$Clade %in% HV_Clades$Clade,1,0))
Common
#################################
### Export common gene table ----
#################################
write_delim(Common,path = "../Atha_NLRome_GeneTable.txt",delim = "\t")

### Plot highly variable gene per genotype-----
hvFreq <- Common %>% filter(Clade %in% HV_Clades$Clade, Allele == 1)%>%group_by(Ecotype)%>%summarise(n=n())%>% print(n=70)
ggplot(hvFreq,aes(x=n))+geom_histogram()

### 14 of the original 65 clades contain hnNLRs
Common %>% filter(Clade %in% HV_Clades$Clade, Allele == 1) %>% select(Clade_0) %>% unique() %>% arrange()

## Final clade assignment has 237 clades
Common %>% select(Clade) %>% unique()

## Export lists of genes of interest-----
AthaHV<-Common %>% filter(Clade %in% HV_Clades$Clade, Allele == 1,Ecotype == "Athaliana")%>%select(Gene)%>%unique()%>% arrange(Gene)%>%print(n=50)
# write_delim(AthaHV,"AthaHV.txt", col_names = F, delim = "\t", na = "")
AthaHVPlus<-Common %>% filter(Clade %in% HV_Clades$Clade, Allele != 1,Ecotype == "Athaliana")%>%select(Gene)%>%unique()%>% arrange(Gene)%>%print(n=50)
# write_delim(AthaHVPlus,"AthaHVPlus.txt", col_names = F, delim = "\t", na = "")

## Map common names to hvNLRs ----
CommonNames <- read_delim("AthaKnownGenes.txt",delim = "\t", col_names = c("Gene","CommonName"))
CommonNames %>% print(n=100)
AthaHV <- left_join(AthaHV,CommonNames)
AthaHV %>% distinct() %>% arrange(Gene) %>% print(n=50)

## Export a table of clade representation in every ecotype--------
CladeEco_long <- Common %>% select("Name","Clade","Ecotype") %>% unique() %>% group_by(Clade,Ecotype) %>% summarise(N=n()) 
CladeEco_wide <- pivot_wider(CladeEco_long,names_from = Ecotype,values_from = N, values_fill = list (N = 0))  %>% ungroup()
write_delim(CladeEco_wide, "../CladeEco.txt", delim = "\t")



##############################################################
### Produce HV Comparator Tables of clades that have hvNLR ---
##############################################################

Common %>% filter(HV==1)%>%select(Clade_0) %>% distinct()->HV_Clades
HV_Tables <- vector(length = nrow(HV_Clades),mode  = "list")
for (ii in seq_along(HV_Clades$Clade_0)){
  clade <- HV_Clades$Clade_0[[ii]]
  
  HV_Tables[[ii]] <- Common %>% filter(Clade_0 == clade, Allele==1, Ecotype == "Athaliana")%>% select(Gene, Clade_0, Clade, HV, CommonName)
  
}
HV_Table <- HV_Tables[[1]]
for (jj in 2:length(HV_Tables)){
  HV_Table <- rbind(HV_Table,HV_Tables[[jj]])
}
HV_Table
getwd()
write_delim(HV_Table,"HV_Comparators.tsv", delim = t, col_names = FALSE,na = '', quote_escape = "double")
Common %>% left_join(genes) -> Common
Common %>% filter(Gene =="Athaliana_AT3G50950.1") %>% select (Clade_1)



Common %>% filter(HV==1, Ecotype == "Athaliana") %>% arrange(Gene) %>% select(Gene) %>% print(n=40) 