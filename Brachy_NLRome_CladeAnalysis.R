## ---------------------------
##
## Script name: BRachy_NLRome_CladeAnalysis.R
##
## Purpose of script: Synthesize clade information into a single "Common tibble.
##
## Author: Daniil Prigozhin
##
## Date Created: 2020-05-20
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

setwd("~/BioInf/Brachy/Protein/RAxML_237aa/Autoclades_70/")

library(tidyverse)

###IMPORT DATA#################
'%ni%' <- Negate('%in%')


(zeroth <-read_delim("../Autoclades_70/zerothassignment.txt",col_names = c("Clade_0","Gene"), delim = "\t"))
(first <-read_delim("../Autoclades_70_Refinement_1/firstassignment.txt",col_names = c("Clade_1","Gene"), delim = "\t"))
(second <-read_delim("../Autoclades_70_Refinement_2/secondassignment.txt",col_names = c("Clade_2","Gene"), delim = "\t"))
(third <-read_delim("../Autoclades_70_Refinement_3/thirdassignment.txt",col_names = c("Clade_3","Gene"), delim = "\t"))
(fourth <-read_delim("../Autoclades_70_Refinement_4/fourthassignment.txt",col_names = c("Clade_4","Gene"), delim = "\t"))

####Remove duplicate assignment ONCE IT IS KNOWN
#first <- first %>% filter(Clade_1 != "Int9836_17_15_L_16")

###JOIN DATA TOGETHER############
Common<-left_join(zeroth, first, by = "Gene")
Common<-left_join(Common,second,by = "Gene")
Common<-left_join(Common,third,by = "Gene")
Common<-left_join(Common,fourth,by = "Gene")

(Common<-unique(Common[c(2,1,3,4,5,6)])) ###Rearrange Columns

#Common %>% select(Gene) %>% unique() %>% write_delim("NLR_Genes.txt", delim = "", col_names = F)## 11479 unique gene names

###FIND Multiple assignments for individual hits###########################
Common %>% select(Gene) %>% unique()
Common  %>% group_by(Gene) %>% filter(n()>1)  %>% summarise(N=n())%>% arrange(Gene) %>% print(n=300)
Common %>% group_by(Gene) %>% filter(n()>1)  %>% arrange(Gene) %>% print(n=300)
# Just one problem gene ends in two different Clades

######MAKE a Column that uniquely establishes final assigned clade######
CladeA<-vector()
for (l in c(1:nrow(Common))){
  b<-Common[l,]
  a<-paste(b[2], sep ="")
  if (!is.na(b[3])){a<-paste(b[3], sep ="")}
  if (!is.na(b[4])){a<-paste(b[4], sep ="")}
  if (!is.na(b[5])){a<-paste(b[5], sep ="")}
  if (!is.na(b[6])){a<-paste(b[6], sep ="")}
  CladeA<-append(CladeA, a, after = length(CladeA))
}
Common<-mutate(Common, Clade = CladeA)
Common

Common %>% select(Gene, Clade) %>% unique() %>% group_by(Gene) %>% filter(n()>1) %>% print(n=39) #Just one gene with non-unique assignment
Common %>% select(Gene, Clade) %>% unique() 
####Find problem Clades#########
Common %>% select(Gene, Clade) %>% unique() %>% group_by(Gene) %>% filter(n()>1) %>% ungroup() %>% select(Clade) %>% unique()

#####Make a column that contains Ecotype ID only######
Annotation<-read_csv("../../NLR_Map.csv", col_names = c("Name","Assembly"))
(Annotation<-mutate(Annotation, Gene = toupper(Name)))
Common <- left_join(Common, Annotation, by = "Gene")

# ###Need to import flag for highly variable clades####
# Visual <- read_tsv("VisualHV.txt", col_names = T)
# Visual <-Visual[,c(2,3)]
# Visual
# Common <- left_join(Common,Visual, by = "Clade")
# Common


#####Find non-unique gene assignments#########
FixList <- Common %>% select(Gene,Clade) %>% distinct() %>% group_by(Gene) %>% filter(n()>1) %>% select(Gene) %>% distinct()
FixList
Common %>% filter(Gene %in% FixList$Gene) %>% arrange(Gene)
####Need to think about how to fix these. Just 1 at the moment########

###FIND Singletons#########
SingletonClades <- Common %>% group_by(Clade) %>% summarise(n=n()) %>% filter(n<2) %>% select(Clade)
Common %>% filter(Gene %in% (Common %>% filter(Clade %in% SingletonClades$Clade))$Gene) %>% arrange(Gene) %>% group_by(Assembly) %>% summarise(n=n()) %>% print(n=100)
Common %>% filter(Gene %in% (Common %>% filter(Clade %in% SingletonClades$Clade))$Gene) %>% arrange(Gene) %>% group_by(Assembly) %>% summarise(n=n()) %>% select(n) %>% sum()
##Currently 26 assemblies with 37 total singletons. 

###CLADE Sizes####
Frequency_table <- Common %>% group_by(Clade) %>% summarise(n=n()) %>% arrange(n)
Frequency_table %>% print(n=500)
ggplot(Frequency_table, aes (x=n))+geom_density()+theme_classic()
ggplot(Frequency_table, aes (x=n))+geom_density()+theme_classic()+xlim(0,85)
ggplot(Frequency_table, aes (x=n))+geom_bar()+theme_classic()#+xlim(0,85)
Frequency_table
###Currently 409 clades based on the 11,479-tip tree
Common %>% select(Gene) %>% unique()


###Clade Characteristics###
AllClades<-levels(as.factor(Common$Clade))
Common$Clade
AllClades
AllEcotypes<-levels(Common$Assembly)
N_Ecotype<-vector()
NDupl_Ecotype<-vector()
for (l in AllClades){
  #l <- "2_AT4G16920.1"
  (a <- Common %>% filter(Clade == l) %>% select (Gene,Clade,Assembly) %>% distinct() %>% group_by(Assembly) %>% filter(n()>1) %>% distinct(Assembly) %>% nrow())
  (b <- Common %>% filter(Clade == l) %>% select(Assembly) %>% distinct()%>% nrow())
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

# EcoTibble <- left_join(EcoTibble,Visual, by = "Clade")
# EcoTibble %>% arrange(-NDupl_Ecotype)%>% print(n=200)


pd <- position_dodge(width = 1)
EcoTibble
ggplot(EcoTibble %>% filter(N_HV_Sites<100), aes(x=N_Ecotype, y=NDupl_Ecotype,color = N_HV_Sites))+geom_point(position = pd)+ylim(-1,30)

EcoTibble %>% filter(NDupl_Ecotype >0) %>%arrange(N_Ecotype) %>% print(n=200)

# 
# ####Finding hv NLRs
# hvNLR <- Common %>% filter(HV == 1, Ecotype == "Athaliana") %>% print(n=400)
# hvNLR %>% print(n=300)
# write_delim(hvNLR, path = "HV_GENES_FullTable.txt", delim = "\t", na = "NA", append = FALSE, quote_escape = "double")


###NLR-ID#####
ID_table <- read_delim("~/Box/NLR_binding_site_prediction/analyses/Brachy/NLR-ID/Bdistachyon_all.NLR-ID.txt", delim = "\t", col_names = F)
colnames(ID_table) <- c("Name","Domains_ID")
ID_table
All_Dom_table <- read_delim("~/Box/NLR_binding_site_prediction/analyses/Brachy/NLR-ID/Bdistachyon_all.NLR.txt", delim = "\t", col_names = F)
colnames(All_Dom_table) <- c("Name","Domains")
RPW8_NLR <- All_Dom_table %>% filter(grepl("RPW8", Domains)) %>% select(Name) %>% mutate(Status = 1)
write_tsv(RPW8_NLR,"RPW8_NLR.list",col_names = F)

Common_ID <- left_join(Common,ID_table, by = "Name")
Common_ID %>% filter(!is.na(Domains_ID)) %>% group_by(Clade) %>% summarise(n = n()) %>%print(n=300)
Common_AD <- left_join(Common,All_Dom_table, by = "Name")
Common_AD %>% filter(!is.na(Domains)) %>% group_by(Clade) %>% summarise(n = n()) %>%print(n=300)
Common_AD %>% filter(is.na(Domains)) ### There are 27 proteins in my tree that likely got filtered out by KVK

#### Create a Pointer table to link individual genes to assemblies. Add column to Common.
dirs <- dir("..", pattern = "^Autoclades_70.*", full.names = T, recursive = F)
dirs2 <- dir(path=dirs, pattern = "Int.*",full.names = T)
files <- list.files(path = dirs2, pattern = "best.fa$",full.names = T)
files
clades <- unlist(strsplit(files, "/"))[4*1:length(files)-1]
Pointer <- tibble(Clade = clades, File = files)
Pointer

Common <- left_join(Common,Pointer, by = "Clade")
Common %>% filter(is.na(File))
#Currently 37 genes without associated files - this must be the real sinlgeton count - confirmed


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
  (Tidy_CM_Gaps <- select(Tidy_CM,(Alph_21)))
  (Tidy_CM_NoGaps <- select(Tidy_CM,(Alph_20)))
  
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
HV_Clades <- stats %>% filter(N_HV_Sites>=10)
Common %>% filter(Clade %in% HV_Clades$Clade)%>%group_by(Assembly)%>%summarise(n=n())%>% print(n=60)
Common %>% filter(Clade %in% HV_Clades$Clade, Assembly == "v2.1") %>% select(Gene,Clade) %>% unique() %>% print(n=200)
Common %>% filter(Clade %in% HV_Clades$Clade, Assembly == "v2.1") %>% group_by(Clade) %>% summarise(N=n()) %>% print(n=200)
levels(as.factor(Common$Assembly))
BradiHV <- Common %>% filter(Clade %in% HV_Clades$Clade, Assembly == "v2.1") %>% select(Gene) %>% mutate(Tr = str_remove(Gene,"...P$"))%>%select (Tr)%>%unique()%>%print(n=50)

write_delim(BradiHV,"BradiHV.txt",col_names = F)
Common %>% select(Gene) %>% unique()
Common %>% select(Clade) %>% unique()
Common %>% filter(grepl("Bradi1g78926",Name)) %>% select(Clade_2)
# ###Output#############
# a<-EcoTibble %>% select(Clade) %>% distinct()
# a
# write_delim(a, path = "HV_table.txt", delim = "\t", na = "NA", append = FALSE, quote_escape = "double")
# write_delim(Common,path = "Clade_Assignment.txt", delim = "\t", na = "NA", append = FALSE, quote_escape = "double")
# 
# HVClades <- Common %>% filter(Ecotype == "Athaliana" & HV ==1 & Allele ==1) %>% group_by(Clade) %>% summarise(n=n())%>% print (n=40)
# write_delim(HVClades, path = "HV_table.txt", delim = "\t", na = "NA", append = FALSE, quote_escape = "double")
# 
# EcoTibble %>% filter(HV==1) %>% arrange(N_Ecotype)
# Common %>% filter(Gene == "Athaliana_AT1G")
# 
# Common %>% filter(Ecotype == "6939" & Allele ==1) %>% arrange(Gene)%>% select(HV) %>% sum()
# 
# ###Label iTOL tree with the HV column##########
# a <- Common %>% filter(Ecotype =="Athaliana" & Allele == 1) %>% arrange(Gene) %>% select(Clade,Gene, HV) %>% distinct() %>%print(n=300)
# Common                                                                                                  
# EcoTibble
# b<-vector()
# for (a in EcoTibble$Clade){
# #a <- "1_1925_T036-R1"
# b <- append(b,Common %>% filter(Allele==1, Clade==a, Ecotype !="Athaliana") %>% distinct(Gene) %>% nrow(),after = length(b))
# }
# N <-b
# EcoTibble<-as_tibble(cbind(EcoTibble,N))
# EcoTibble %>% filter(HV==1) %>% arrange(N)
# 
# Common %>% filter(Gene =="Athaliana_AT1G31540.1")
# Common
# levels(as.factor(Common$Clade)) %>% length()
# CladeList<-as.factor(Common$Clade)
# for (clade in CladeList){
# #  print(clade)
# #  print(Common %>% filter(Clade == clade) %>%select(Gene))  
#   write_delim(Common %>% filter(Clade == clade) %>%select(Name), path = paste0("~/Box/NLR_binding_site_prediction/analyses/ArabidopsisRENSEQ/Phylogeny/CurrentClades/Full_Gene_Sequences/",clade,".GeneName.txt"), delim = "\t", na = "NA", append = FALSE, quote_escape = "double", col_names =F)
# }
# Common %>% filter(HV==1) %>%select(Clade)%>%distinct()
# write_delim(Common %>% filter(HV==0) %>%select(Clade)%>%distinct(), path = paste0("~/Box/NLR_binding_site_prediction/analyses/ArabidopsisRENSEQ/Phylogeny/CurrentClades/LV_Clade.list"), delim = "\t", na = "NA", append = FALSE, quote_escape = "double", col_names =F)
# write_delim(Common %>% select(Clade,HV) %>%distinct(), path = "~/Box/NLR_binding_site_prediction/analyses/ArabidopsisRENSEQ/Phylogeny/CurrentClades/HV_map.txt", delim = "\t", na = "NA", append = FALSE, quote_escape = "double", col_names =F)
# 
# Number <- Common %>% filter(Ecotype != "Athaliana" & Allele == 1) %>% group_by(Clade) %>% select(Clade,Gene) %>% distinct() %>% summarise(Number=n())
# Number
# write_delim(Number, path = "Number.txt", delim = "\t", na = "NA", append = FALSE, quote_escape = "double", col_names =F)
# 
