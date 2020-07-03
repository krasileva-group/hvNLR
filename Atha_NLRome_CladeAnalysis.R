## ---------------------------
##
## Script name: Atha_NLRome_CladeAnalysis.R
##
## Purpose of script: Synthesize info on clade assingnement, ecotypes, alleles, etc. in one "Common" table.
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
## Notes:
##   
##
## ---------------------------
## load packages


setwd("~/BioInf/ArabidopsisRENSEQ/Phylogeny/Autoclades_70/")

library(tidyverse)

###IMPORT DATA#################
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

###ADD Column that creates unique ID's for individual HMM hits#############
# NBARC<-paste(Common$Gene,Common$Residues,sep = ":")
# Common<-cbind(Common,NBARC)
# Common <-as_tibble(Common)

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
Common %>% filter(Clade == "1_Int8532_350_545_R_17") %>% select(Gene) %>% unique()


#####Make a column that contains Ecotype ID only######
Ecotype<-vector()
for (l in c(1:nrow(Common))){
  b<-Common[l,]
  a<-strsplit(b$Gene,"_")[[1]][1]
  Ecotype<-append(Ecotype, a, after = length(Ecotype))
}
(Common<-mutate(Common,Ecotype = as.factor(Ecotype)))
Common

###GET allele number##############
Allele <-vector()
for (a in Common$Gene){
  if (grepl("AT.G",a)){b<-strsplit(a,".",fixed=T)[[1]][2]}else if(grepl("-R",a)){b<-strsplit(a,"-R",fixed=T)[[1]][2]}else{b<-vector()}
  Allele <- append(Allele,b,after = length(Allele))
}
(Common <-mutate(Common, Allele = as.integer(Allele)))
Common %>% filter(Allele !=1) %>%print(n=200)
Common %>% filter(Allele ==1, Ecotype != "Athaliana") ###Number of unique genes estimated ~7566

# ###Need to import flag for highly variable clades####
# Visual <- read_tsv("VisualHV.txt", col_names = T)
# Visual <-Visual[,c(2,3)]
# Visual
# Common <- left_join(Common,Visual, by = "Clade")
# Common


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
Common <- left_join(Common,LookUp, by = "Gene")
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
164+73 ## 73 clades with 50 or more genes responsible for 
2998+4886 ##under 50 /over 50
2296+5588 ##under 40 /over 40
795/7884

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

#EcoTibble <- left_join(EcoTibble,Visual, by = "Clade")
EcoTibble %>% arrange(-NDupl_Ecotype)%>% print(n=200)


pd <- position_dodge(width = 1)

ggplot(EcoTibble, aes(x=N_Ecotype, y=NDupl_Ecotype))+geom_point(position = pd)

EcoTibble %>% filter(NDupl_Ecotype >0) %>%arrange(N_Ecotype) %>% print(n=200)

# 
# ####Finding hv NLRs
# hvNLR <- Common %>% filter(HV == 1, Ecotype == "Athaliana") %>% print(n=400)
# hvNLR %>% print(n=300)
# write_delim(hvNLR, path = "HV_GENES_FullTable.txt", delim = "\t", na = "NA", append = FALSE, quote_escape = "double")


###NLR-ID#####
ID_table <- read_delim("../IntegratedDomains/Ath-NLRome.protein.fa_pfamscan-08-26-2019.1e-3.parsed.verbose.NLR-ID.txt", delim = "\t", col_names = F)
colnames(ID_table) <- c("Gene","Domains")
ID_table
All_Dom_table <- read_delim("../IntegratedDomains/Ath-NLRome.protein.fa_pfamscan-08-26-2019.1e-3.parsed.verbose.NLR.txt", delim = "\t", col_names = F)
colnames(All_Dom_table) <- c("Gene","Domains")
All_Dom_table
cur_gene <-"45050"
All_Dom_table %>% filter(grepl(cur_gene,Gene)) %>% select(Domains)
y %>%  filter(grepl(cur_gene,label))
Common_ID <- left_join(Common,ID_table, by = "Gene")
Common_ID %>% filter(!is.na(Domains)) %>% group_by(Clade) %>% summarise(n = n()) %>%print(n=300)



#### Create a Pointer table to link individual genes to assemblies. Add column to Common.
dirs <- dir("..", pattern = "^Autoclades_70.*", full.names = T, recursive = F)
dirs2 <- dir(path=dirs, pattern = "Int.*",full.names = T)
files <- list.files(path = dirs2, pattern = "best.fa$",full.names = T)
files
(clades <- unlist(strsplit(files, "/"))[4*1:length(files)-1])
Pointer <- tibble(Clade = clades, File = files)
Pointer

Common <- left_join(Common,Pointer, by = "Clade")
Common %>% filter(is.na(File))%>% select(Clade) %>% unique()

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
stats

N_Tips <- Common %>% group_by(Clade)%>%select(Gene)%>%summarise(N=n())
EcoTibble <- left_join(EcoTibble,N_Tips)
EcoTibble <- left_join(EcoTibble,stats)
HV_Clades <- stats %>% filter(N_HV_Sites>9)

Common %>% filter(Clade %in% HV_Clades$Clade, Allele == 1)%>%group_by(Ecotype)%>%summarise(n=n())%>% print(n=70)->hvFreq
ggplot(hvFreq,aes(x=n))+geom_histogram()
Common %>% filter(Clade %in% HV_Clades$Clade, Allele == 1) %>% select(Clade_0) %>% unique() %>% arrange()
Common %>% select(Clade) %>% unique()
AthaHV<-Common %>% filter(Clade %in% HV_Clades$Clade, Allele == 1,Ecotype == "Athaliana")%>%select(Gene)%>%unique()%>% arrange(Gene)%>%print(n=50)
write_delim(AthaHV,"AthaHV.txt", col_names = F, delim = "\t", na = "")
AthaHVPlus<-Common %>% filter(Clade %in% HV_Clades$Clade, Allele != 1,Ecotype == "Athaliana")%>%select(Gene)%>%unique()%>% arrange(Gene)%>%print(n=50)
write_delim(AthaHVPlus,"AthaHVPlus.txt", col_names = F, delim = "\t", na = "")

CommonNames <- read_delim("AthaKnownGenes.txt",delim = "\t", col_names = c("Gene","CommonName"))
CommonNames %>% print(n=100)
AthaHV <- left_join(AthaHV,CommonNames)


AthaHV %>% arrange(Gene) %>% print(n=50)
Common %>% select("Clade")%>%unique()

Common %>% filter(grepl("AT1G63360.1",Gene))%>%select(Clade)
HV_Clades

mutate(Common, HV = ifelse(Common$Clade %in% HV_Clades$Clade,1,0)) ->Common
Common %>% select(Clade,HV) %>% unique() %>% arrange(Clade) %>% write_delim("AllCladesHV.txt",delim = "\t")



## Would like to look for recent duplications--------
Common %>% select("Name","Clade","Ecotype") %>% unique() %>% group_by(Clade,Ecotype) %>% summarise(N=n()) ->CladeEco_long

pivot_wider(CladeEco_long,names_from = Ecotype,values_from = N, values_fill = list (N = 0))  %>% ungroup()-> CladeEco_wide


ggplot(CladeEco_Ent,aes(x = Ent))+geom_histogram()
filter(CladeEco_Ent, Ent ==0)
write_delim(CladeEco_wide, "CladeEco.txt", delim = "\t")
Common %>% filter(Clade =="Int14314_286_309_R_81_56_L_80_97_R_56", Ecotype=="Athaliana") %>% filter(Allele ==1) %>% select(Gene)
Common
AthaHV
#







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
