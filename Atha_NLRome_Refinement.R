## ---------------------------
##
## Script name: Atha_NLRome_Refinement.R
##
## Purpose of script: Refine clades by splitting into smaller sub-clades
##
## Author: Daniil Prigozhin
##
## Date Created: 2020-03-26
##
## Copyright (c) Daniil Prigozhin, 2020
## Email: daniilprigozhin@lbl.gov
##
## ---------------------------
##
## Notes: This is a script where variable parameters are set on top including input/output directories. 
##        This allows the same script to be used in successive refinements.
##   
##
## ---------------------------
##Installing Packages for alignment manipulation-----------
# install.packages("entropy")
# install.packages("tidyverse")
# if (!requireNamespace("BiocManager"))
#   install.packages("BiocManager")
# BiocManager::install("msa")
# BiocManager::install("odseq")
# BiocManager::install("ggtree")

#Loading libraries-----------------------------------------
library("tidyverse")
library("ggtree")
library("treeio")
library("msa")
library("entropy")
#library(odseq)


##########################################
### Parameters for first refinement-------

setwd("~/BioInf/ArabidopsisRENSEQ/Phylogeny/Autoclades_70/")
MinGapFraction <- 0.9
MinGapBlockWidth <- 1
hvSiteEntCutoff <-  1.5
OutputDirectory <- "../Autoclades_70_Refinement_1/"
##########################################

##########################################
# ### Parameters for second refinement-------
# 
# setwd("~/BioInf/ArabidopsisRENSEQ/Phylogeny/Autoclades_70_Refinement_1/")
# MinGapFraction <- 1
# MinGapBlockWidth <- 1
# hvSiteEntCutoff <-  1.5
# OutputDirectory <- "../Autoclades_70_Refinement_2/"
##########################################

##########################################
# ### Parameters for third refinement-------
# 
# setwd("~/BioInf/ArabidopsisRENSEQ/Phylogeny/Autoclades_70_Refinement_2/")
# MinGapFraction <- 1
# MinGapBlockWidth <- 1
# hvSiteEntCutoff <-  1.5
# OutputDirectory <- "../Autoclades_70_Refinement_3/"
##########################################

##########################################
# ### Parameters for fourth refinement-------
# 
# setwd("~/BioInf/ArabidopsisRENSEQ/Phylogeny/Autoclades_70_Refinement_3/")
# MinGapFraction <- 1
# MinGapBlockWidth <- 1
# hvSiteEntCutoff <-  1.5
# OutputDirectory <- "../Autoclades_70_Refinement_4/"
##########################################

if (!dir.exists(OutputDirectory)){dir.create(OutputDirectory)}

##Collect subdirectories and aligment files ending in *.best.fa
dirs <- dir(path = ".", pattern = "^Int*", recursive = F)
files <- list.files(path = dirs ,pattern = "*.best.fa$", recursive = F, full.names = TRUE) ### This now includes all subsequent subclades

##Calculate entropy with and without gap character, gather alignment stats---------------------
Alph_21 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","-")
Alph_20 <- c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V")

EntropyNG <- vector("list",length(files))
Entropy <- vector("list",length(files))
stats<-vector()
## Entropy calculation
for (i in seq_along(files)) {
  ## Read alignment
  maa <- readAAMultipleAlignment(files[[i]])
  
  ## Extracting folder name
  folder <- str_split(files[[i]], "/")[[1]][1]
  
  ## Filter out gappy columns
  if (MinGapFraction <1){autoMasked <- maskGaps(maa, min.fraction = MinGapFraction, min.block.width = MinGapBlockWidth) ##KEY FILTERING PARAMETERS
  MinAli <- as(autoMasked, "AAStringSet")}else{MinAli<-as(maa, "AAStringSet")}
  ## Calculating Consensus Matrix
  (Tidy_CM<-as_tibble(t(consensusMatrix(MinAli, baseOnly = T))))
  ## Compensating for consensus matrix not keeping full alphabet in output
  for (a in setdiff(Alph_21,colnames(Tidy_CM))){
      vec <- as_tibble(0*(1:nrow(Tidy_CM)))
      colnames(vec) <- paste(a)
      Tidy_CM <- as_tibble(cbind(Tidy_CM,vec))
  } 
  ##Selecting relevant columns
  (Tidy_CM_Gaps <- select(Tidy_CM,all_of(Alph_21)))
  (Tidy_CM_NoGaps <- select(Tidy_CM,all_of(Alph_20)))
    
  ##Entropy Calculation
  ent <- apply(Tidy_CM_Gaps, 1, entropy,unit="log2") %>% as_tibble()
  colnames(ent)<-paste0("Entropy_",folder)
  #ent
  
  ##Entropy Calculation Ignoring Gaps
  entNG <- apply(Tidy_CM_NoGaps, 1, entropy,unit="log2") %>% as_tibble()
  colnames(entNG)<-paste0("EntropyNoGaps_",folder)
  #entNG
  
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
stats %>% arrange(FractionZeroNG) %>% print(n=500)
write_delim(stats, path = paste0(OutputDirectory,"/Stats.txt"), delim = "\t", col_names = T)


##Plot distributions of invariant sites with/without gaps and of highly variable sites---------
ggplot(stats, aes(FractionZeroNG))+geom_histogram()
ggplot(stats, aes(N_HV_Sites))+geom_histogram()#+xlim(-1,50)

ggplot(stats, aes(x=FractionZeroNG, y=N_HV_Sites))+geom_point()+ylim(0,100)
ggplot(stats, aes(x=FractionZeroNG, y=FractionZero))+geom_point()+xlim(0,1)+ylim(0,1)

##Use cutoff value of 90% invariant sites (without gaps) for clades not to be split
DoNotSplit <- stats %>% filter(FractionZeroNG > 0.9) %>% arrange(Clade) %>% print(n=nrow(stats))
stats %>% filter(FractionZeroNG < 0.9) %>% arrange(FractionZeroNG) %>% print(n=nrow(stats))

##Plot Entropy figures for all Clades, save in ./Plots/----------------------
## Set plot hight in bits
m <- max(2.2,na.omit(unlist(Entropy)))
## Make Plots directory if it is not there already
if (!dir.exists("./Plots")){dir.create("./Plots")}
## Entropy plotting
for (l in seq_along(Entropy)){
    folder <- str_remove(colnames(Entropy[[l]]),"Entropy_")
    Ent <- as_tibble(cbind(1:nrow(Entropy[[l]]),Entropy[[l]],EntropyNG[[l]]))
    colnames(Ent)<-c("Position","Entropy", "EntropyNG")

    ggplot(Ent, aes(x = Position))+
      geom_line(aes(y = Entropy), color = "blue")+
      ylim(0,m)+
      ggtitle(paste(folder)) +
      xlab("Position") + 
      ylab("Shannon Entropy with Gaps")+
      theme_classic()
    ggsave(paste0("./Plots","/",folder,"_Entropy_Masked",".pdf"))
    
    ggplot(Ent, aes(x = Position))+
      geom_line(aes(y = EntropyNG), color = "red")+
      ylim(0,m)+
      ggtitle(paste(folder)) +
      xlab("Position") + 
      ylab("Shannon Entropy")+
      theme_classic()
    ggsave(paste0("./Plots","/",folder,"_Entropy_MaskedNG",".pdf"))
}

########################################
###Analyse clade trees------------------
###Import trees, calculate tree statistics, store all trees in BigTable object
'%ni%' <- Negate('%in%')
trees <- list.files(path = dirs,pattern = "RAxML_bipartitionsBranchLabels.*.First.out$", recursive = F, full.names = TRUE)

BigTable<-vector("list",length = length(trees))
for (i in seq_along(trees)){
  #i <- 1
  adress <- trees[[i]]  
  folder <- str_split(adress, "/")[[1]][1]
  tree <- read.raxml(file = paste(adress))
  tree@phylo$tip.label <-  str_remove(tree@phylo$tip.label,"_pacid.*$")
  (table <- as_tibble(tree))
  (table<-mutate(table,Clade = folder))
  ###Get number of unique ecotypes left and right of every node
  ###Get number of duplicated ecotypes left and right of every node 
  x<-table
  alltips<-tree@phylo$tip.label
  CountTable<-vector()
  for (nd in x$node){
  #  nd<-434
    (R_tips <- offspring(x,nd, tiponly=T, self_include=T) %>% select(label) )
    (L_tips <- x[which(alltips %ni% R_tips$label),] %>% select(label))

    if (nrow(L_tips) == 0){L_NE <- 0; L_ND <-0}else{
      L_cc <- strsplit(L_tips$label,'_')
      (L_ecotypes <- unlist(L_cc)[2*(1:length(L_tips$label))-1])
      (L_ecotypes <- na.omit(L_ecotypes))
      L_NE <- length(unique(L_ecotypes))
      L_ND <- length(unique(L_ecotypes[which(duplicated(L_ecotypes))]))
    }

    R_cc <- strsplit(R_tips$label,'_')
    R_ecotypes <- unlist(R_cc)[2*(1:nrow(R_tips))-1]
    (R_ecotypes <- na.omit(R_ecotypes))
    R_NE <- length(unique(R_ecotypes))
    R_ND <- length(unique(R_ecotypes[which(duplicated(R_ecotypes))]))
    
    All_ecotypes<-c(unique(R_ecotypes),unique(L_ecotypes))
    if (nrow(L_tips) == 0){N_OverlapEco <- 0}else{N_OverlapEco <- length(unique(All_ecotypes[which(duplicated(All_ecotypes))]))}
    
    vec <- c(nd,R_NE,R_ND,nrow(R_tips),L_NE,L_ND,nrow(L_tips),N_OverlapEco)
    CountTable <- rbind(CountTable, vec)
  }
  CountTable
  ###Merge with tree data
  colnames(CountTable) <- c("node","R_Eco","R_DuplEco","R_tips","L_Eco","L_DuplEco","L_tips","N_OverlapEco")
  x <- left_join(x,as_tibble(CountTable), by = "node")
  ###Add current tree with node statistics to the BigTable object
  BigTable[[i]] <- x      
}

## Flatten the list of tables into one table 
Full_table<-vector()
for (i in seq_along(BigTable)){Full_table <- rbind(Full_table,BigTable[[i]])} ##This way runs faster than doing rbind within loop above
BigTable <- Full_table

## Check that the table is populated correctly
BigTable %>% filter(!is.na(label)) %>% select(label) %>% distinct() ###Number of genes in the clades under refinement

###Selecting branches in BigTable to cut clades into subclades----------
curr <- "Int8532_350"
tr_tbl <- BigTable %>% filter(Clade == curr)

#Cut_Nodes_10 <- BigTable %>% filter(R_Eco >p, L_Eco>p, bootstrap >94)
a <- BigTable %>% filter(Clade ==curr,node==545) 

###Add branches over 0.5 to the Cut list
LongBranches <- BigTable %>% filter(branch.length > 0.5)
Cut_Nodes_10 <- rbind(Cut_Nodes_10,LongBranches)

Cut_Nodes_10 %>% filter(grepl(Clade,pattern = "Int10720_62$")) %>% arrange(-branch.length)
##To Add to the table
(a <- BigTable %>% filter(grepl(Clade,pattern = "Int9492_85$"), branch.length >0.3)%>% arrange(-branch.length))
Cut_Nodes_10 <- rbind(Cut_Nodes_10,a)

##To Remove from the table
Cut_Nodes_10 <- Cut_Nodes_10 %>% filter(Clade != "Int8392_60")

## Import prior cut list table
Cut_Nodes_10<-read_delim("CutList_Refinement_1.txt",delim = " ",col_names = T)
Cut_Nodes_10<-read_delim("../Autoclades_70/Refined/CutListSplit2_70.txt",delim = " ",col_names = T)
Cut_Nodes_10<-read_delim("../Autoclades_70/Refined/Refined2/CutListSplit3_70.txt",delim = " ",col_names = T)
Cut_Nodes_10 <- Cut_Nodes_10 %>% mutate(label = as.character(label))


Cut_Nodes_10<-left_join(Cut_Nodes_10,Full_table)
Cut_Nodes_10<-left_join(SplitClade,Full_table)
Cut_Nodes_10 %>% arrange(branch.length) %>% print(n=100)
### Write new cut list table. Do Not Use below without redoing subsequent refinements
# write_delim(Cut_Nodes_10, paste0(OutputDirectory,"/CutList_Refinement_1.txt"),quote_escape = "double")
# write_delim(Cut_Nodes_10 ,paste0(OutputDirectory,"/CutList_Refinement_2.txt"),quote_escape = "double")
# write_delim(Cut_Nodes_10 ,paste0(OutputDirectory,"/CutList_Refinement_3.txt"),quote_escape = "double")

## The trees from raxml are unrooted, therefore the program traverses the trees, choosing the ""biggest"" node to split on
## Takes everything behind it as a clade, then takes individual leaves and walks them up to the first Cut_Nodes_10 node
## A new column "Split_1" with the name of that node appended to the original Clade name works as new ubclade name.
## Clade lists can then be exported based on that column. 
## Run over all the nodes in Big Table (the concatenated tree collection), and for each tip assign a new Split_1 clade
## For non-tip lines, give NA. For tips trees with no nodes in the cut list, give _NS_N to signify lack of further splitting.
## For other tips figure out if the current Tip is in the offspring of the Top Clade (one with the most tips)
## If no, assing the name CurClade_TopClade_L
## If yes, find the smallest ancestor in the SplitNodes
## Assign the name CurClade_SmallestAncestor_R

## Append a column to BigTable that contains the TRUE values for nodes in the cut list.
Cut_Nodes_10 <- mutate(Cut_Nodes_10, Split_Node_1 = T)
Cut_Nodes_10
BigTable <- Full_table
BigTable <- left_join(BigTable, Cut_Nodes_10 %>% select(label,node,Clade,Split_Node_1)%>%distinct(), by = c("label","node","Clade"))
BigTable %>% filter(Split_Node_1)
Cut_Nodes_10 %>% filter(Split_Node_1)



## Generate Split_1 column with new clade assignment for every tip of every tree
Split_1 <-vector(length = nrow(BigTable))
for (i in 1:nrow(BigTable)){if (!is.na(BigTable[i,]$label)){                    #works
  Tip <- BigTable[i,]
  CurClade <- BigTable[i,]$Clade
  TreeData <- filter(BigTable,Clade == CurClade)
  SplitNodes <- filter(TreeData,Split_Node_1)
  if (nrow(SplitNodes)>0){
    TopClade <- SplitNodes %>% filter(R_tips == max(SplitNodes$R_tips))
    if (nrow(dplyr::intersect(Tip,offspring(TreeData,TopClade$node,self_include = T))) == 1){
      if (TopClade$R_tips > 1){
        SplitAncestors <- dplyr::intersect(ancestor(TreeData,Tip$node),SplitNodes)
        SmallestAncestor <- SplitAncestors %>% filter(R_tips == min(SplitAncestors$R_tips))
        BestNode <- paste0(CurClade,'_',SmallestAncestor$node,'_R')
      }else{BestNode <- paste0(CurClade,'_',TopClade$node,'_R')}
    }else{BestNode <- paste0(CurClade,'_',TopClade$node,'_L')}
  }else{BestNode <- paste0(CurClade,"_NS_N")}
  Split_1[[i]] <- BestNode
}else{Split_1[[i]] <- NA}}

as_tibble(Split_1) %>% print(n=3000)
(BigTable_1 <- mutate(BigTable, Split_1 = as.factor(Split_1)))

###Comparing results to expectations
BigTable_1 %>% filter(Clade == "Int8532_350", !is.na(Split_1)) %>% group_by(Split_1) %>% summarise(n=n()) %>% print(n=1000)
BigTable_1 %>% filter(Clade == "Int8532_350", !is.na(Split_1)) %>% group_by(Split_1) %>% summarise(n=n()) %>% select(n) %>% sum()
BigTable_1 %>% filter(Clade == "Int12040_327_355_R_173", !is.na(Split_1)) %>% group_by(Split_1) %>% summarise(n=n()) %>% print(n=1000)
BigTable_1 %>% filter(Clade == "Int12040_327_355_R_173", !is.na(Split_1)) %>% group_by(Split_1) %>% summarise(n=n()) %>% select(n) %>% sum()


Cut_Nodes_10 %>% filter(Clade == "Int14919_84") %>% print(n=1000)
stats %>% filter(Clade == "Int13451_126")

## Count numbers of tips in each new clade
CladeStat <- BigTable_1 %>% filter(!is.na(Split_1)) %>% group_by(Clade,Split_1) %>% summarise(n=n()) %>% print(n=300)
CladeStat %>% ungroup() %>% select(Clade) %>% distinct()
CladeStat %>% ungroup() %>% select(Split_1) %>% distinct()
BigTable  %>% select(label) %>% distinct()
write_delim(CladeStat, paste0(OutputDirectory,"/CladeStat.txt"), delim = "\t", append = F, col_names = T)

ggplot(CladeStat,aes(x=n))+geom_density()
ggplot(CladeStat,aes(x=n))+geom_density()+xlim(0,100)
ggplot(CladeStat,aes(x=n))+geom_histogram()

# ## First Refinement. Print out new clade assignment for every tip
# FirstAssignment <- BigTable_1 %>% filter(!is.na(Split_1)) %>% select(label,Split_1)
# write_delim(x = FirstAssignment, path = paste0(OutputDirectory,"/firstassignment.txt") , delim = "\t",quote_escape = "double",append = F,col_names = T)

# ## Second Refinement. Print out new clade assignment for every tip
# SecondAssignment <- BigTable_1 %>% filter(!is.na(Split_1)) %>% select(label,Split_1)
# write_delim(x = SecondAssignment, path =  paste0(OutputDirectory,"/secondassignment.txt") , delim = "\t",quote_escape = "double",append = F,col_names = T)

# ## Third Refinement. Print out new clade assignment for every tip
# ThirdAssignment <- BigTable_1 %>% filter(!is.na(Split_1)) %>% select(label,Split_1)
# write_delim(x = ThirdAssignment, path =  paste0(OutputDirectory,"/thirdassignment.txt") , delim = "\t",quote_escape = "double",append = F,col_names = T)

###Find duplicate tips
Duplicates <- BigTable_1 %>% group_by(label) %>%summarise(n=n()) %>% filter(n>1)
Duplicates <- Duplicates %>% filter(!is.na(label))
Duplicates
BigTable_1 %>% filter(label %in% Duplicates$label) %>% group_by(Split_1) %>% summarise(n=n())
BigTable_1 %>% filter(label %in% Duplicates$label) %>% arrange(label) %>% print(n=50) #%>% group_by(Clade) %>% summarise(n=n())

BigTable %>% filter(!is.na(Clade)) %>% select(Clade) %>% distinct()
BigTable %>% filter(Split_1 =="Int8532_350_354_R", !is.na(label))
Duplicates

###Write Clade Lists - DO NOT RUN without repeating refinement
for (n in 1:(nrow(CladeStat))) {
  clade <- CladeStat[n,]$Split_1
  tips <- BigTable_1 %>% filter(Split_1 == clade)
  tipnames <- tips$label
  write_delim(x = as.data.frame(tipnames), path = paste0(OutputDirectory,"/",clade, "_",length(tips$label),".txt"), delim = "\t",quote_escape = "double",append = F,col_names = F)
}