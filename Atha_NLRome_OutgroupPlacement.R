## ---------------------------
##
## Script name: Atha_NLRome_OutgroupPlacement.R
##
## Purpose of script: Place outgroup sequences in the NLRome tree, output tree annotations and propagate into clade trees 
##
## Author: Daniil Prigozhin
##
## Date Created: 2020-11-02
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
require(tidytree)
require(treeio)
## set working directory

setwd("~/BioInf/ArabidopsisRENSEQ/Phylogeny/")  
  
a<-read.jplace("/Users/prigozhin/BioInf/ArabidopsisRENSEQ/Phylogeny/RAxML_NLRome_237min_Outgroup/EPA_z/epa_result.jplace")
big_table <- left_join(as_tibble(a@phylo),a@placements)%>% filter(!is.na(name))

placement <- big_table

### import NLRome tree rooted in iTOL
placement %>% print(n=500)
tree_pl <- left_join(x,placement)

placement %>% filter(grepl("CRUB",name))->crub_place
placement %>% filter(grepl("ALYR",name))->alyr_place
alyr_place %>% filter(label == "Int11609")
### Print symbol dataset ----
sink("alyr_placement.txt", append = F)
cat("DATASET_SYMBOL
SEPARATOR SPACE

#label is used in the legend table (can be changed later)
DATASET_LABEL A.lyrataPlacement

#dataset color (can be changed later)
COLOR #FFFF00
MAXIMUM_SIZE 10
DATA
")
ii<-1
for (ii in 1:length(alyr_place$label)){
  cat(alyr_place[ii,4] %>% unlist())
  cat(" 2 10 #ffff00 1 0.25 ")
  cat(as.character(alyr_place$name)[ii] %>% unlist())
  cat("\n")
}
sink()
#
### Print symbol dataset ----
sink("crub_placement.txt", append = F)
cat("DATASET_SYMBOL
SEPARATOR SPACE

#label is used in the legend table (can be changed later)
DATASET_LABEL C.rubellaPlacement

#dataset color (can be changed later)
COLOR #00FF00
MAXIMUM_SIZE 10
DATA
")
ii<-1
for (ii in 1:length(crub_place$label)){
  cat(crub_place[ii,4] %>% unlist())
  cat(" 2 10 #00ff00 1 0.75 ")
  cat(as.character(crub_place$name)[ii] %>% unlist())
  cat("\n")
}
sink()
#

## for every name in placement find an ancestor in partition
placement %>% print(n=400)
cname <- "CRUB_R0I9X3/50-376"
find_clade_for_og <- function(cname){
  lines <- placement %>% filter(name == cname) %>% arrange(-like_weight_ratio)
  lab <- lines[1,4]
  node <- x %>% filter(label %in% lab) %>% select(node) %>% unlist
  ancestors <- ancestor(x,node) %>% union(x[node,])
  label <- intersect(ancestors,partition) %>% select(label)
  if(nrow(label) == 0){return(NA)} else if(nrow(label ==1)){return(unlist(label))}else {return("MULTI")}
}

test <- placement %>% mutate(Clade = lapply(placement$name, find_clade_for_og)%>%unlist)
test %>% print(n=400)

test %>% filter(grepl("CRUB_R0GTC0",name)) %>% select(Clade) %>% unlist
Summary <- test %>% select(Clade,name) %>% distinct %>% group_by(Clade)%>%summarise(Number = n())
Alyr <- test %>% filter(grepl("ALYR",name)) %>% select(Clade,name) %>% distinct %>% group_by(Clade)%>%summarise(Number_alyr = n())
Crub <- test %>% filter(grepl("CRUB",name)) %>% select(Clade,name) %>% distinct %>% group_by(Clade)%>%summarise(Number_crub = n())

All<-left_join(left_join(Summary,Alyr),Crub) %>% print(n=100)
All[54,1]<-"None"
All[is.na(All)]<-0
All %>% print(n=100)
All_clades <- filter(partition, label %ni% All$Clade)%>%select(label) %>% mutate(Clade = label) %>% select (Clade) %>% mutate(Number=0, Number_alyr=0, Number_crub=0) %>% rbind(All)
All_clades %>% arrange(Clade)->All_clades
All_clades %>% print(n=100)
test %>% select(Clade,name) %>% distinct %>% filter(Clade == "Int14015") %>% print(n=300)
Common %>% select(Clade_0) %>% arrange %>% distinct %>% write_delim("~/Desktop/Clade_names.txt", delim = "\t")
All_clades %>% arrange(Clade) %>% write_delim("~/Desktop/Clade_numbers.txt", delim = "\t")
Common %>% filter(Allele ==1, Ecotype== "Athaliana") %>% group_by(Clade_0) %>% summarise(n=n()) %>% write_delim("~/Desktop/Col0_numbers.txt", delim = "\t")
Named_genes<-read_delim("/Users/prigozhin/Box/NLR_binding_site_prediction/analyses/ArabidopsisRENSEQ/Phylogeny/Autoclades_70/AthaKnownGenes.txt", delim = "\t", col_names = c("Gene","CommonName"))
left_join(Common,Named_genes)%>%filter(!is.na(CommonName))%>%select(Clade_0,CommonName)%>% arrange(Clade_0)%>%write_delim("~/Desktop/Repres.txt", delim = "\t")

Common %>% select(Clade_0,Clade) %>% distinct %>% filter(Clade_0 =="Int10969_308")

### for every node, make a list of outgroup sequences placed at or below the node
(placement$label %in% x$label)
x_pl<- left_join(x,placement %>% mutate(OG = name) %>% select(OG,label))
x_pl %>% filter(!is.na(OG)) %>% print(n=200)

placement %>% filter(label == "Int9553")
x %>% filter(label == "Int9553")

clade_og <- vector()
for (clade in partition$label){
  nd <- x %>% filter(label == clade)
  node_set <- offspring(x,nd$node,tiponly = F, self_include = T) %>% select(node,label)
  line <- placement %>% filter(label %in% node_set$label) %>% select(name) %>% distinct() %>% mutate(Clade=clade)
  clade_og <- rbind(line,clade_og)
}
clade_og
(clade_og$Clade %in% partition$label)
(clade_og$Clade %in% Summary$Clade)
partition %>% filter(label %in% clade_og$Clade)
partition %>% filter(label %ni% clade_og$Clade)

partition %>% filter(label %ni% clade_og$Clade) %>% select(N_tips) %>% sum() ## 635 tips in initial clades without Capsella or lyrata sequences
partition %>% filter(label %ni% Summary$Clade) %>% select(N_tips) %>% sum() ## 635 tips in initial clades without Capsella or lyrata sequences

clade_og %>% mutate(Name = str_remove(name,"\\/.*$"), Clade = as.factor(Clade)) %>% select(-name) -> export_og
export_og
(export_og$Clade %in% str_remove(Common$Clade_0,"_.*"))

## Use Common table to output outgroup ID's for every node

Common %>% mutate(Clade = str_remove(Common$Clade_0,"_.*")) %>% 
  select(Clade_0, Clade) %>% 
  distinct() %>% 
  left_join(export_og) %>% select(-Clade) -> export_og
export_og %>% mutate(Clade_0 = as.factor(Clade_0)) %>% filter(!is.na(Name))->export_og

for (clade in levels(export_og$Clade_0)){
    tips <- export_og %>% filter(Clade_0 ==clade) %>% select(Name)
    if (nrow(tips) >0){print(clade)
    print(tips)
    write_delim(x = as.data.frame(tips), 
                path = paste0("./Autoclades_70/",clade,".outgroup.txt"), 
                delim = "\t",quote_escape = "double",
                append = F,col_names = F)
    }
}
### Import jplace files for refined clades -----
jplace_files <- list.files("Autoclades_70/Outgroup", pattern = ".*jplace",recursive = F,full.names = T)
jplace<- read.jplace(jplace_files[[1]])
jplace@placements
as_tibble(jplace) %>% print(n=2000)
as_tibble(jplace) 
