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
  
sink()
sink()
a<-read.jplace("/Users/prigozhin/BioInf/ArabidopsisRENSEQ/Phylogeny/RAxML_NLRome_237min_Outgroup/EPA_z/epa_result.jplace")
big_table <- left_join(as_tibble(a@phylo),a@placements)%>% filter(!is.na(name))
# b<-read_delim("/Users/prigozhin/BioInf/ArabidopsisRENSEQ/Phylogeny/RAxML_NLRome_237min_Outgroup/EPA_z/Int_to_EPA.tsv","\t",col_names = c("label","Number"))
# b<-b[1:7815,]
# left_join(big_table,b)%>%print(n=1000)

### limit placements to >95% likelihood ratio cutoff
#placement <- big_table %>% filter(like_weight_ratio>.95) %>% arrange(name) %>% print(n=200)
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
partition %>% filter(label %in% clade_og$Clade)
partition %>% filter(label %ni% clade_og$Clade)

partition %>% filter(label %ni% clade_og$Clade) %>% select(N_tips) %>% sum() ## 635 tips in initial clades without Capsella or lyrata sequences

clade_og %>% group_by(Clade) %>% summarise(n=n()) %>% arrange(n) %>%print(n=100) -> og_summary
partition %>% filter(label %ni% clade_og$Clade) %>% select(label) %>% distinct()%>% mutate(Clade = label, n=0)%>% select(-label) %>% rbind(og_summary) -> og_summary
og_summary %>% print(n=70)

crub_og <- clade_og %>% filter(grepl("CRUB",name))
crub_og %>% group_by(Clade) %>% summarise(n=n()) %>% arrange(Clade)->crub_summary
partition %>% filter(label %ni% crub_og$Clade) %>% select(label) %>% distinct()%>% mutate(Clade = label, n=0)%>% select(-label)%>%rbind(crub_summary)->crub_summary
crub_summary

alyr_og <- clade_og %>% filter(grepl("ALYR",name))
alyr_og %>% group_by(Clade) %>% summarise(n=n()) %>% arrange(Clade)->alyr_summary
partition %>% filter(label %ni% alyr_og$Clade) %>% select(label) %>% distinct()%>% mutate(Clade = label, n=0)%>% select(-label)%>%rbind(alyr_summary)->alyr_summary

og_summary %>% mutate(Total = n, Crub = crub_summary$n, Alyr = alyr_summary$n) %>% select(-n)->og_summary
og_summary %>% print(n=100)
clade_og %>% mutate(Name = str_remove(name,"\\/.*$"), Clade = as.factor(Clade)) %>% select(-name) -> export_og
export_og
(export_og$Clade %in% str_remove(Common$Clade_0,"_.*"))
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
### Import results -----
jplace_files <- list.files("Autoclades_70/Outgroup", pattern = ".*jplace",recursive = F,full.names = T)
jplace<- read.jplace(jplace_files[[1]])
jplace@placements
as_tibble(jplace)
