## ---------------------------
##
## Script name: Brachy_NLRome_Initial_Assignment.R
##
## Purpose of script: Produce initial clade assignements based on an NLR tree cut off at 70% NB-ARC coverage.
##
## Author: Daniil Prigozhin
##
## Date Created: 2020-04-19
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


#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
#BiocManager::install("ggtree")
#BiocManager::install("treeio")
#BiocManager::install("tidytree")
#install.packages("tidyverse")

library("tidyverse")
library("ggtree")
library("treeio")
library("tidytree")

'%ni%' <- Negate('%in%')
setwd("~/BioInf/Zenodo/hvNLR/Brachy_NLR_Phylogeny/")

########Import RAXML tree---------------
raxml <- read.raxml("RAxML_Brachy_237min/RAxML_bipartitionsBranchLabels.Raxml.out")
y <- as_tibble(raxml)
which(!is.na(y$label)) 
########Annotate Internal Nodes and save Bootstrap values---------------
z <- mutate(y, label=if_else(is.na(label),paste0("Int",node),label))
bs <- z %>% select(label,bootstrap)

########Export to iTol for rooting---------------
# write.tree(as.phylo(z),"ztree.txt")

########Import from iTol after rooting---------------
tree<-read.newick("YPO1up61DFSutQOmDIgGwg_newick.txt")
w <-as_tibble(tree)
x <- left_join(w,bs)
x
x %>% filter(label == "Int13000")
z %>% filter(label == "Int13000")

###Get number of leaves per node------
N_tips<-vector(length = nrow(x))
for (i in 1:nrow(x)){
  nd <- x$node[[i]]
  N_tips[[i]] <- length(offspring(x,nd, tiponly = T,self_include = T)$node)
}
x <- mutate(x, N_tips = N_tips)
x

## Find clades that are a good size and have best available support------
good_size_clades<-vector()
for (m in tree$tip.label){
  #m <- messy[5,]$node
  #m <- "7416_T436-R1/25-406"
  print(m)
  (ancestry <- ancestor(x,m) %>% filter(N_tips >40, N_tips <500))
  n <- ancestry %>% filter(bootstrap == max(ancestry$bootstrap))
  good_size_clades<-rbind(good_size_clades,n)
}
good_size_clades <- good_size_clades %>% distinct()
good_size_clades
good_size_clades %>% filter(bootstrap <70) %>% arrange(N_tips) %>% print(n=50) ## 10 poor clades

##check that the assignment is unique
a <- good_size_clades$node
anc_pool<-vector()
offs_pool<-vector()
for (nd in a){
  anc <- ancestor(x,nd)$node
  anc_pool <- c(anc_pool,anc)
  offs <- offspring(x,nd)$node
  offs_pool <- c(offs_pool,offs)
}
anc_pool <- unique(anc_pool)
anc_pool
offs_pool <- unique(offs_pool)
offs_pool

partition <- good_size_clades
partition
(partition <- x[a[which(a %ni% offs_pool)],])
partition %>% select(N_tips) %>% sum()
x %>% filter(is.na(bootstrap))
partition %>% arrange(bootstrap) %>% print(n=300)

###Where are the missing tips and why are they missing?
tips<-tree$tip.label
missing<-vector()
for (a in 1:length(tips)){ if (x[a,]$node %ni% offs_pool){missing<-rbind(missing,x[a,])}}
missing 
ancestor(x,4702)

compl<-vector()
for (n in missing$node){
  (ancestry <- ancestor(x,n) %>% filter(N_tips <51))
  c <- ancestry %>% filter(N_tips == max(ancestry$N_tips))
  compl<-rbind(compl,c)
}
compl <- unique(compl)
compl
#compl <- missing
partition <-rbind(partition,compl)

##check that the assignment is unique
a <- partition$node
anc_pool<-vector()
offs_pool<-vector()
for (nd in a){
  anc <- ancestor(x,nd)$node
  anc_pool <- c(anc_pool,anc)
  offs <- offspring(x,nd)$node
  offs_pool <- c(offs_pool,offs)
}
anc_pool <- unique(anc_pool)
offs_pool <- unique(offs_pool)

partition 
x[a[which(a %ni% offs_pool)],]
partition <- x[a[which(a %ni% offs_pool)],]
partition %>% select(N_tips) %>% sum()
length(tips)
partition %>% arrange(N_tips) %>% print(n=300)

ggplot(partition, aes(x=N_tips))+geom_histogram(binwidth = 1)
ggplot(partition, aes(x=N_tips))+geom_histogram(binwidth = 25)
ggplot(partition, aes(x=bootstrap))+geom_histogram(binwidth = 1)

###Export text files for every label in partition40 and populate with properly formated gene id's for automatic retreaval---------
# for (n in 1:(nrow(partition))) {
#   clade <- partition[n,]$label
#   node <- partition[n,]$node
#   tips <- offspring(x,node, tiponly = T, self_include = T)
#   tipnames <- unlist(strsplit(tips$label,"/",fixed = T))[2*1:nrow(tips)-1]
#   write_delim(x = as.data.frame(tipnames), path = paste0("./Autoclades_70/",clade, "_",length(tips$label),".txt"), delim = "\t",quote_escape = "double",append = F,col_names = F)
# }
