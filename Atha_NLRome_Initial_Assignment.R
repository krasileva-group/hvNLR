## ---------------------------
##
## Script name: Atha_NLRome_Initial_Assignment.R
##
## Purpose of script: Split NB-ARC-based NLRome tree into clades to prepare for de novo alignment
##
## Author: Daniil Prigozhin
##
## Date Created: 2020-03-25
##
## Copyright (c) Daniil Prigozhin, 2020
## Email: daniilprigozhin@lbl.gov
##
## ---------------------------
##
## Notes:
##   import NLRome tree
##   name internal nodes in preparation for rooting in iTol, export
##   import rooted and add bootstrap values back
##   get node statistics - number of tips, ecotypes, and duplicated ecotypes
##   for each tip, select best clade between 40 and 500 tips
##   fill in clades for tips that were missed due to size gate
##   check that each tip is uniquely assigned to a clade
##   output lists of genes for each clade
## ---------------------------

###Setup------
library("tidyverse")
library("treeio")
library("tidytree")
library("ggplot2")

'%ni%' <- Negate('%in%')
setwd("~/BioInf/ArabidopsisRENSEQ/Phylogeny/")

###Import RAXML tree---------------
raxml <- read.raxml("./RAxML_NLRome_237min/RAxML_bipartitionsBranchLabels.Raxml.out")
y <- as_tibble(raxml)
y %>% filter(grepl("66910",label))

###Annotate Internal Nodes and save Bootstrap values---------------
z <- mutate(y, label=if_else(is.na(label),paste0("Int",node),label))
bs <- z %>% select(label,bootstrap)

###Export to iTol for rooting---------------
#write.tree(as.phylo(z),"Autoclades_70/ztree.txt")

###Import from iTol after rooting---------------
tree<-read.newick("7UC1-p5MYQ-WWvafzz-S-Q_newick.txt")
w <-as_tibble(tree)
x <- left_join(w,bs)
###Check that bootstrap values transferred correctly-----
x %>% filter(label == "Int13000")
z %>% filter(label == "Int13000")

###Perform Clade analysis---------------
###Get number of unique ecotypes for every node
###Get number of duplicated ecotypes for every node 

CountTable<-vector()
for (nd in x$node){
  tips <- offspring(x,nd, tiponly=T, self_include=T) %>% select(label) 
  cc <- strsplit(tips$label,'_')
  ecotypes <- unlist(cc)[2*(1:length(tips$label))-1]
  NE <- length(unique(ecotypes))
  ND <- length(unique(ecotypes[which(duplicated(ecotypes))]))
  vec <- c(nd,NE,ND)
  CountTable <- rbind(CountTable, vec)
}
###Merge with tree data
colnames(CountTable) <- c("node","N_Eco","N_DuplEco")
x <- left_join(x,as_tibble(CountTable))

###Get number of leaves per node
count_tips<-vector()
for (nd in x$node){
  n <- c(nd,length(offspring(x,nd, tiponly = T,self_include = T)$node))
  count_tips <-rbind(count_tips,n)
}
colnames(count_tips) <- c("node", "N_tips")
x<-left_join(x, as_tibble(count_tips), by = "node")
x
#################################################################################
####Select Clades between 40 and 500 leaves with the highest bootstrap score-----
good_size_clades<-vector()
for (m in tree$tip.label){
  #print(m)
  (ancestry <- ancestor(x,m) %>% filter(N_tips >40, N_tips <500))
  n <- ancestry %>% filter(bootstrap == max(ancestry$bootstrap))
  good_size_clades<-rbind(good_size_clades,n)
}
good_size_clades <- good_size_clades %>% distinct()
good_size_clades %>% filter(bootstrap <95) %>% arrange(N_tips) %>% print(n=50)

##Collect offspring nodes
offs_pool<-vector()
a <- good_size_clades$node
for (nd in a){offs_pool <- c(offs_pool,offspring(x,nd)$node)}
offs_pool <- unique(offs_pool)

###Remove nesting nodes, check total number of tips
partition <- x[a[which(a %ni% offs_pool)],]
partition %>% select(N_tips) %>% sum()
partition %>% arrange(bootstrap) %>% print(n=300)

###Find missing tips, find biggest clade under 41 tips, add these to the partition
tips<-tree$tip.label
missing<-vector()
for (a in 1:length(tips)){ if (x[a,]$node %ni% offs_pool){missing<-rbind(missing,x[a,])}}
missing ### 187 missing from the first partition

compl<-vector()
for (n in missing$node){
  (ancestry <- ancestor(x,n) %>% filter(N_tips <41))
  c <- ancestry %>% filter(N_tips == max(ancestry$N_tips))
  compl<-rbind(compl,c)
}
compl <- unique(compl)
compl ### 8 additional small clades are enough to cover all sequences

partition <-rbind(partition,compl)

##check that the assignment is unique
a <- partition$node
offs_pool<-vector()
for (nd in a){offs_pool <- c(offs_pool,offspring(x,nd)$node)}
offs_pool <- unique(offs_pool)

partition ###Current partition
x[a[which(a %ni% offs_pool)],] ###Partition checked for nesting clades
partition <- x[a[which(a %ni% offs_pool)],]
partition %>% select(N_tips) %>% sum() ###Check is all tips are there
partition %>% arrange(bootstrap) %>% print(n=300)

###Find missing tips
tips<-tree$tip.label
missing<-vector()
for (a in 1:length(tips)){ if (x[a,]$node %ni% offs_pool){missing<-rbind(missing,x[a,])}}
missing ### shouldn't be any missing tips at this point

###Plot partition----------
ggplot(partition, aes(x=N_tips))+geom_histogram(binwidth = 1)
ggplot(partition, aes(x=N_tips))+geom_histogram(binwidth = 25)
ggplot(partition, aes(x=bootstrap))+geom_histogram(binwidth = 1)

###Export partition--------
write_delim(partition,"Partition.tsv", delim = ";")

###Export node annotation for iTol------
annotation <- partition %>%select(label,bootstrap) %>% mutate(position = 0.5, color = "rgb(0, 0, 0)"	,style = "normal", size = 	2	)
write_tsv(annotation, "Node_annotation.txt",col_names = F)

# ###Export text files for every label in partition40 and populate with properly formated gene id's for automatic retreaval---------
# for (n in 1:(nrow(partition))) {
#   clade <- partition[n,]$label
#   node <- partition[n,]$node
#   tips <- offspring(x,node, tiponly = T, self_include = T)
#   tipnames <- unlist(strsplit(tips$label,"/",fixed = T))[2*1:nrow(tips)-1]
#   write_delim(x = as.data.frame(tipnames), path = paste0("./Autoclades_70/",clade, "_",length(tips$label),".txt"), delim = "\t",quote_escape = "double",append = F,col_names = F)
# }
# 