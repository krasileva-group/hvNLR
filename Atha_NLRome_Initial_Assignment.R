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
##   Last Tested on R version 3.6.3
## ---------------------------
# sessionInfo()

## Install packages. This needs to be done once for each R installation.
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("treeio")
# install.packages("tidyverse")
# install.packages("tidytree")

###Setup------

library("tidyverse")
library("treeio")
library("tidytree")

'%ni%' <- Negate('%in%')
setwd("~/BioInf/Zenodo/hvNLR/Athaliana_NLR_Phylogeny/")

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
x

###Check that bootstrap values transferred correctly-----
x %>% filter(label == "Int13000")
z %>% filter(label == "Int13000")

###Perform Clade analysis---------------
###Get number of unique ecotypes for every node
###Get number of duplicated ecotypes for every node 

CountTable<-vector(length = length(x$node),mode = "list")
for (ii in seq_along(CountTable)){
  tips <- offspring(x,ii, tiponly=T, self_include=T) %>% select(label) 
  cc <- strsplit(tips$label,'_')
  ecotypes <- unlist(cc)[2*(1:length(tips$label))-1]
  NE <- length(unique(ecotypes))
  ND <- length(unique(ecotypes[which(duplicated(ecotypes))]))
  CountTable[[ii]] <- c(ii,NE,ND)
}
###Merge with tree data
CT <- vector()
for (jj in seq_along(CountTable)){CT <- rbind(CT,CountTable[[jj]])}
(x <- mutate(x,N_Eco = CT[,2],N_DuplEco=CT[,3]))

###Get number of leaves per node
count_tips<-vector(length = length(x$node))
for (ii in seq_along(x$node)){count_tips[[ii]]  <- length(offspring(x,ii, tiponly = T,self_include = T)$node)}
(x<-mutate(x, N_tips = count_tips))

#################################################################################
####Select Clades between 40 and 500 leaves with the highest bootstrap score-----
good_size_clades<-vector(mode = "list",length(tree$tip.label))
for (ii in seq_along(tree$tip.label)){
  m<-tree$tip.label[[ii]]
  b<- length(tree$tip.label)
  print(paste0("Node ",ii, " of ",b))
  (ancestry <- ancestor(x,m) %>% filter(N_tips >40, N_tips <500))
  good_size_clades[[ii]] <- ancestry %>% filter(bootstrap == max(ancestry$bootstrap))
}
G_S_C <-vector()
for (jj in seq_along(good_size_clades)){G_S_C <- rbind(G_S_C,good_size_clades[[jj]])}
good_size_clades <- distinct(as_tibble(G_S_C))
good_size_clades %>% filter(bootstrap <70) %>% arrange(N_tips) %>% print(n=50) #12 clades with low support

##Collect offspring nodes
offs_pool<-vector()
a <- good_size_clades$node
for (nd in a){offs_pool <- c(offs_pool,offspring(x,nd)$node)}
offs_pool <- unique(offs_pool)
length(which(offs_pool <= length(tree$tip.label))) ## these are tips covered by nodes in a

###Remove nesting nodes, check total number of tips
partition <- x[a[which(a %ni% offs_pool)],] ##subsets tree x by nodes in a that are not in their own offspring pool
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

### Add the small clades to the partition
partition <-rbind(partition,compl)

##check that the assignment is unique
a <- partition$node
offs_pool<-vector()
for (nd in a){offs_pool <- c(offs_pool,offspring(x,nd)$node)}
offs_pool <- unique(offs_pool)

partition ###Current partition
x[a[which(a %ni% offs_pool)],] ###Partition checked for nesting clades
partition <- x[a[which(a %ni% offs_pool)],]
partition %>% select(N_tips) %>% sum() ###Check if all tips are there (7818)
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

########################################################
###Export partition--------
# write_delim(partition,"Partition.tsv", delim = "\t")

###Export node annotation for iTol------
# annotation <- partition %>%select(label,bootstrap) %>% mutate(position = 0.5, color = "rgb(0, 0, 0)"	,style = "normal", size = 	2	)
# write_tsv(annotation, "Node_annotation.txt",col_names = F)

# ###Export text files for every label in partition40 and populate with properly formatted gene id's for automatic retreaval---------
# for (n in 1:(nrow(partition))) {
#   clade <- partition[n,]$label
#   node <- partition[n,]$node
#   tips <- offspring(x,node, tiponly = T, self_include = T)
#   tipnames <- unlist(strsplit(tips$label,"/",fixed = T))[2*1:nrow(tips)-1]
#   write_delim(x = as.data.frame(tipnames), path = paste0("./Autoclades_70/",clade, "_",length(tips$label),".txt"), delim = "\t",quote_escape = "double",append = F,col_names = F)
# }
# 