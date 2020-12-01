library(tidyverse)

setwd("~/Box/NLR_binding_site_prediction/figures/Figure_RPP1/WSB_Entropy")


## Import data from Chimera printouts
distances <- read_delim("RR_distances.tsv",col_names = F, comment = '#', delim = "\t")
names <- read_delim("currentres.txt", col_names = F, delim  = "\t")
colnames(distances) <- names
nrow(distances)
ncol(distances)
nrow(names)
distances
names %>% print(n=100000)

## Format data

names %>% mutate(name = str_remove(names$X1,"..$"))->names
names %>% mutate(Residue = map(str_split(names$name, " "),2)%>%unlist%>%as.integer, AA = map(str_split(names$name, " "),1)%>%unlist) -> names
names %>% mutate(Protein = ifelse(Residue < 400, "ATR1","RPP1")) %>% select(Residue, AA, Protein)-> names
names %>% filter(Protein =="ATR1") %>% print(n=1000)
names %>% mutate(Index = 1:nrow(names))-> names


RPP1_res <- names %>% filter(Protein =="RPP1")
RPP1_res %>% print(n=1000)


find_min_dist <- function(x){
  return(min(distances[x,774:1010]))
}

RPP1_res %>% mutate(ATR_Dist = lapply(Index,find_min_dist)%>%unlist) -> RPP1_res
RPP1_res %>% filter(ATR_Dist < 5) %>% print(n=100)
RPP1_res %>% ggplot(aes(x=Residue, y=1/ATR_Dist))+geom_point()
RPP1_res %>% ggplot(aes(x=1/ATR_Dist))+geom_histogram()
1/.2
## Conclusion - these are not inter-atomic distances. i guess these are centroid or some such

install.packages("bio3d", dependencies=TRUE)
library(bio3d)
pdb <- read.pdb("../RPP1-ATR1-Final.pdb")
pdb <- read.pdb("minimal.pdb")
?dist.xyz
?dm
dm_minimal <- dm(pdb,ncore = 8, mask.lower = F)
dm_minimal
filled.contour(dm_minimal, nlevels = 10)
dm_minimal[1:773,3]
distances <- dm_minimal
?dm


entropy <- read_delim("6981_T283-R1_.ChimeraEntropy.txt",delim = "\t",skip = 3, col_names = c("Empty","Res","Entropy"))
entropy %>% mutate(Residue = as.integer(str_remove(Res,":")))%>%select(Residue,Entropy) -> entropy
left_join(RPP1_res,entropy) -> RPP1_res
RPP1_res %>% ggplot(aes(x=Entropy, y= 1/ATR_Dist))+geom_point()

RPP1_res %>% filter(Entropy > 1.5) %>% ggplot(aes(x=ATR_Dist))+geom_histogram()


dist_cut = 6
Entropy <- 0+1:30*.1

get_recall <- function(x){
  ent_cut <- x 
  tp <- RPP1_res %>% filter(ATR_Dist <= dist_cut, Entropy >= ent_cut) %>% nrow
  fp <- RPP1_res %>% filter(ATR_Dist > dist_cut, Entropy >= ent_cut) %>% nrow
  precision = tp/(tp+fp)
  fn <- RPP1_res %>% filter(ATR_Dist <= dist_cut, Entropy < ent_cut) %>% nrow
  recall <- tp/(tp+fn)
  return(recall)
}
get_precision <- function(x){
  ent_cut <- x 
  tp <- RPP1_res %>% filter(ATR_Dist <= dist_cut, Entropy >= ent_cut) %>% nrow
  fp <- RPP1_res %>% filter(ATR_Dist > dist_cut, Entropy >= ent_cut) %>% nrow
  precision = tp/(tp+fp)
  fn <- RPP1_res %>% filter(ATR_Dist <= dist_cut, Entropy < ent_cut) %>% nrow
  recall <- tp/(tp+fn)
  return(precision)
}

Table <- tibble(Entropy = Entropy)
Table %>% mutate(Recall = lapply(Entropy, get_recall)%>%unlist, Precision = lapply(Entropy, get_precision)%>%unlist) -> Table
Table <- Table %>% filter(Entropy <=2.3)
Table %>% ggplot(aes(x= Precision, y= Recall, color = Entropy))+geom_point()
Table %>% ggplot(aes(x= -Entropy, y= Recall, color = Precision))+geom_point()
Table %>% ggplot(aes(x= Entropy, y= Precision, color = Recall))+geom_point()

Table %>% filter(Entropy>.6, Entropy <1.9) %>% ggplot(aes(x= Precision, y= Recall, label = Entropy))+
  geom_point()+geom_text(hjust = 0, nudge_x = 0.005)+ labs(title = "Precision vs Recall at shown Entropy cutoffs")

Table %>% print(n=100)
