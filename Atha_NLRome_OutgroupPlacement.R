a<-read.jplace("/Users/prigozhin/BioInf/ArabidopsisRENSEQ/Phylogeny/RAxML_NLRome_237min_Outgroup/EPA_z/epa_result.jplace")
big_table <- left_join(as_tibble(a@phylo),a@placements)%>% filter(!is.na(name))
# b<-read_delim("/Users/prigozhin/BioInf/ArabidopsisRENSEQ/Phylogeny/RAxML_NLRome_237min_Outgroup/EPA_z/Int_to_EPA.tsv","\t",col_names = c("label","Number"))
# b<-b[1:7815,]
# left_join(big_table,b)%>%print(n=1000)

### limit placements to >95% likelihood ratio cutoff
#placement <- big_table %>% filter(like_weight_ratio>.95) %>% arrange(name) %>% print(n=200)
placement <- big_table

### import NLRome tree rooted in iTOL
placement
tree_pl <- left_join(x,placement)

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
partition %>% filter(label %ni% clade_og$Clade)
partition %>% filter(label %ni% clade_og$Clade) %>% select(N_tips) %>% sum() ## 635 tips in initial clades without Capsella or lyrata sequences

clade_og %>% group_by(Clade) %>% summarise(n=n()) %>% arrange(Clade) -> og_summary
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
