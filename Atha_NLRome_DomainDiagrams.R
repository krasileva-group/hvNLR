## ---------------------------
##
## Script name: Atha_NLRome_DomainDiagrams.R
##
## Purpose of script:
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

## set working directory

setwd("~/Downloads/Atha_NLRome_Domains/")  
domains_wide <- read_delim("Ath.tsv", delim = "\t", col_names = c("Gene", paste0("Domain_",as.character(1:7))))  
dom_long <- gather(domains_wide, value = "Domain", key = "Key",paste0("Domain_",as.character(1:7)), factor_key=TRUE) %>% select(-Key)
domains <- dom_long %>% mutate(
                    Dom = str_remove(Domain,"\\(.*$") %>% as.factor(),
                    Start = str_remove(Domain,"^.*start=") %>% str_remove(", stop.*$") %>%as.integer(),
                    Stop = str_remove(Domain,"^.*stop=") %>% str_remove(", eval.*$") %>%as.integer(),
                    Eval = str_remove(Domain,"^.*evalue=") %>% str_remove("\\)") %>%as.numeric()
                    ) %>% select(-Domain)
domains %>% select(Gene) %>% distinct() %>% arrange(Gene)
write_delim(domains,delim = "\t", col_names = TRUE, path = "Atha_NLRome_Domains.tsv")
domains %>% filter(!is.na(Dom)) ->domains
domains %>% group_by(Dom) %>% summarise(n=n()) %>% arrange (n) %>% print(n=300)
domains %>% group_by(Dom,Gene) %>% summarise(n=n()) %>% arrange (n) %>% filter (n>1, !grepl("LRR", Dom), Dom != "NB-ARC", Dom != "TIR") %>% print(n=300)

### iTOL export----
# need to assign a color and a shape to every domain
# this can then be joined to the domains table and subset on the leaves in the tree
# finally sink+cat can be used to export iTOL-readable annotation file

### find all domains, order by prevalence
domains %>% group_by(Dom) %>% summarise(n=n()) %>% arrange (n) %>% print(n=300)
### 177 domains right now.
dist_doms <- domains %>% group_by(Dom) %>% filter (!grepl(x = Dom, pattern = "LRR_")) %>% summarise(n=n()) %>% arrange (n) %>% print(n=300)
### 169 excluding LRR (LRRs from LRR predictor can be added later)
### RPPW, Rx_N, TIR, and NB-ARC are the most popular, The first three should share shape and have different colors
### NB-ARC and LRR should have different colors and shapes
### NB-ARC - yellow hor hexagon
### LRR - red rectangle
### Post-LRR - red diamond
### TIR - orange ellipse
### RPP8 - purple ellipse
### Rx_n - dark green ellipse
### Others  - randomly assign 12 divergent colors and shapes from the list: HV, TR, TL, PL, PR, PU, PD, OC
shapes <- c("HV", "TR", "TL", "PL", "PR", "PU", "PD", "OC")
sample(shapes, 165, replace = T)
#
### 12 Colors to use from Color Brewer----
colors <- c("#a6cee3",
  "#1f78b4",
  "#b2df8a",
  "#33a02c",
  "#fb9a99",
  "#e31a1c",
  "#fdbf6f",
  "#ff7f00",
  "#cab2d6",
  "#6a3d9a",
  "#ffff99",
  "#b15928")
sample(colors, 165, replace = T)

#
### Building the color tables ---------
color_doms <- tibble(Dom = dist_doms$Dom[1:165], 
                     Shape = sample(shapes, 165, replace = T),
                     Color = sample(colors, 165, replace = T)
                     )
color_doms
main_doms <- tibble(Dom = c("NB-ARC", "RPW8", "Rx_N", "TIR"), 
                    Shape = c("HH", "EL", "EL", "EL"),
                    Color = c("#ffff99", "#6a3d9a", "#ff7f00", "#33a02c")
                    )

lrr_doms <- tibble(Dom = c("LRR_1", "LRR_2", "LRR_3", "LRR_4", "LRR_5", "LRR_6", "LRR_8", "LRR_9"), 
                   Shape = rep_len("RE",8),
                   Color = rep_len("#e31a1c",8)
                    )

all_doms <- rbind(color_doms,main_doms,lrr_doms)
## got 177 combinations of domains and shapes/colors 
write_delim(all_doms, "All_Domains_with_colors.tsv", delim = "\t", col_names = T)

#
#### Domain shapes for reference ----
#RE  rectangle
#HH  horizontal hexagon
#HV  vertical hexagon
#EL  ellipse
#DI  rhombus (diamond)
#TR  right pointing triangle
#TL  left pointing triangle
#PL  left pointing pentagram
#PR  right pointing pentagram
#PU  up pointing pentagram
#PD  down pointing pentagram
#OC  octagon
#GP  rectangle (gap; black filled rectangle with 1/3 normal height)


#
### Export to iTOL ------------
all_color_doms <- left_join(domains, all_doms) 
all_color_doms <- all_color_doms %>% mutate(iTOL = paste(Shape, Start, Stop, Color, Dom, sep = "|"))
### Need to leftjoin a table of protein lengths to this and will be ready for export!!!
### For the original tree will also need to bind names/start-stop to match leaf IDs
setwd("~/BioInf/ArabidopsisRENSEQ/Proteomes/")
getwd()
files <- list.files(pattern = ".fa")
require(Biostrings)
lengths<-vector("list",length = length(files))
ii<-1
for (ii in seq_along(files)){
  print(paste0("Looking at file ",ii, " of ", length(files), " (",files[[ii]],")..."))
  a<-readAAStringSet(files[[ii]])
  lengths[[ii]] <- tibble(Gene = a@ranges@NAMES, Length = a@ranges@width)
}
rm(a,b,c)

prot_l<-lengths[[1]]
for (jj in 2:length(lengths)){prot_l<-rbind(prot_l,lengths[[jj]])}
prot_l %>% mutate(Gene = str_remove(Gene, " ")) ->prot_l
prot_l %>% mutate(Gene = str_remove(Gene,"pacid.*$")) -> prot_l

all_color_doms <- left_join(all_color_doms,prot_l) %>% filter(!is.na(Length)) %>% print(n=1000)
all_color_doms %>% filter(grepl(x = Gene,pattern = "AT"))


setwd("~/Downloads/Atha_NLRome_Domains/")  
genreg <- read_csv("Leaves.txt", col_names = "GeneRegion")
genreg %>% mutate(Gene = str_remove(GeneRegion,"\\/.*$")) ->genreg
left_join(genreg,all_color_doms) %>% mutate(Gene = GeneRegion) %>% select(-GeneRegion) ->genreg

export<-all_color_doms
export<-genreg
sink("genreg_domains.txt", append = F)
cat(
"DATASET_DOMAINS
SEPARATOR COMMA
DATASET_LABEL,Domains
COLOR,#ff0000
DATA
")
for (gene in levels(as.factor(export$Gene))){
  tbl <- export %>% filter(Gene == gene)
  cat(paste(gene, tbl[[1,9]], sep = ","))
  
  for (ii in seq_along(tbl$Dom)){
    cat(",")
    cat(tbl[[ii,8]])
      }
  cat("\n")
}
sink()
