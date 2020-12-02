## ggtree demo
# modified by Shoshana Williams, Dec 1, 2020

#install.packages("BiocManager"), then instsall all necessary libraries
#BiocManager::install("ggtree")
library(ggtree)
library(ape)
library(ggplot2) # ggplot() for plotting
library(tidyr) # data reformatting
library(stringr) # string manipulation
set.seed(2017-02-16)

# set working directory
pth <- "/Users/ShoshanaWilliams/Desktop/Austin Lab/final_fused-rd_tree/"
wd_pth <- pth
setwd(pth)

# load in metadata and newick tree
tree_pth <- paste(pth,"RAxML_bestTree.output.tree", sep = "")
tree <- read.tree(tree_pth) 
info_pth <- paste(pth,"curated_RDS_genus.csv",sep="")
info <- read.csv(info_pth)

# set colors for info.csv category 'location'
p <- ggtree(tree, layout='circular', size = 0.5, branch.length = 'none') %<+% info + 
    geom_treescale(x=0, y=60, width=1, fontsize=0) +
    geom_tiplab(aes(label=name, color=genus), align=T,
             size=2, offset=0, hjust=0) + xlim(0, 300)
   ggsave('test.pdf', width = 50, height =50, limitsize = FALSE)
p

