# Ruan van Mazijk, 2020

library(phytools)  # Importing trees
library(stringr)   # String manipulation
library(ggtree)    # Multi-phylo plots
# (Intalled with BiocManager::install("ggtree"))

tree <- read.tree("data/phylogenies/2020-07-11_RAxML-HPC-reconstruction_03/RAxML_bestTree.result")
#tree <- read.tree("data/phylogenies/2020-07-11_RAxML-HPC-reconstruction_03/RAxML_bipartitions.result")
trees <- read.tree("data/phylogenies/2020-07-11_RAxML-HPC-reconstruction_03/RAxML_bootstrap.result")

Schoenus <- tree %>%
  drop.tip(.$tip.label[!str_detect(.$tip.label, "Schoenus")]) %>%
  ladderize()

Schoenus$tip.label <- str_replace(Schoenus$tip.label, "Schoenus", "S.")

plotTree(Schoenus, fsize = 0.75, ftype = "i")

ggdensitree(trees, alpha = 0.1) +
  geom_tiplab(size = 2)
