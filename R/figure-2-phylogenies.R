# Ruan van Mazijk, 2020

library(phytools)  # Importing trees
library(tidyverse)
library(ggtree)    # Multi-phylo plots
# (Intalled with BiocManager::install("ggtree"))

tree <- read.tree("data/phylogenies/2020-07-11_RAxML-HPC-reconstruction_03/RAxML_bipartitions.result")

tree <- read.tree("data/phylogenies/2020-07-14_RAxML-HPC-reconstruction_04/RAxML_bipartitions.result")

Schoenus <- tree %>%
  drop.tip(.$tip.label[!str_detect(.$tip.label, "Schoenus")]) %>%
  ladderize()

Schoenus$tip.label <- str_replace(Schoenus$tip.label, "Schoenus_", "S. ")
Schoenus$node.label <- as.numeric(Schoenus$node.label)

Schoenus_BS_plot <- ggtree(Schoenus) +
  xlim(0, 0.33) +
  geom_tiplab(aes(label = paste0('italic(\"', label, '\")')), parse = TRUE, size = 2.5) +
  geom_nodepoint(aes(fill = as.numeric(label)), pch = 21, size = 2.5) +
  scale_fill_distiller(
    name      = "BS (%)",
    palette   = "RdYlGn",
    direction = 1,
    na.value  = "white"
  ) +
  theme(legend.position = c(0.1, 0.9))

Schoenus_simpler <- Schoenus
Schoenus_simpler$node.label <- Schoenus_simpler$node.label %>%
  as.numeric() %>%
  {case_when(
    . == 100 ~   "100",
    . >=  75 ~ ">= 75",
    TRUE     ~  "< 75"
  )}

Schoenus_simpler_BS_plot <- ggtree(Schoenus_simpler) +
  xlim(0, 0.33) +
  geom_tiplab(aes(label = paste0('italic(\"', label, '\")')), parse = TRUE, size = 2.5) +
  geom_nodepoint(aes(fill = factor(label, levels = c("100", ">= 75", "< 75"))), pch = 21, size = 2) +
  scale_fill_grey(name = "BS (%)", start = 0, end = 1) +
  theme(legend.position = c(0.1, 0.9))

tree_sample <- read.tree("data/phylogenies/2020-07-14_RAxML-HPC-reconstruction_04/RAxML_bootstrap.result")

set.seed(1234)
tree_sample <- sample(tree_sample, 100)

Schoenus_sample <- tree_sample
for (i in 1:length(tree_sample)) {
  Schoenus_sample[[i]] <- tree_sample[[i]] %>%
    drop.tip(.$tip.label[!str_detect(.$tip.label, "Schoenus")]) %>%
    compute.brlen()
  Schoenus_sample[[i]]$tip.label <- str_replace(Schoenus_sample[[i]]$tip.label, "Schoenus_", "S. ")
}

Schoenus_multitree_plot <- ggdensitree(Schoenus_sample, alpha = 0.05, jitter = 0) +
  xlim(0, 1.25) +
  geom_tiplab(aes(label = paste0('italic(\"', label, '\")')), parse = TRUE, size = 2.5)

ggsave("figures/Schoenus_BS_plot.pdf", Schoenus_BS_plot, width = 5, height = 10)
ggsave("figures/Schoenus_simpler_BS_plot.pdf", Schoenus_simpler_BS_plot, width = 5, height = 10)
ggsave("figures/Schoenus_multitree_plot.pdf", Schoenus_multitree_plot, width = 5, height = 10)
