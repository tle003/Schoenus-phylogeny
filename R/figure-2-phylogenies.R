# Ruan van Mazijk, 2020

library(phytools)  # Importing trees
library(tidyverse)
library(ggtree)    # Multi-phylo plots
# (Intalled with BiocManager::install("ggtree"))

#remotes::install_github("joelnitta/jntools")
library(jntools)  # ::get_tips_in_ape_plot_order()

tree <- read.tree("data/phylogenies/2020-07-14_RAxML-HPC-reconstruction_04/RAxML_bipartitions.result")

Schoenus <- tree %>%
  drop.tip(.$tip.label[!str_detect(.$tip.label, "Schoenus")]) %>%
  ladderize()

Schoenus$tip.label <- str_replace(Schoenus$tip.label, "Schoenus_", "S. ")
Schoenus$node.label <- as.numeric(Schoenus$node.label)

Schoenus_BS_plot <- ggtree(Schoenus, ladderize = TRUE, right = TRUE, root.position = 0.01) +
  geom_rootedge(rootedge = 0.01) +
  xlim(0, 0.33) +
  geom_tiplab(aes(label = paste0('italic(\"', label, '\")')), parse = TRUE, size = 2.5) +
  geom_nodepoint(aes(fill = as.numeric(label)), pch = 21, size = 2.5) +
  scale_fill_distiller(
    limits = c(50, 100),
    name      = "BS (%)",
    palette   = "Greys",
    direction = 1,
    na.value  = "white",
    labels = c("< 50", "60", "70", "80", "90", "100")
  ) +
  theme(legend.position = c(0.9, 0.15), legend.text.align = 1)

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

Schoenus_multitree_plot <- ggdensitree(Schoenus_sample, alpha = 0.05, jitter = 0, tip.order = get_tips_in_ape_plot_order(Schoenus_simpler)) +
  xlim(0, 1.25) +
  geom_tiplab(aes(label = paste0('italic(\"', label, '\")')), parse = TRUE, size = 2.5)

ggsave("figures/Schoenus_BS_plot.pdf", Schoenus_BS_plot, width = 5, height = 10)
ggsave("figures/Schoenus_simpler_BS_plot.pdf", Schoenus_simpler_BS_plot, width = 5, height = 10)
ggsave("figures/Schoenus_multitree_plot.pdf", Schoenus_multitree_plot, width = 5, height = 10)
