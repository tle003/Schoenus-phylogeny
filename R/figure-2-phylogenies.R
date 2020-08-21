# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(tidyverse)
library(phytools)   # Importing trees
library(phyloch)    # For ::read.beast()
                    # (Installed with devtools::install_github("fmichonneau/phyloch"))
library(ggtree)     # Multi-phylo plots
                    # (Installed with BiocManager::install("ggtree"))
library(jntools)    # For ::get_tips_in_ape_plot_order()
                    # (Installed with remotes::install_github("joelnitta/jntools"))
library(patchwork)  # Figure panelling

# Import data ------------------------------------------------------------------

MCC_tree <- treeio::read.beast("data/phylogenies/Cyperaceae-all-taxa-6calib-max-clad-AUG12.tre")
posterior_sample <- read.tree("data/phylogenies/Cyperaceae.trees.100.trees")

# Tidy data --------------------------------------------------------------------

# MCC tree:
# Extract Schoenus
Schoenus_MRCA_node <- MCC_tree@phylo %>%
  MRCA(.$tip.label[str_detect(.$tip.label, "Schoenus")])
# (Look at node numbers to confirm correct node recovered)
#plotTree(ladderize(MCC_tree@phylo), fsize = 1, node.numbers = TRUE)
Schoenus_MCC <- MCC_tree %>%
  treeio::tree_subset(Schoenus_MRCA_node, levels_back = 0)
Schoenus_MCC@phylo <- Schoenus_MCC@phylo %>%
  ladderize(right = TRUE) %>%
  force.ultrametric(method = "extend")

# Tidy tip labels
Schoenus_MCC@phylo$tip.label <- str_replace(
  Schoenus_MCC@phylo$tip.label,
  "Schoenus_", "S. "
)

# Get node numbers for clades to highlight
Clade_A_node <- getMRCA(Schoenus_MCC@phylo, c(
  "S. insolitus",
  "S. sculptus"
))
Clade_B_node <- getMRCA(Schoenus_MCC@phylo, c(
  "S. falcatus",
  "S. australis"
))
Cape_clade_node <- getMRCA(Schoenus_MCC@phylo, c(
  "S. dregeanus",
  "S. australis"
))

# Posterior sample (already thinned to 100 trees):
Schoenus_posterior <- list(length = length(posterior_sample))
for (i in seq_along(posterior_sample)) {
  # Extract Schoenus from each tree
  Schoenus_posterior[[i]] <- posterior_sample[[i]] %>%
    drop.tip(.$tip.label[!str_detect(.$tip.label, "Schoenus")]) %>%
    force.ultrametric(method = "extend")
  # Tidy tip labels
  Schoenus_posterior[[i]]$tip.label <- str_replace(
    Schoenus_posterior[[i]]$tip.label,
    "Schoenus_", "S. "
  )
}
# Convert list of pruned trees to multiphylo
Schoenus_posterior_multiphylo <- Schoenus_posterior[[1]]
for (i in 2:length(Schoenus_posterior)) {
  Schoenus_posterior_multiphylo <- c(
    Schoenus_posterior_multiphylo,
    Schoenus_posterior[[i]]
  )
}
Schoenus_posterior <- Schoenus_posterior_multiphylo

# Plots ------------------------------------------------------------------------

# X-axis scaling things:
tree_height <- max(nodeHeights(Schoenus_MCC@phylo))
my_labels <- c(50, 40, 30, 20, 10, 0)
label_positions <- tree_height - my_labels

Schoenus_MCC_plot <-
  ggtree(Schoenus_MCC) +
  geom_tiplab(
    aes(label = paste0('italic(\"', label, '\")')),
    parse = TRUE,
    size = 2.5,
    offset = 2
  ) +
  geom_hilight(Clade_A_node,    fill = "black", alpha = 0.1) +
  geom_hilight(Clade_B_node,    fill = "black", alpha = 0.2) +
  geom_hilight(Cape_clade_node, fill = "black", alpha = 0.2) +
  geom_nodepoint(aes(fill = posterior), pch = 21, size = 2.5) +
  scale_fill_gradient(name = "PP",
    na.value  = "white", low = "white", high = "darkgreen",
    limits = c(0.5, 1),
    labels = c("<= 0.5", "0.6", "0.7", "0.8", "0.9", "1.0")
  ) +
  theme_tree2() +
  scale_x_continuous(name = "Ma",
    limits = c(0, 68),
    breaks = label_positions,
    labels = my_labels
  ) +
  theme(
    legend.position   = c(0.150, 0.875),
    legend.text.align = 1,
    plot.margin       = unit(c(0, 0, 0, 0), "cm")
  )

# X-axis scaling things:
max_tree_height <- Schoenus_posterior %>%
  map_dbl(~max(nodeHeights(.))) %>%
  max()
label_positions <- max_tree_height - my_labels

Schoenus_posterior_plot <-
  ggdensitree(Schoenus_posterior,
    alpha     = 0.03,
    tip.order = rev(get_tips_in_ape_plot_order(ladderize(Schoenus_MCC@phylo)))
    # NOTE: Don't why this works, but it's the only way to get the 2 panels
    # to have the same tip-order...
    # Must be something to do with ape, phytools and ggtree ladderizing in
    # slightly different ways.
  ) +
  theme_tree2() +
  scale_x_reverse(name = "Ma",
    breaks = label_positions,
    labels = my_labels
  ) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

Schoenus_tree_plots <- Schoenus_MCC_plot + Schoenus_posterior_plot

# Save plots -------------------------------------------------------------------

ggsave(
  "figures/Schoenus_tree_plots.pdf",
  Schoenus_tree_plots,
  width = 10, height = 12
)

ggsave(
  "figures/Schoenus_tree_plots.png",
  Schoenus_tree_plots,
  width = 10, height = 12, dpi = 300
)
