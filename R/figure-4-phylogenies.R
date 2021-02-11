# Ruan van Mazijk, 2021

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
library(lemon)      # For ::coord_capped_cart()

# Import data ------------------------------------------------------------------

MCC_tree <- treeio::read.beast("data/phylogenies/Cyperaceae-all-taxa-6-calib-st6-Dec26.tre")
posterior_sample <- read.tree("data/phylogenies/Schoenus.trees.100.trees")

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
  ladderize(right = FALSE) %>%
  force.ultrametric(method = "extend")

# Tidy tip labels
Schoenus_MCC@phylo$tip.label <- str_replace(
  Schoenus_MCC@phylo$tip.label,
  "Schoenus_", "S. "
)

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
  ggtree(Schoenus_MCC, ladderize = FALSE) +  # (already ladderized above!)
  geom_rootedge(rootedge = 5) +
  geom_tiplab(
    aes(label = paste0('italic(\"', label, '\")')),
    parse = TRUE,
    size = 2.5,
    offset = 2
  ) +
  geom_nodepoint(aes(fill = posterior), pch = 21, size = 2.5) +
  scale_fill_gradient(name = "PP",
    na.value  = "white", low = "white", high = "darkgreen",
    limits = c(0.9, 1),
    breaks = c(0.9, 0.95, 1.0),
    labels = c(expression(phantom(x) < "0.90"), "0.95", "1.00")
  ) +
  theme_tree2() +
  scale_x_continuous(name = "Ma",
    limits = c(-5, 68),
    breaks = label_positions,
    labels = -my_labels
  ) +
  # Remove empty space above, below tree
  scale_y_continuous(
    limits = c(0, Ntip(Schoenus_MCC@phylo) + 1),
    expand = c(0, 0)
  ) +
  # Remove extra line at right of time axis
  coord_capped_cart(bottom = "right") +
  theme(
    legend.position   = c(0.08, 0.80),
    legend.text.align = 1,
    legend.background = element_blank()
  )

# X-axis scaling things:
max_tree_height <- Schoenus_posterior %>%
  map_dbl(~max(nodeHeights(.))) %>%
  max()
label_positions <- max_tree_height - c(80, 70, 60, my_labels)

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
    labels = -c(80, 70, 60, my_labels)
  ) +
  # Remove empty space above, below tree
  scale_y_continuous(
    limits = c(0, Ntip(Schoenus_MCC@phylo) + 1),
    expand = c(0, 0)
  ) +
  # Remove extra line at left of time axis
  coord_capped_cart(bottom = "left") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

# Panel these pltos together using patchwork::
Schoenus_tree_plots <- Schoenus_MCC_plot + Schoenus_posterior_plot +
  plot_annotation(tag_levels = "A")

# Save plots -------------------------------------------------------------------

ggsave(
  "figures/Schoenus-phy-bio-Fig4-RvM.pdf",
  Schoenus_tree_plots,
  width = 10, height = 12
)

ggsave(
  "figures/Schoenus-phy-bio-Fig4-RvM.png",
  Schoenus_tree_plots,
  width = 10, height = 12, dpi = 300
)
