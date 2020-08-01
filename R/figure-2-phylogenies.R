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

MCC_tree <- read.beast("data/phylogenies/2020-07-29_BEAST-reconstruction/Cyperaceae-all-taxa-6-calib-comb-29JUL.tre")
posterior_sample <- read.nexus("data/phylogenies/2020-07-29_BEAST-reconstruction/Cyperaceae-all-taxa-6-calib-comb-29JUL-thinned.trees")

# Tidy data --------------------------------------------------------------------

# MCC tree:
# Extract Schoenus
Schoenus_MRCA_node <- MCC_tree@phylo %>%
  MRCA(.$tip.label[str_detect(.$tip.label, "Schoenus")])
# (Look at node numbers to confirm correct node recovered)
#plotTree(ladderize(MCC_tree@phylo), fsize = 1, node.numbers = TRUE)
Schoenus_MCC <- MCC_tree %>%
  tree_subset(Schoenus_MRCA_node, levels_back = 0)
Schoenus_MCC@phylo <- ladderize(Schoenus_MCC@phylo, right = TRUE)

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
    drop.tip(.$tip.label[!str_detect(.$tip.label, "Schoenus")])
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

max_tree_height <- max(nodeHeights(Schoenus_MCC@phylo))
label_positions <- c(max_tree_height-50, max_tree_height-40, max_tree_height-30, max_tree_height-20, max_tree_height-10, max_tree_height)
Schoenus_MCC_plot <-
  ggtree(Schoenus_MCC) +
  geom_tiplab(
    aes(label = paste0('italic(\"', label, '\")')),
    parse = TRUE,
    size = 2.5,
    offset = 2
  ) +
  geom_nodepoint(aes(fill = posterior), pch = 21, size = 2.5) +
  scale_fill_gradient(name = "PP",
    na.value  = "white", low = "white", high = "darkgreen",
    limits = c(0.5, 1),
    labels = c("<= 0.5", "0.6", "0.7", "0.8", "0.9", "1.0")
  ) +
  theme_tree2() +
  scale_x_continuous(name = "Ma", limits = c(0, 68), breaks = label_positions, labels = c(50, 40, 30, 40, 10, 0)) +
  theme(
    legend.position   = c(0.150, 0.875),
    legend.text.align = 1,
    plot.margin       = unit(c(0, 0, 0, 0), "cm")
  )

max_tree_height <- max(purrr::map_dbl(Schoenus_posterior, ~max(nodeHeights(.))))
#mean_tree_height <- mean(purrr::map_dbl(Schoenus_posterior, ~max(nodeHeights(.))))
label_positions <- c(max_tree_height-50, max_tree_height-40, max_tree_height-30, max_tree_height-20, max_tree_height-10, max_tree_height)
Schoenus_posterior_plot <-
  ggdensitree(Schoenus_posterior,
    alpha     = 0.03,
    tip.order = rev(get_tips_in_ape_plot_order(ladderize(Schoenus_MCC@phylo)))
  ) +
  theme_tree2() +
  scale_x_reverse(name = "Ma", breaks = label_positions, labels = c(50, 40, 30, 40, 10, 0)) +
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
