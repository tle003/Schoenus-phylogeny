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

root_length <- 0.01

Schoenus_RAxML_plot <-
  ggtree(Schoenus_RAxML,
    ladderize     = TRUE,
    right         = TRUE,
    root.position = root_length
  ) +
  geom_rootedge(rootedge = root_length) +
  xlim(0, 0.31) +
  geom_tiplab(
    label = "",
    size = 2.5,
    align = TRUE,
    colour = "grey50",
  ) +
  geom_nodepoint(aes(fill = as.numeric(label)), pch = 21, size = 2.5) +
  scale_fill_gradient(name = "BS (%)",
    na.value  = "white", low = "white", high = "darkgreen",
    limits = c(50, 100),
    labels = c("<= 50", "60", "70", "80", "90", "100")
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    legend.position = c(0.1, 0.125), legend.text.align = 1,
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

Schoenus_posterior_plot <-
  ggdensitree(Schoenus_posterior,
    alpha     = 0.025,
    tip.order = get_tips_in_ape_plot_order(Schoenus_RAxML)
  ) +
  geom_tiplab(
    aes(label = paste0('italic(\"', label, '\")')),
    parse = TRUE,
    size = 2.5,
    offset = -20
  ) +
  scale_x_reverse(expand = c(0, 0)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

Schoenus_tree_plots <- Schoenus_RAxML_plot + Schoenus_posterior_plot

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
