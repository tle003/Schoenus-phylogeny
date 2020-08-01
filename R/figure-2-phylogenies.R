# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(tidyverse)
library(phytools)   # Importing trees
library(ggtree)     # Multi-phylo plots
                    # (Installed with BiocManager::install("ggtree"))
library(jntools)    # For ::get_tips_in_ape_plot_order()
                    # (Installed with remotes::install_github("joelnitta/jntools"))
library(patchwork)  # Figure panelling

# Import data ------------------------------------------------------------------

# RAxML-HPC reconstruction:
# Best tree with nodes' bootstrap support values
RAxML_tree <- read.tree("data/phylogenies/2020-07-14_RAxML-HPC-reconstruction_04/RAxML_bipartitions.result")
# Bootstrap sample (1000 + 8 trees)
BS_trees <- read.tree("data/phylogenies/2020-07-14_RAxML-HPC-reconstruction_04/RAxML_bootstrap.result")

# BEAST reconstruction:
MCC_tree <- read.nexus("data/phylogenies/2020-07-29_BEAST-reconstruction/Cyperaceae-all-taxa-6-calib-comb-29JUL.tre")
posterior_sample <- read.nexus("data/phylogenies/2020-07-29_BEAST-reconstruction/Cyperaceae-all-taxa-6-calib-comb-29JUL-thinned.trees")

# Tidy data --------------------------------------------------------------------

# .... RAxML-HPC reconstruction ------------------------------------------------

# Best tree:
# Extract Schoenus
Schoenus_RAxML <- RAxML_tree %>%
  drop.tip(.$tip.label[!str_detect(.$tip.label, "Schoenus")]) %>%
  ladderize()
# Tidy tip labels
Schoenus_RAxML$tip.label <- str_replace(
  Schoenus_RAxML$tip.label,
  "Schoenus_", "S. "
)
# Tidy node labels
Schoenus_RAxML$node.label <- as.numeric(Schoenus_RAxML$node.label)

####
library(treeio)
well_supported_nodes <- Schoenus %>%
  as.treedata() %>%
  as_tibble() %>%
  filter(!str_detect(label, "S")) %>%
  mutate(BS = as.numeric(label)) %>%
  select(parent, node, BS) %>%
  filter(BS >= 80, parent != node) %>%
  mutate(n_spp = map_int(node, ~length(offspring(Schoenus, ., tiponly = TRUE)))) %>%
  filter(n_spp < 40, n_spp > 3) %>%
  #filter(!(parent %in% map(node, ~offspring(Schoenus, .)))) %>%
  #filter(
  #  map2_lgl(parent, node, ~ .y %in% offspring(Schoenus, .x))
  #) %>%
  pull(node)
p <- ggtree(Schoenus, ladderize = TRUE, right = TRUE)
for (i in seq_along(well_supported_nodes)) {
  p <- p + geom_hilight(node = well_supported_nodes[[i]])
}
p
####

# Bootstrap sample:
# Sub-sample 100 trees
set.seed(1234)
BS_sample <- sample(BS_trees, 100)
Schoenus_BS_sample <- BS_sample
for (i in seq_along(BS_sample)) {
  Schoenus_BS_sample[[i]] <- BS_sample[[i]] %>%
    # Extract Schoenus from each tree
    drop.tip(.$tip.label[!str_detect(.$tip.label, "Schoenus")]) %>%
    # Ultrametricise each tree and make tree depth 1
    compute.brlen()
  # Tidy tip labels
  Schoenus_BS_sample[[i]]$tip.label <- str_replace(
    Schoenus_BS_sample[[i]]$tip.label,
    "Schoenus_", "S. "
  )
}

# .... BEAST reconstruction ----------------------------------------------------

# MCC tree:
# Extract Schoenus
Schoenus_MCC <- MCC_tree %>%
  drop.tip(.$tip.label[!str_detect(.$tip.label, "Schoenus")]) %>%
  ladderize()
# Tidy tip labels
Schoenus_MCC$tip.label <- str_replace(
  Schoenus_MCC$tip.label,
  "Schoenus_", "S. "
)

# Posterior sample:
Schoenus_posterior <- list(length = length(posterior_sample))
for (i in seq_along(posterior_sample)) {
  Schoenus_posterior[[i]] <- posterior_sample[[i]] %>%
    # Extract Schoenus from each tree
    drop.tip(.$tip.label[!str_detect(.$tip.label, "Schoenus")])
  # Tidy tip labels
  Schoenus_posterior[[i]]$tip.label <- str_replace(
    Schoenus_posterior[[i]]$tip.label,
    "Schoenus_", "S. "
  )
}
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
