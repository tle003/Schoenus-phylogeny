# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(tidyverse)
library(phytools)  # Importing trees
library(ggtree)    # Multi-phylo plots
                   # (Intalled with BiocManager::install("ggtree"))
library(jntools)   # For ::get_tips_in_ape_plot_order()
                   # (Installed with remotes::install_github("joelnitta/jntools"))

# Import data ------------------------------------------------------------------

# RAxML-HPC reconstruction:
# Best tree with nodes' bootstrap support values
tree <- read.tree("data/phylogenies/2020-07-14_RAxML-HPC-reconstruction_04/RAxML_bipartitions.result")
# Bootstrap sample (1000 + 8 trees)
tree_sample <- read.tree("data/phylogenies/2020-07-14_RAxML-HPC-reconstruction_04/RAxML_bootstrap.result")

# Tidy data --------------------------------------------------------------------

# RAxML-HPC reconstruction:

# Best tree:
# Extract Schoenus
Schoenus <- tree %>%
  drop.tip(.$tip.label[!str_detect(.$tip.label, "Schoenus")]) %>%
  ladderize()
# Tidy tip labels
Schoenus$tip.label <- str_replace(
  Schoenus$tip.label,
  "Schoenus_", "S. "
)
# Tidy node lables
Schoenus$node.label <- as.numeric(Schoenus$node.label)

# Bootstrap sample:
# Sub-sample 100 tree
set.seed(1234)
tree_sample <- sample(tree_sample, 100)
Schoenus_sample <- tree_sample
for (i in 1:length(tree_sample)) {
  Schoenus_sample[[i]] <- tree_sample[[i]] %>%
    # Extract Schoenus from each tree
    drop.tip(.$tip.label[!str_detect(.$tip.label, "Schoenus")]) %>%
    # Ultrametricise each tree and make tree depth 1
    compute.brlen()
  # Tidy tip labels
  Schoenus_sample[[i]]$tip.label <- str_replace(
    Schoenus_sample[[i]]$tip.label,
    "Schoenus_", "S. "
  )
}

# Plots ------------------------------------------------------------------------

root_length <- 0.01

Schoenus_BS_plot <-
  ggtree(Schoenus,
    ladderize     = TRUE,
    right         = TRUE,
    root.position = root_length
  ) +
  geom_rootedge(rootedge = root_length) +
  xlim(0, 0.33) +
  geom_tiplab(
    aes(label = paste0('italic(\"', label, '\")')),
    parse = TRUE,
    size = 2.5
  ) +
  geom_nodepoint(aes(fill = as.numeric(label)), pch = 21, size = 2.5) +
  scale_fill_gradient(name = "BS (%)",
    na.value  = "white", low = "white", high = "black",
    limits = c(50, 100),
    labels = c("< 50", "60", "70", "80", "90", "100")
  ) +
  theme(legend.position = c(0.9, 0.15), legend.text.align = 1)

Schoenus_multitree_plot <-
  ggdensitree(Schoenus_sample,
    alpha     = 0.05,
    tip.order = get_tips_in_ape_plot_order(Schoenus)
  ) +
  xlim(0, 1.25) +
  geom_tiplab(
    aes(label = paste0('italic(\"', label, '\")')),
    parse = TRUE,
    size = 2.5
  )

# Save plots -------------------------------------------------------------------

ggsave(
  "figures/Schoenus_BS_plot.pdf",
  Schoenus_BS_plot,
  width = 6, height = 12
)
ggsave(
  "figures/Schoenus_multitree_plot.pdf",
  Schoenus_multitree_plot,
  width = 5, height = 10
)
