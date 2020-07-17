# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(tidyverse)
library(phytools)  # Importing trees
library(ggtree)    # Multi-phylo plots
                   # (Intalled with BiocManager::install("ggtree"))
library(jntools)   # For ::get_tips_in_ape_plot_order()
                   # (Installed with remotes::install_github("joelnitta/jntools"))
library(patchwork)  # Figure panelling

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
  xlim(0, 0.31) +
  geom_tiplab(
    label = "",
    size = 2.5,
    align = TRUE,
    colour = "grey50",
  ) +
  geom_nodepoint(aes(fill = as.numeric(label)), pch = 21, size = 2.5) +
  scale_fill_gradient(name = "BS (%)",
    na.value  = "white", low = "white", high = "black",
    limits = c(50, 100),
    labels = c("< 50", "60", "70", "80", "90", "100")
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    legend.position = c(0.1, 0.15), legend.text.align = 1,
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )

Schoenus_multitree_plot <-
  ggdensitree(Schoenus_sample,
    alpha     = 0.05,
    tip.order = get_tips_in_ape_plot_order(Schoenus)
  ) +
  geom_tiplab(
    aes(label = paste0('italic(\"', label, '\")')),
    parse = TRUE,
    size = 2.5,
    offset = -0.25
  ) +
  scale_x_reverse(expand = c(0, 0)) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

Schoenus_RAxML_plot <- Schoenus_BS_plot + Schoenus_multitree_plot

# Save plot -------------------------------------------------------------------

ggsave(
  "figures/Schoenus_RAxML_plot.pdf",
  Schoenus_RAxML_plot,
  width = 10, height = 10
)
