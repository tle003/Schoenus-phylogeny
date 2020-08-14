# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(tidyverse)
library(phytools)
library(ggtree)  # Multi-phylo plots
                 # (Installed with BiocManager::install("ggtree"))
library(treeio)  # For ::read.beast()

# Import data ------------------------------------------------------------------

MCC_tree <- read.beast("data/phylogenies/Cyperaceae-all-taxa-6calib-max-clad-AUG12.tre")

MCC_tree@phylo <- MCC_tree@phylo %>%
  force.ultrametric(method = "extend") %>%
  ladderize(right = TRUE)

Schoenus_MRCA_node <- MCC_tree@phylo %>%
  getMRCA(.$tip.label[str_detect(.$tip.label, "Schoenus")])

Schoeneae_MRCA_node <- MCC_tree@phylo %>%
  getMRCA(c(
    .$tip.label[str_detect(.$tip.label, "Schoenus")],
    "Gymnoschoenus_sphaerocephalus"
  ))

MCC_tree@phylo$tip.label <- MCC_tree@phylo$tip.label %>%
  str_replace("Schoenus_", "S._") %>%
  str_replace("_", " ")

# Plot -------------------------------------------------------------------------

# X-axis scaling things:
tree_height <- max(nodeHeights(MCC_tree@phylo))
my_labels <- c(90, 80, 70, 60, 50, 40, 30, 20, 10, 0)
label_positions <- tree_height - my_labels

nodes_with_future_HPDs <- MCC_tree@data %>%
  select(node, height_0.95_HPD) %>%
  unnest() %>%
  filter(height_0.95_HPD > tree_height) %>%
  pull(node)

MCC_tree@data %>%
  filter(node %in% nodes_with_future_HPDs) %>%
  unnest()

Cyperaceae_tree_plot <- ggtree(MCC_tree) +
  geom_tiplab(
    aes(label = paste0('italic(\"', label, '\")')),
    parse = TRUE,
    size = 2.5,
    offset = 10
  ) +
  geom_range("height_0.95_HPD",
    center = "height_median",
    size   = 1.5,
    alpha  = 0.33
  ) +
  #geom_hilight(Schoenus_MRCA_node,  fill = "darkblue",  alpha = 0.25) +
  #geom_hilight(Schoeneae_MRCA_node, fill = "blue", alpha = 0.25) +
  scale_x_continuous(name = "Ma",
    limits = c(0, 150),
    breaks = label_positions,
    labels = my_labels
  ) +
  theme_tree2() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

# Save plot --------------------------------------------------------------------

ggsave(
  "figures/Cyperaceae_tree_plot.pdf",
  Cyperaceae_tree_plot,
  width = 7, height = 17
)

ggsave(
  "figures/Cyperaceae_tree_plot.png",
  Cyperaceae_tree_plot,
  width = 7, height = 17, dpi = 300
)

