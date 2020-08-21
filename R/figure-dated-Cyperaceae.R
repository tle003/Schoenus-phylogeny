# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(tidyverse)
library(phytools)
library(ggtree)   # Multi-phylo plots
                  # (Installed with BiocManager::install("ggtree"))
library(treeio)   # For ::read.beast()
library(jntools)  # For ::get_tips_in_ape_plot_order()
library(lemon)    # For ::coord_capped_cart()

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

Mapanioid_genera <- c(
  "Capitularina",
  "Chorizandra",
  "Chrysitrix",
  "Diplasia",
  "Exocarya",
  "Hypolytrum",
  "Lepironia",
  "Mapania",
  "Paramapania",
  "Principina",
  "Scirpodendron"
)

Mapanioideae_MRCA_node <- MCC_tree@phylo %>%
  getMRCA(.$tip.label[str_detect(.$tip.label,
    paste0("(", paste(Mapanioid_genera, collapse = "|"), ")")
  )])

# Get node numbers for clades to highlight
Clade_A_node    <- getMRCA(MCC_tree@phylo, c("Schoenus_insolitus", "Schoenus_sculptus"))
Clade_B_node    <- getMRCA(MCC_tree@phylo, c("Schoenus_falcatus",  "Schoenus_australis"))
Cape_clade_node <- getMRCA(MCC_tree@phylo, c("Schoenus_dregeanus", "Schoenus_australis"))

MCC_tree@phylo$tip.label <- str_replace(MCC_tree@phylo$tip.label, "_", " ")

# Plot -------------------------------------------------------------------------

# X-axis scaling things:
tree_height <- max(nodeHeights(MCC_tree@phylo))
my_labels <- c(90, 80, 70, 60, 50, 40, 30, 20, 10, 0)
label_positions <- tree_height - my_labels

# Try and fix node HPDs plotting backwards
nodes_with_future_HPDs <- MCC_tree@data %>%
  select(node, height_0.95_HPD) %>%
  unnest() %>%
  filter(height_0.95_HPD > tree_height) %>%
  pull(node)
MCC_tree@data %>%
  filter(node %in% nodes_with_future_HPDs) %>%
  as.data.frame()

# This doesn't work:
#MCC_tree@data <- MCC_tree@data %>%
#  mutate(height_0.95_HPD2 = purrr::map(height_0.95_HPD, rev))

# Neither does this:
#MCC_tree@data$height_0.95_HPD <- MCC_tree@data$height_0.95_HPD %>%
#  purrr::map(rev)
#  purrr::map(~{. - tree_height})

# Temporary solution used below:
# plot tree right-to-left instead of left-to-right
# (Don't know why, but this seems to work best)

# Alternatives:
# - FigTree
# - phytools:: + base::

# Make data for grey and white blocks for timescale-background of tree
my_panel_grid <- MCC_tree@phylo %>%
  get_tips_in_ape_plot_order() %>%
  map_dfr(~ tibble(
    x       = label_positions - 5,
    species = .x,
    alpha   = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
  )) %>%
  mutate(species = species %>%
    factor(levels = get_tips_in_ape_plot_order(MCC_tree@phylo)) %>%
    as.numeric()
  )

Cyperaceae_tree_plot <-
  ggtree(MCC_tree, ladderize = FALSE) +  # (already ladderized above!)
  geom_rootedge(rootedge = -10) +
  geom_tile(
    data = my_panel_grid,
    aes(x, species, alpha = alpha),
    fill = "black"
  ) +
  scale_alpha_manual(values = c(0, 0.2), guide = FALSE) +
  geom_hilight(node = Clade_A_node,    fill = "black", alpha = 0.125) +
  geom_hilight(node = Clade_B_node,    fill = "black", alpha = 0.250) +
  geom_hilight(node = Cape_clade_node, fill = "blue",  alpha = 0.125) +
  #annotate(geom = "text",
  #  label = "Clade A",
  #  x = 11, y = Ntip(Schoenus_MCC@phylo)   - 1
  #) +
  #annotate(geom = "text",
  #  label = "Clade B",
  #  x = 11, y = Ntip(Schoenus_MCC@phylo)/2 - 1
  #) +
  #annotate(geom = "text",
  #  label = "Cape clade",
  #  x = 31, y = Ntip(Schoenus_MCC@phylo)/2 - 1
  #) +
  geom_tiplab(
    aes(label = paste0('italic(\"', label, '\")')),
    parse = TRUE,
    size = 2.5,
    offset = -40
  ) +
  geom_hilight(node = Mapanioideae_MRCA_node) +
  geom_range("height_0.95_HPD",
    center = "height_median",
    size   = 1.5,
    alpha  = 0.4,
    colour = "darkblue"
  ) +
  #geom_hilight(Schoenus_MRCA_node,  fill = "darkblue",  alpha = 0.25) +
  #geom_hilight(Schoeneae_MRCA_node, fill = "lightblue", alpha = 0.25) +
  theme_tree2() +
  scale_x_reverse(name = "Ma",
    limits   = c(135, -10),
    breaks   = label_positions,
    labels   = my_labels
  ) +
  # Remove empty space above, below tree
  scale_y_continuous(limits = c(0, Ntip(MCC_tree@phylo) + 1), expand = c(0, 0)) +
  # Remove extra line at left of time axes
  coord_capped_cart(bottom = "left", top = "left") +
  # Move time axes' titles to the left
  theme(axis.title.x = element_text(hjust = 0.65))

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
