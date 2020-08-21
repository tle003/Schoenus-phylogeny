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
  ladderize(right = FALSE)

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

# Clades to collapse
subtribes <- read_csv("data/Schoeneae-subtribes.csv")
as.data.frame(subtribes)

# TODO: funcionalise
Caustiinae_node <- MCC_tree@phylo %>%
  getMRCA(c(
    .$tip.label[str_detect(.$tip.label,
      paste0("(", paste(subtribes$Genus[1:2], collapse = "|"), ").+")
    )],
    paste(subtribes$Genus[[3]], subtribes$Species[[3]])
  ))
Gahniinae_node <- MCC_tree@phylo %>%
  getMRCA(c(
    .$tip.label[str_detect(.$tip.label,
      paste0("(", paste(subtribes$Genus[4:7], collapse = "|"), ").+")
    )]
  ))
Lepidosperminae_node <- MCC_tree@phylo %>%
  getMRCA(c(
    paste(subtribes$Genus[[8]], subtribes$Species[[8]]),
    .$tip.label[str_detect(.$tip.label,
      paste0("(", paste(subtribes$Genus[9:11], collapse = "|"), ").+")
    )]
  ))
Oreobolus_node <- MCC_tree@phylo %>%
  getMRCA(c(
    .$tip.label[str_detect(.$tip.label,
      paste0("(", paste(subtribes$Genus[12:16], collapse = "|"), ").+")
    )]
  ))
#Tricostulariinae_species <- MCC_tree@phylo %>%
#  {c(
#    paste(subtribes$Genus[[19]], subtribes$Species[[19]]),
#    .$tip.label[str_detect(.$tip.label,
#      paste0("(", paste(subtribes$Genus[c(18, 20:25)], collapse = "|"), ").+")
#    )] %>%
#    (function(x) {
#      x[x != paste(subtribes$Genus[[3]], subtribes$Species[[3]])]
#    })() %>%
#    (function(x) {
#      x[x != paste(subtribes$Genus[[8]], subtribes$Species[[8]])]
#    })()
#  )}
#foo <- MCC_tree@phylo
#foo$tip.label <- ifelse(foo$tip.label %in% Tricostulariinae_species, foo$tip.label, " ")
#plotTree(foo, fsize = 0.25)
Tricostulariinae_node <- MCC_tree@phylo %>%
  getMRCA(c(
    paste(subtribes$Genus[[19]], subtribes$Species[[19]]),
    .$tip.label[str_detect(.$tip.label,
      paste0("(", paste(subtribes$Genus[c(18, 20:25)], collapse = "|"), ").+")
    )] %>%
    (function(x) {
      x[x != paste(subtribes$Genus[[3]], subtribes$Species[[3]])]
    })() %>%
    (function(x) {
      x[x != paste(subtribes$Genus[[8]], subtribes$Species[[8]])]
    })() %>%
    (function(x) {
      x[x != "Anthelepis paludosa"]
    })()
  ))
Gymnoschoeninae_node <- MCC_tree@phylo %>%
  getMRCA(c(
    .$tip.label[str_detect(.$tip.label,
      paste0("(", paste(subtribes$Genus[26:27], collapse = "|"), ").+")
    )]
  ))

clades_to_collapse <- list(
  Caustiinae       = Caustiinae_node,
  Gahniinae        = Gahniinae_node,
  Lepidosperminae  = Lepidosperminae_node,
  Oreobolus        = Oreobolus_node,
  Tricostulariinae = Tricostulariinae_node,
  Gymnoschoeninae  = Gymnoschoeninae_node
)

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

clade_label_offset <- 85
clade_bar_extension <- 0.2

Cyperaceae_tree_plot <-
  ggtree(MCC_tree, ladderize = FALSE) +  # (already ladderized above!)
  geom_rootedge(rootedge = -10) +
  geom_cladelabel(node = clades_to_collapse$Caustiinae, label = "Caustiinae", offset = clade_label_offset, extend = clade_bar_extension) +
  geom_cladelabel(node = clades_to_collapse$Gahniinae, label = "Gahniinae", offset = clade_label_offset, extend = clade_bar_extension) +
  geom_cladelabel(node = clades_to_collapse$Lepidosperminae, label = "Lepidosperminae", offset = clade_label_offset, extend = clade_bar_extension) +
  geom_cladelabel(node = clades_to_collapse$Oreobolus, label = paste0('italic("Oreobolus")~clade'), offset = clade_label_offset, extend = clade_bar_extension, parse = TRUE) +
  geom_cladelabel(node = clades_to_collapse$Tricostulariinae, label = "Tricostulariinae", offset = clade_label_offset, extend = clade_bar_extension) +
  geom_cladelabel(node = clades_to_collapse$Gymnoschoeninae, label = "Gymnoschoeninae", offset = clade_label_offset, extend = clade_bar_extension) +
  geom_tile(
    data = my_panel_grid,
    aes(x, species, alpha = alpha),
    fill = "black"
  ) +
  scale_alpha_manual(values = c(0, 0.2), guide = FALSE) +
  geom_cladelabel(node = Clade_A_node, label = "Clade A", offset = clade_label_offset - 20, extend = clade_bar_extension) +
  geom_cladelabel(node = Clade_B_node, label = "Clade B", offset = clade_label_offset - 20, extend = clade_bar_extension) +
  geom_cladelabel(node = Cape_clade_node, label = "Cape clade", offset = clade_label_offset - 45, extend = clade_bar_extension) +
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
  geom_range("height_0.95_HPD",
    center = "height_median",
    size   = 1.5,
    alpha  = 0.4,
    colour = "darkblue"
  ) +
  geom_cladelabel(node = Schoenus_MRCA_node, label = paste0('italic("Schoenus")'), offset = clade_label_offset, extend = clade_bar_extension, parse = TRUE) +
  geom_cladelabel(node = Schoeneae_MRCA_node, label = "Schoeneae", offset = clade_label_offset + 45, extend = clade_bar_extension) +
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
