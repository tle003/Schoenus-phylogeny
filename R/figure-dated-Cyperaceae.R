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

subtribes <- read_csv("data/Schoeneae-subtribes.csv")

# Define helper functions ------------------------------------------------------

find_node <- function(tree, tip_pattern,
                            additional_taxa = NULL,
                            taxa_to_exclude = NULL) {
  pattern_matches <- tree$tip.label[str_detect(tree$tip.label, tip_pattern)]
  taxa <- c(pattern_matches, additional_taxa)
  taxa <- taxa[!(taxa %in% taxa_to_exclude)]
  getMRCA(tree, taxa)
}

vector2regex <- function(...) {
  x <- c(...)
  paste0("(",
    paste(x, collapse = "|"),
  ").+")
}

find_subtribe <- function(tree,
                          subtribe_name, subtribes_df,
                          additional_taxa_to_exclude = NULL) {
  genera <- subtribes_df %>%
    filter(subtribe == subtribe_name, is.na(species)) %>%
    pull(taxa)
  additional_taxa <- subtribes_df %>%
    filter(subtribe == subtribe_name, !is.na(species)) %>%
    pull(taxa)
  taxa_to_exclude <- subtribes_df %>%
    filter(subtribe != subtribe_name, genus %in% genera) %>%
    pull(taxa)
  find_node(tree,
    tip_pattern = vector2regex(genera),
    additional_taxa,
    c(taxa_to_exclude, additional_taxa_to_exclude)
  )
}

# Tidy data --------------------------------------------------------------------

MCC_tree@phylo <- force.ultrametric(MCC_tree@phylo, method = "extend")
MCC_tree@phylo$tip.label <- str_replace(MCC_tree@phylo$tip.label, "_", " ")

# Prune tree to Schoeneae-only
Schoeneae_node <-
  find_node(MCC_tree@phylo, "Schoenus", "Gymnoschoenus sphaerocephalus")
Schoenoid_taxa <- treeio::offspring(MCC_tree, Schoeneae_node)
non_Schoenoid_taxa <- MCC_tree@data$node %>%
  {.[!(. %in% Schoenoid_taxa)]} %>%
  as.numeric()
Schoeneae_tree <- drop.tip(MCC_tree, non_Schoenoid_taxa)

# Adjust height-related node-data by the tree-height
# to get HPDs to plot correctly (not backwards)
tree_height <- max(nodeHeights(Schoeneae_tree@phylo))
Schoeneae_tree@data <- Schoeneae_tree@data %>%
  mutate(
    height          = -height,
    height_median   = -height_median,
    height_0.95_HPD = purrr::map(height_0.95_HPD, ~-.)
  )
#Schoeneae_tree@data <- Schoeneae_tree@data %>%

tree_height <- -tree_height

plotTree(ladderize(Schoeneae_tree@phylo), node.numbers = TRUE, fsize = 0.5)
Schoeneae_tree %>%
  as_tibble() %>%
  mutate(
    height_min = map_dbl(height_0.95_HPD, min),
    height_max = map_dbl(height_0.95_HPD, max)
  ) %>%
  filter(is.na(label)) %>%
  select(node, height, height_min, height_median, height_max) %>%
  filter(height < 0)

subtribes <- subtribes %>%
  mutate(taxa = str_remove(paste(genus, species), " NA"))

# Identify nodes for clades to highlight ---------------------------------------

# .... Other Schoeneae subtribes -----------------------------------------------

subtribe_names <- unique(na.exclude(subtribes$subtribe))
subtribe_names <- subtribe_names[subtribe_names != "Schoeninae"]
clades_to_collapse <- purrr::map(subtribe_names, ~find_subtribe(
  Schoeneae_tree@phylo,
  subtribe_name = .x,
  subtribes_df = subtribes,
  additional_taxa_to_exclude = ifelse(.x == "Tricostulariinae",
    "Anthelepis paludosa",
    NA
  )
))
names(clades_to_collapse) <- subtribe_names

# .... Ingroup clades ----------------------------------------------------------

Schoenus_node <-
  find_node(Schoeneae_tree@phylo, "Schoenus")
Clade_A_node <-
  find_node(Schoeneae_tree@phylo, "Schoenus insolitus", "Schoenus sculptus")
Clade_B_node <-
  find_node(Schoeneae_tree@phylo, "Schoenus falcatus",  "Schoenus australis")
Cape_clade_node <-
  find_node(Schoeneae_tree@phylo, "Schoenus dregeanus", "Schoenus australis")

# Plot -------------------------------------------------------------------------

# .... X-axis scaling things ---------------------------------------------------

my_labels <- c(70, 60, 50, 40, 30, 20, 10, 0)
label_positions <- tree_height - my_labels

# .... Make data for grey and white blocks for timescale-background of tree ----

my_panel_grid <- Schoeneae_tree@phylo %>%
  get_tips_in_ape_plot_order() %>%
  map_dfr(~ tibble(
    x       = label_positions - 5,
    species = .x,
    alpha   = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE)
  )) %>%
  mutate(species = species %>%
    factor(levels = get_tips_in_ape_plot_order(Schoeneae_tree@phylo)) %>%
    as.numeric()
  )

# .... Main plot ---------------------------------------------------------------

Cyperaceae_tree_plot <-
  ggtree(Schoeneae_tree, ladderize = TRUE) +
  geom_rootedge(rootedge = 5) +
  geom_tile(
    data = my_panel_grid,
    aes(x, species, alpha = alpha),
    fill = "black"
  ) +
  scale_alpha_manual(values = c(0, 0.1), guide = FALSE) +
  geom_tiplab(
    aes(label = paste0('italic(\"', label, '\")')),
    parse = TRUE,
    size = 2.5,
    offset = 2
  ) +
  geom_range(
    "height_0.95_HPD", center = "height_median",
    colour = "darkblue", alpha = 0.5,
    size = 1.5
  ) +
  theme_tree2() +
  scale_x_continuous(name = "Ma",
    limits   = c(-50, 105),  # very wide to make space for annotations (below)
    breaks   = label_positions,
    labels   = my_labels
  ) +
  # Remove empty space above, below tree
  scale_y_continuous(
    limits = c(0, Ntip(Schoeneae_tree@phylo) + 1),
    expand = c(0, 0)
  ) +
  # Remove extra line at left of time axes
  coord_capped_cart(bottom = "right") +
  # Move time axes' titles to the left
  theme(axis.title.x = element_text(hjust = 0.35))

# .... Annotations -------------------------------------------------------------

clade_label_offset <- 25
clade_bar_extension <- 0.2

# Label other Schoeneae subtribes
for (subtribe in subtribe_names) {
  Cyperaceae_tree_plot <- Cyperaceae_tree_plot +
    geom_cladelabel(clades_to_collapse[[subtribe]],
      label =
        if (subtribe == "Oreobolus") {
          paste0('italic("Oreobolus")~clade')
        } else {
          subtribe
        },
      parse  = TRUE,
      offset = clade_label_offset,
      extend = clade_bar_extension
    )
}

# Label ingroup clades
Cyperaceae_tree_plot <- Cyperaceae_tree_plot +
  geom_cladelabel(Schoenus_node, paste0('italic("Schoenus")'), parse  = TRUE,
    offset = clade_label_offset,
    extend = clade_bar_extension
  ) +
  stat_hilight(
    node = Clade_A_node,
    colour = "darkblue", xmin = 20,
    fill = NA, alpha = NA
  ) +
  stat_hilight(
    node = Clade_B_node,
    colour = "darkblue", xmin = 20,
    fill = NA, alpha = NA
  ) +
  stat_hilight(
    node = Cape_clade_node,
    colour = "darkblue",
    fill = NA, alpha = NA
  ) +
  annotate(
    label = "Clade A", geom = "text",
    x = 25, y = Ntip(Schoeneae_tree@phylo) - 2
  ) +
  annotate(
    label = "Clade B", geom = "text",
    x = 25, y = (0.6 * Ntip(Schoeneae_tree@phylo)) - 2
  ) +
  annotate(
    label = "Cape clade", geom = "text",
    x = 42, y = (0.6 * Ntip(Schoeneae_tree@phylo)) - 2
  )

# Save plot --------------------------------------------------------------------

ggsave(
  "figures/Cyperaceae_tree_plot.pdf",
  Cyperaceae_tree_plot,
  width = 10, height = 15
)

ggsave(
  "figures/Cyperaceae_tree_plot.png",
  Cyperaceae_tree_plot,
  width = 10, height = 15, dpi = 300
)
