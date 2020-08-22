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

# Tidy data --------------------------------------------------------------------

MCC_tree@phylo <- MCC_tree@phylo %>%
  force.ultrametric(method = "extend") %>%
  ladderize(right = FALSE)

MCC_tree@phylo$tip.label <- str_replace(MCC_tree@phylo$tip.label, "_", " ")

subtribes <- subtribes %>%
  mutate(taxa = str_remove(paste(genus, species), " NA"))

# Identify nodes for clades to highlight ---------------------------------------

# .... Define helper functions -------------------------------------------------

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

# .... Outgroup taxa -----------------------------------------------------------

Mapanioid_genera <- vector2regex(  # Source: Wikispecies (Accessed 2020-08-21)
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
Mapanioideae_node <- find_node(MCC_tree@phylo, Mapanioid_genera)

# TODO: cont.

# .... Other Schoeneae subtribes -----------------------------------------------

subtribe_names <- unique(na.exclude(subtribes$subtribe))
subtribe_names <- subtribe_names[subtribe_names != "Schoeninae"]
clades_to_collapse <- purrr::map(subtribe_names, ~find_subtribe(
  MCC_tree@phylo,
  subtribe_name = .x,
  subtribes_df = subtribes,
  additional_taxa_to_exclude = ifelse(.x == "Tricostulariinae",
    "Anthelepis paludosa",
    NA
  )
))
names(clades_to_collapse) <- subtribe_names

# .... Ingroup clades ----------------------------------------------------------

Schoeneae_node <-
  find_node(MCC_tree@phylo, "Schoenus", "Gymnoschoenus sphaerocephalus")
Schoenus_node <-
  find_node(MCC_tree@phylo, "Schoenus")
Clade_A_node <-
  find_node(MCC_tree@phylo, "Schoenus insolitus", "Schoenus sculptus")
Clade_B_node <-
  find_node(MCC_tree@phylo, "Schoenus falcatus",  "Schoenus australis")
Cape_clade_node <-
  find_node(MCC_tree@phylo, "Schoenus dregeanus", "Schoenus australis")

# Plot -------------------------------------------------------------------------

# .... X-axis scaling things ---------------------------------------------------

tree_height <- max(nodeHeights(MCC_tree@phylo))
my_labels <- c(90, 80, 70, 60, 50, 40, 30, 20, 10, 0)
label_positions <- tree_height - my_labels

# .... Make data for grey and white blocks for timescale-background of tree ----
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

# .... Main plot ---------------------------------------------------------------

Cyperaceae_tree_plot <-
  ggtree(MCC_tree, ladderize = FALSE) +  # (already ladderized above!)
  geom_cladelabel(node = clades_to_collapse$Caustiinae, label = "Caustiinae", offset = clade_label_offset, extend = clade_bar_extension) +
  geom_cladelabel(node = clades_to_collapse$Gahniinae, label = "Gahniinae", offset = clade_label_offset, extend = clade_bar_extension) +
  geom_cladelabel(node = clades_to_collapse$Lepidosperminae, label = "Lepidosperminae", offset = clade_label_offset, extend = clade_bar_extension) +
  geom_cladelabel(node = clades_to_collapse$Oreobolus, label = paste0('italic("Oreobolus")~clade'), offset = clade_label_offset, extend = clade_bar_extension, parse = TRUE) +
  geom_cladelabel(node = clades_to_collapse$Tricostulariinae, label = "Tricostulariinae", offset = clade_label_offset, extend = clade_bar_extension) +
  geom_cladelabel(node = clades_to_collapse$Gymnoschoeninae, label = "Gymnoschoeninae", offset = clade_label_offset, extend = clade_bar_extension) +
  geom_rootedge(rootedge = 10) +
  geom_tile(
    data = my_panel_grid,
    aes(x, species, alpha = alpha),
    fill = "black"
  ) +
  scale_alpha_manual(values = c(0, 0.2), guide = FALSE) +
  geom_cladelabel(node = Clade_A_node, label = "Clade A", offset = clade_label_offset - 20, extend = clade_bar_extension) +
  geom_cladelabel(node = Clade_B_node, label = "Clade B", offset = clade_label_offset - 20, extend = clade_bar_extension) +
  geom_cladelabel(node = Cape_clade_node, label = "Cape clade", offset = clade_label_offset - 45, extend = clade_bar_extension) +
  geom_tiplab(
    aes(label = paste0('italic(\"', label, '\")')),
    parse = TRUE,
    size = 2.5,
    offset = 2
  ) +
  geom_cladelabel(node = Schoenus_node, label = paste0('italic("Schoenus")'), offset = clade_label_offset, extend = clade_bar_extension, parse = TRUE) +
  geom_cladelabel(node = Schoeneae_node, label = "Schoeneae", offset = clade_label_offset + 45, extend = clade_bar_extension) +
  theme_tree2() +
  scale_x_continuous(name = "Ma",
    limits   = c(-10, 235),  # very wide to make space for annotations (below)
    breaks   = label_positions,
    labels   = my_labels
  ) +
  # Remove empty space above, below tree
  scale_y_continuous(limits = c(0, Ntip(MCC_tree@phylo) + 1), expand = c(0, 0)) +
  # Remove extra line at left of time axes
  coord_capped_cart(bottom = "right") +
  # Move time axes' titles to the left
  theme(axis.title.x = element_text(hjust = 0.35))

# .... Annotations -------------------------------------------------------------

clade_label_offset <- 85
clade_bar_extension <- 0.2

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
