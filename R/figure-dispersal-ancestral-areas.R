# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(tidyverse)
library(phytools)   # Importing trees
library(ggtree)     # Multi-phylo plots
                    # (Intalled with BiocManager::install("ggtree"))
library(treeio)     # For ::read.beast()
library(jntools)    # For ::get_tips_in_ape_plot_order()
                    # (Installed with remotes::install_github("joelnitta/jntools"))
library(lemon)      # For ::coord_capped_cart()

# Import data ------------------------------------------------------------------

# BEAST reconstruction (MCC tree), pruned to Schoeneae,
# with ancestral area probabilities for nodes (from DEC analysis):
Schoeneae_tree <- read.beast("BioGeoBEARS/tr_w_ancestral_areas2.tre")

# Biogeographical coding for extant species (used in DEC analysis):
Schoeneae_DEC_areas <- read_delim(
  "data/occurence-data/Schoeneae-DEC-9areas.txt",
  delim = " "
)

# Tidy data --------------------------------------------------------------------

# Extend tips that didn't reach time 0
Schoeneae_tree@phylo <-
  force.ultrametric(Schoeneae_tree@phylo, method = "extend")

# TODO(?):
#Schoeneae_tree@data <- Schoeneae_tree@data %>%
#  mutate(state = ifelse(relative_prob >= 0.5, state, NA))

# X-axis scaling things:
# (for both tree's time-axis labels and
# region-tiles' x-axis being on that same scale)
tree_height <- max(nodeHeights(Schoeneae_tree@phylo))
my_labels <- c(70, 60, 50, 40, 30, 20, 10, 0)
label_positions <- tree_height - my_labels

# Turn single CFWAZNPTH-column into multiple columns (complicated):
colnames(Schoeneae_DEC_areas)[[1]] <- "species"
write_file(
  paste0(colnames(Schoeneae_DEC_areas)[[2]], "\n"),
  "data/occurence-data/Schoeneae-DEC-9areas.tmp"
)
Schoeneae_DEC_areas %>%
  select(-species) %>%
  as.matrix() %>%
  paste(collapse = "\n") %>%
  write_file(
    "data/occurence-data/Schoeneae-DEC-9areas.tmp",
    append = TRUE
  )
Schoeneae_DEC_areas_only <- read_fwf(
  "data/occurence-data/Schoeneae-DEC-9areas.tmp",
  col_positions = fwf_widths(rep(1, times = 9))
)
colnames(Schoeneae_DEC_areas_only) <- Schoeneae_DEC_areas_only[1, ]
Schoeneae_DEC_areas_only <- Schoeneae_DEC_areas_only[-1, ]
Schoeneae_DEC_areas_only <- purrr::map_df(Schoeneae_DEC_areas_only, as.numeric)
Schoeneae_DEC_areas_tidy <- Schoeneae_DEC_areas_only %>%
  cbind(species = Schoeneae_DEC_areas$species) %>%
  as_tibble()
# TODO: export Schoeneae_DEC_areas_tidy as CSV

# To make the tiles etc. in the region panel easier to position
# (on the same scale as the tree) and then make narrower
my_scale_factor <-  5

# Tidy region data nicely and include x-axis position's column
Schoeneae_DEC_areas_tidy <- Schoeneae_DEC_areas_tidy %>%
  gather(area, present, -species) %>%
  filter(species != "Schoenus_adnatus") %>%
  mutate(
    species = factor(species, levels =
      get_tips_in_ape_plot_order(Schoeneae_tree@phylo)
    ),
    x = case_when(
      # TODO: refactor
      area == "C" ~ (label_positions[[1]] / my_scale_factor) - 12.5,
      area == "F" ~ (label_positions[[2]] / my_scale_factor) - 12.5,
      area == "W" ~ (label_positions[[3]] / my_scale_factor) - 12.5,
      area == "A" ~ (label_positions[[4]] / my_scale_factor) - 12.5,
      area == "Z" ~ (label_positions[[5]] / my_scale_factor) - 12.5,
      area == "N" ~ (label_positions[[6]] / my_scale_factor) - 12.5,
      area == "P" ~ (label_positions[[7]] / my_scale_factor) - 12.5,
      area == "T" ~ (label_positions[[8]] / my_scale_factor) - 12.5,
      area == "H" ~ ((tree_height + 10)   / my_scale_factor) - 12.5
    ),
    area =
      # TODO: functionalise
      case_when(
        area == "C" ~ "Cape",
        area == "F" ~ "Africa",
        area == "W" ~ "Western Australia",
        area == "A" ~ "Australia",
        area == "Z" ~ "New Zealand",
        area == "N" ~ "Neotropics",
        area == "P" ~ "Pacific",
        area == "T" ~ "Tropical Asia",
        area == "H" ~ "Holarctic"
      ) %>%
      factor(levels = c(
        "Cape",
        "Africa",
        "Western Australia",
        "Australia",
        "New Zealand",
        "Neotropics",
        "Pacific",
        "Tropical Asia",
        "Holarctic"
      )),
    present = as.logical(present)
  )

# Plots ------------------------------------------------------------------------

# Make data for grey and white blocks for timescale-background of tree
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

Schoeneae_tree_plot <-
  ggtree(Schoeneae_tree, ladderize = TRUE) +
  geom_rootedge(rootedge = 5) +
  geom_tile(
    data = my_panel_grid,
    aes(x, species, alpha = alpha),
    fill = "black"
  ) +
  scale_alpha_manual(values = c(0, 0.1), guide = FALSE) +
  geom_tiplab(
    aes(label = label %>%
      str_replace("Schoenus_", "S._") %>%
      str_replace("Gymnoschoenus_", "G._") %>%  # (otherwise very long name!)
      str_replace("_", " ") %>%
      {paste0('italic(\"', ., '\")')}
    ),
    parse  = TRUE,
    size   = 2.5,
    offset = 3
  ) +
  # Add time axis
  theme_tree2() +
  scale_x_continuous(
    limits   = c(-15, tree_height + 28),  # 28 is the min needed to fit labels
    breaks   = label_positions,
    labels   = -my_labels
  ) +
  # Remove empty space above, below tree
  scale_y_continuous(
    limits = c(0, Ntip(Schoeneae_tree@phylo) + 1),
    expand = c(0, 0)
  ) +
  # Remove extra line at right of time axes
  coord_capped_cart(bottom = "right", top = "right") +
  theme(axis.title.x = element_blank())

my_palette <- scales::brewer_pal(palette = "Paired")(
  n = length(unique(
    Schoeneae_DEC_areas_tidy$area
  ))
)
# Darken purple
my_palette[[9]] <- "#ab71c7"
my_palette2 <- vector("character", length = 2*length(my_palette))
for (i in 1:length(my_palette)) {
  pos2 <- 2*i
  pos1 <- pos2 - 1
  my_palette2[[pos1]] <- my_palette[[i]]
  my_palette2[[pos2]] <- "white"
}

Schoeneae_DEC_areas_plot <-
  facet_plot(Schoeneae_tree_plot,
    geom = "geom_tile",
    data = Schoeneae_DEC_areas_tidy,
    panel = "DEC areas",
    aes(
      x = x,
      colour = area,
      fill = factor(paste(area, !present), levels = c(
        "Cape FALSE",              "Cape TRUE",
        "Africa FALSE",            "Africa TRUE",
        "Western Australia FALSE", "Western Australia TRUE",
        "Australia FALSE",         "Australia TRUE",
        "New Zealand FALSE",       "New Zealand TRUE",
        "Neotropics FALSE",        "Neotropics TRUE",
        "Pacific FALSE",           "Pacific TRUE",
        "Tropical Asia FALSE",     "Tropical Asia TRUE",
        "Holarctic FALSE",         "Holarctic TRUE"
      ))
    ),
    width = 10 / my_scale_factor
  ) +
  scale_fill_manual(values = my_palette2, guide = FALSE) +
  scale_colour_manual(values = rep(NA, times = 9)) +
  guides(colour = guide_legend(
    title = "Region",
    override.aes = list(fill = my_palette)
  )) +
  theme(strip.text = element_blank())

# Manually remove region panel's "time"-axis
Schoeneae_DEC_areas_plot <- gridExtra::arrangeGrob(Schoeneae_DEC_areas_plot)
#plot(Schoeneae_DEC_areas_plot)
#str(Schoeneae_DEC_areas_plot, max.level = 1)
#Schoeneae_DEC_areas_plot$grobs[[1]]
Schoeneae_DEC_areas_plot$grobs[[1]]$grobs[[7]] <- zeroGrob()

# Manually remove grey and white blocks from region panel
Schoeneae_DEC_areas_plot$grobs[[1]]$grobs[[3]]$children[[6]] <-  zeroGrob()

# Check:
plot(Schoeneae_DEC_areas_plot)

# Save plot --------------------------------------------------------------------

ggsave(
  "figures/Schoeneae_DEC_areas_plot.pdf",
  Schoeneae_DEC_areas_plot,
  width = 10, height = 13
)
ggsave(
  "figures/Schoeneae_DEC_areas_plot.png",
  Schoeneae_DEC_areas_plot,
  width = 10, height = 13, dpi= 300
)
