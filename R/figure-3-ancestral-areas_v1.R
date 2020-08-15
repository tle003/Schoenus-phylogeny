# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(tidyverse)
library(phytools)   # Importing trees
library(ggtree)     # Multi-phylo plots
                    # (Intalled with BiocManager::install("ggtree"))
library(treeio)     # For ::read.beast()
library(jntools)    # For ::get_tips_in_ape_plot_order()
                    # (Installed with remotes::install_github("joelnitta/jntools"))
library(patchwork)  # Figure panelling

# Import data ------------------------------------------------------------------

# RAxML-HPC reconstruction:
# Best tree with nodes' bootstrap support values
MCC_tree <- read.beast("data/phylogenies/Cyperaceae-all-taxa-6calib-max-clad-AUG12.tre")

# Biogeographical coding for extant species (used in DEC analysis)
Schoeneae_DEC_areas <- read_delim("data/occurence-data/Schoeneae-DEC-9areas.txt", delim = " ")

# Tidy data --------------------------------------------------------------------


MCC_tree@phylo <- MCC_tree@phylo %>%
  force.ultrametric(method = "extend") %>%
  ladderize(right = TRUE)
Schoeneae_MRCA_node <- MCC_tree@phylo %>%
  getMRCA(c(
    .$tip.label[str_detect(.$tip.label, "Schoenus")],
    "Gymnoschoenus_sphaerocephalus"
  ))
Schoeneae_tree <- extract.clade(MCC_tree@phylo, Schoeneae_MRCA_node)

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
Schoeneae_DEC_areas_tidy <-
  cbind(species = Schoeneae_DEC_areas$species, Schoeneae_DEC_areas_only) %>%
  as_tibble()

Schoeneae_DEC_areas_tidy %>%
  gather(area, present, -species) %>%
  mutate(
    species = factor(species, levels = get_tips_in_ape_plot_order(Schoenus)),
    area = factor(area, levels = c(
      "Cape", "Africa",
      "Western Australia", "Australia",
      "New Zealand", "Neotropics",
      "Pacific", "Tropical Asia", "Holarctic"
    )),
    area_group =
      case_when(
        area %in% c("Cape", "Western Australia") ~ "Cape + Western Australia",
        area %in% c("Africa", "Australia")       ~ "Africa + Australia",
        area %in% c("New Zealand", "Neotropics") ~ "New Zealand + Neotropics",
        TRUE                                     ~ as.character(area)
      ) %>%
      factor(levels = c(
        "Cape + Western Australia",
        "Africa + Australia",
        "New Zealand + Neotropics",
        "Pacific", "Tropical Asia", "Holarctic"
      )),
    present = as.logical(present),
  )

# Plots ------------------------------------------------------------------------

Schoeneae_tree_plot <-
  ggtree(Schoeneae) +
  geom_tiplab(
    aes(label = label %>%
      str_replace("Schoenus_", "S._") %>%
      str_replace("_", " ") %>%
      {paste0('italic(\"', ., '\")')}
    ),
    parse = TRUE,
    size = 2.5
  ) +
  theme_tree2() +
  scale_x_continuous(name = "Ma",
    limits = c(-15, tree_height + 20),
    breaks = label_positions,
    labels = my_labels
  )

Schoenus_DEC_areas_plot <-
  facet_plot(Schoenus_BS_plot,
    geom = "geom_tile",
    data = Schoenus_DEC_areas_tidy,
    panel = "DEC areas",
    aes(x = 0, group = area_group, fill = area, alpha = present),
    width = tile_width, position = position_dodge(width = tile_width)
  ) +
  scale_fill_brewer(palette = "Paired") +
  scale_alpha_manual(values = c(0, 1), guide = FALSE) +
  theme(strip.text = element_blank())

ggsave("figures/Schoenus_DEC_areas_plot.pdf", Schoenus_DEC_areas_plot, width = 10, height = 12)

####

#theme_set(theme_classic())
#Schoenus_areas_plot <- ggplot(Schoenus_DEC_areas_tidy) +
#  aes(area_group, species, fill = area, alpha = present) +
#  geom_tile() +
#  scale_alpha_manual(values = c(0, 1), guide = FALSE) +
#  theme(
#    axis.title   = element_blank(),
#    axis.text.x  = element_blank(), axis.text.y = element_blank(),
#    axis.ticks.x = element_blank(),
#    axis.line.x  = element_blank(),
#    plot.margin  = unit(c(2, 0, 2, 0), "cm")
#  )
#
#library(patchwork)
#{Schoenus_BS_plot + xlim(0, 0.5) + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))} + Schoenus_areas_plot
