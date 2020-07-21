# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(tidyverse)
library(phytools)   # Importing trees
library(ggtree)     # Multi-phylo plots
                    # (Intalled with BiocManager::install("ggtree"))
library(jntools)    # For ::get_tips_in_ape_plot_order()
                    # (Installed with remotes::install_github("joelnitta/jntools"))
library(patchwork)  # Figure panelling

# Import data ------------------------------------------------------------------

# RAxML-HPC reconstruction:
# Best tree with nodes' bootstrap support values
tree <- read.tree("data/phylogenies/2020-07-14_RAxML-HPC-reconstruction_04/RAxML_bipartitions.result")

species_DEC_areas <- read_csv("data/occurence-data/Schoenus-DEC-9areas.csv")

# Tidy data --------------------------------------------------------------------

# RAxML-HPC reconstruction:

# Best tree:
# Extract Schoenus
Schoenus <- tree %>%
  drop.tip(.$tip.label[!str_detect(.$tip.label, "Schoenus")]) %>%
  ladderize()

colnames(species_DEC_areas)[[1]] <- "id"
species_DEC_areas <- map_df(species_DEC_areas, as.character)
for (i in 1:nrow(species_DEC_areas)) {
  for (j in 2:ncol(species_DEC_areas)) {
    species_DEC_areas[i, j] <-
      if (species_DEC_areas[i, j] == "1") {
        colnames(species_DEC_areas)[[j]]
      } else {
        ""
      }
  }
}
species_DEC_areas$areas <- NA
species_DEC_areas$n_areas <- NA
for (i in 1:nrow(species_DEC_areas)) {
  areas <- sort(species_DEC_areas[i, 2:ncol(species_DEC_areas)])
  areas <- areas[areas != ""]
  species_DEC_areas$n_areas[[i]] <- length(areas)
  species_DEC_areas$areas[[i]] <- paste(areas, collapse = ", ")
}
species_DEC_areas_summary <- species_DEC_areas %>%
  select(id, areas, n_areas) %>%
  arrange(desc(n_areas)) %>%
  mutate(areas = case_when(
    n_areas >= 4                                               ~ "Cosmopolitan",
    str_detect(areas, "Western Australia")                     ~ "Western Australia",
    str_detect(areas, "(Australia|New Zealand|Tropical Asia)") ~ "Oceania",
    str_detect(areas, "Cape")                                  ~ "Cape",
    TRUE                                                       ~ areas
  )) %>%
  mutate(areas = factor(areas, unique(areas)))

# Plots ------------------------------------------------------------------------

Schoenus_BS_plot <- ggtree(Schoenus, ladderize = TRUE, right = TRUE)

facet_plot(Schoenus_BS_plot,
  geom = "geom_tile",
  data = species_DEC_areas_summary,
  panel = "DEC areas",
  aes(x = 0, fill = areas)
)
