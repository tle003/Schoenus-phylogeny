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
species_DEC_areas <- species_DEC_areas %>%
  gather(area, present, -id) %>%
  mutate(present = as.logical(present))

# Plots ------------------------------------------------------------------------

Schoenus_BS_plot <- ggtree(Schoenus, ladderize = TRUE, right = TRUE) #+
  #geom_tiplab(
  #  aes(label = paste0('italic(\"', label, '\")')),
  #  parse = TRUE,
  #  size  = 2.5,
  #  align = TRUE
  #) +
  #xlim(0, 0.5)

facet_plot(Schoenus_BS_plot,
  geom = "geom_tile",
  data = species_DEC_areas,
  panel = "DEC areas",
  aes(x = as.numeric(as.factor(area)), fill = area, alpha = present)
) +
  scale_alpha_manual(values = c(0, 1))
