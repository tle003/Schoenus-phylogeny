# Ruan van Mazijk, 2020

library(patchwork)  # Figure panelling

source("R/map-Schoenus-worldwide.R")
source("R/map-Schoenus-richness-South-Africa.R")
source("R/map-Schoenus-richness-Oceania.R")

richness_maps <-
  worldwide_plot +
  (south_africa_plot | oceania_plot) +
  plot_layout(nrow = 2, heights = c(1, 0.75))

ggsave(
  "figures/01-maps.pdf",
  richness_maps,
  width = 10, height = 7
)

ggsave(
  "figures/01-maps.png",
  richness_maps,
  width = 10, height = 7, dpi = 300
)

sampling_maps <- n_spp_map / p_spp_map

ggsave(
  "figures/01-sampling.pdf",
  sampling_maps,
  width = 10, height = 7
)

ggsave(
  "figures/01-sampling.png",
  sampling_maps,
  width = 10, height = 7, dpi = 300
)
