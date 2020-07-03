# Ruan van Mazijk, 2020

library(patchwork)  # Figure panelling

source("R/map-Schoenus-worldwide.R")
source("R/map-Schoenus-richness-South-Africa.R")
source("R/map-Schoenus-richness-Oceania.R")

maps <-
  worldwide_plot /
  proportion_sampled_plot /
  (south_africa_plot | oceania_plot) +
  plot_layout(heights = c(1, 1, 0.75))

ggsave(
  "figures/01-maps.pdf",
  maps,
  width = 7, height = 7
)

ggsave(
  "figures/01-maps.png",
  maps,
  width = 7, height = 7, dpi = 300
)
