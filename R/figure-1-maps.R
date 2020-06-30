library(patchwork)

source("R/map-Schoenus-richness-Oceania.R")
source("R/map-Schoenus-richness-South-Africa.R")

both_plots <- south_africa_plot + oceania_plot

ggsave(
  "figures/01-maps.pdf",
  both_plots,
  width = 10, height = 4
)

ggsave(
  "figures/01-maps.png",
  both_plots,
  width = 10, height = 4, dpi = 300
)
