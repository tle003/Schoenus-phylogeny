library(patchwork)

source("map-Schoenus-richness-Oceania.R")
source("map-Schoenus-richness-South-Africa.R")

both_plots <- south_africa_plot + oceania_plot

ggsave(
  "figure-1-maps.pdf",
  both_plots,
  width = 10, height = 4
)

ggsave(
  "figure-1-maps.png",
  both_plots,
  width = 10, height = 4, dpi = 300
)
