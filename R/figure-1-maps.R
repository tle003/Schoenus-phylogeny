# Arrange Schoenus (Cyperaceae, Tribe Schoeneae) species richness and
# sampling effort maps (worlwide, in South Africa, Australia and New Zealand)

# Ruan van Mazijk, 2020

# Run other scripts  -----------------------------------------------------------

source("R/map-Schoenus-worldwide.R")
source("R/map-Schoenus-richness-South-Africa.R")
source("R/map-Schoenus-richness-Oceania.R")

# Load extra package -----------------------------------------------------------

library(patchwork)  # Figure panelling

# Panel the figures using patchwork syntax(!) ----------------------------------

maps <-
  worldwide_plot /
  proportion_sampled_plot /
  (south_africa_plot | oceania_plot) +
  plot_layout(heights = c(1, 1, 0.75))

# Save plot --------------------------------------------------------------------

ggsave(
  "figures/Schoenus-phy-bio-Fig1-RvM.pdf",
  maps,
  width = 7, height = 7
)
ggsave(
  "figures/Schoenus-phy-bio-Fig1-RvM.png",
  maps,
  width = 7, height = 7, dpi = 300
)
