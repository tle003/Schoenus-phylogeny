# Arrange Schoenus (Cyperaceae, Tribe Schoeneae) species richness and
# sampling effort maps (worlwide, in South Africa, Australia and New Zealand)

# Ruan van Mazijk, 2021

# Run other scripts  -----------------------------------------------------------

source("R/map-Schoenus-worldwide.R")
source("R/map-Schoenus-richness-South-Africa.R")
source("R/map-Schoenus-richness-Oceania.R")

# Load extra package -----------------------------------------------------------

library(patchwork)  # Figure panelling

# Panel the figures using patchwork syntax(!) ----------------------------------

maps <-
  richness_map /
  proportion_sampled_map /
  (south_africa_plot | oceania_plot) +
  plot_layout(heights = c(1, 1, 0.75))

# Save plot --------------------------------------------------------------------

# As PDF
ggsave(
  "figures/Schoenus-phy-bio-Fig3-RvM.pdf",
  maps,
  width = 7, height = 7
)
# As PNG
ggsave(
  "figures/Schoenus-phy-bio-Fig3-RvM.png",
  maps,
  width = 7, height = 7, dpi = 300
)
