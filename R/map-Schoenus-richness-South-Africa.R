# Map Schoenus (Cyperaceae, Tribe Schoeneae) species richness in
# quarter degree grid cells in South Africa

# Ruan van Mazijk, 2021

# Load packages ----------------------------------------------------------------

library(rgdal)      # Shapefile I/O
library(raster)     # Country borders, other GIS functions
library(tidyverse)  # Data wrangling, figures

# Import data ------------------------------------------------------------------

# South African Schoenus richness in QDGC (from Tammy Elliott)
richness <- read_csv("data/Schoenus-South-Africa-richness.csv")

# South African border
border_ZAF <- getData("GADM", country = "ZAF", level = 0)

# African continental border
border_world<- readOGR("data/shapefiles/World_Continents/v10/continent.gdb/")
# (Will zoom in on southern Africa in ggplot2:: code below)

# Plot map ---------------------------------------------------------------------

# Set ggplot2 theme
theme_set(theme_bw() + theme(panel.grid = element_blank()))

# Note, the panels of Figure 2 are as follows:
#
#     A--------------------\
#     |                    |
#     | richness           |
#     |                    |
#     \--------------------/
#     B--------------------\
#     |                    |
#     | proportion sampled |
#     |                    |
#     \--------------------/
#     C--------\ D---------\
#     |        | |         |
#     | S Afr  | | W Aus   |
#     |        | |         |
#     \--------/ \---------/
#
# Consequently, will omit the legend for C (shared with D)

# Panel C
south_africa_plot <- ggplot(richness) +
  # Plot coloured tiles representing QDGC-scale richness
  aes(QDGC_lon, QDGC_lat, fill = richness) +
  geom_tile() +
  # Plot African and South African borders
  geom_polygon(
    data = border_world,
    aes(x = long, y = lat, group = group),
    fill = NA, colour = "grey30", size = 0.15
  ) +
  geom_polygon(
    data = border_ZAF,
    aes(x = long, y = lat, group = group),
    fill = NA, colour = "grey30", size = 0.15
  ) +
  # Panel label
  annotate("text", label = "C", x = 17, y = -23) +
  # Colour scheme customisation
  scale_fill_viridis_c(
    direction = -1,
    breaks = c(5, 15, 25, 34),
    limits = c(0, 34)
  ) +
  # Lat/lon scale customisation
  coord_equal(xlim = c(16, 33), ylim = c(-36, -22)) +
  scale_x_continuous(
    breaks = c(20, 25, 30),
    labels = scales::label_math(expr = .x*"º")
  ) +
  scale_y_continuous(
    breaks = c(-25, -30, -35),
    labels = scales::label_math(expr = .x*"º")
  ) +
  theme(
    axis.title.x    = element_blank(),
    axis.title.y    = element_blank(),
    legend.position = "none"  # Omit the legend
  )
