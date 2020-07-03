# Map Schoenus (Cyperaceae, Tribe Schoeneae) species richness in
# quarter degree grid cells in South Africa

# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(rgdal)      # Shapefile I/O
library(raster)     # Country borders, other GIS functions
library(tidyverse)  # Data wrangling, figures

# Import data ------------------------------------------------------------------

# South African Schoenus richness in QDGC (from Tammy Elliott)
richness <- read_csv("data/Schoenus-South-Africa-richness.csv")

# South African border
border_ZAF <- getData("GADM", country = "ZAF", level = 0)

# Plot map ---------------------------------------------------------------------

# Set ggplot2 theme
theme_set(theme_bw() + theme(panel.grid = element_blank()))

south_africa_plot <- ggplot(richness) +
  aes(QDGC_lon, QDGC_lat, fill = richness) +
  geom_tile() +
  geom_polygon(
    data = border_ZAF,
    aes(x = long, y = lat, group = group),
    fill = NA, colour = "grey30", size = 0.2
  ) +
  annotate("text", label = "(c)", x = 17, y = -23) +
  coord_equal() +
  scale_x_continuous(
    breaks = c(20, 25, 30),
    labels = scales::label_math(expr = .x*"º")
  ) +
  scale_y_continuous(
    breaks = c(-25, -30, -35),
    labels = scales::label_math(expr = .x*"º")
  ) +
  scale_fill_viridis_c(
    direction = -1,
    breaks = c(5, 15, 25, 34),
    limits = c(0, 34)
  ) +
  theme(
    axis.title.y    = element_blank(),
    axis.title.x    = element_blank(),
    legend.position = "none"
  )
