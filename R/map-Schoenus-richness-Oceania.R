# Map Schoenus (Cyperaceae, Tribe Schoeneae) species richness in
# quarter degree grid cells in Australia and New Zealand

# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(rgdal)      # Shapefile I/O
library(raster)     # Country borders, other GIS functions
library(tidyverse)  # Data wrangling, figures

# Import data ------------------------------------------------------------------

# Schoenus richness in QDGC
richness <- read_csv("data/Schoenus-Oceania-richness.csv")

# Borders
border_AUS <- readOGR("data/shapefiles/border_AUS", layer = "GID_0")
border_NZL <- readOGR("data/shapefiles/border_NZL", layer = "GID_0")

# Tidy borders -----------------------------------------------------------------

border <- rbind(border_AUS, border_NZL)
# Have a look (WARNING: takes a while)
#plot(border)
# Lots of stray islands...

# Crop to focus on the mainlands (NOTE: takes a while)
main_extent <-extent(110, 180, -50, -10)
border_cropped <- crop(border, main_extent)

# Plot map ---------------------------------------------------------------------

# Set ggplot2 theme
theme_set(theme_bw() + theme(panel.grid = element_blank()))

oceania_plot <- ggplot(richness) +
  aes(QDGC_lon, QDGC_lat, fill = richness) +
  geom_tile() +
  geom_polygon(
    data = border_cropped,
    aes(x = long, y = lat, group = group),
    colour = "black", fill = NA, size = 0.25
  ) +
  coord_equal() +
  scale_x_continuous(breaks = c(115, 125, 135, 145, 155, 165, 175), labels = scales::label_math(expr = .x*"º")) +
  scale_y_continuous(breaks = c(-15, -25, -35, -45),                labels = scales::label_math(expr = .x*"º")) +
  scale_fill_viridis_c(
    name = "No. species\nper QDGC",
    direction = -1,
    breaks = c(5, 15, 25, 34),
    limits = c(0, 34)
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  )
