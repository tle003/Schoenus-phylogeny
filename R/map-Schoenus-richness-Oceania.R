# Map Schoenus (Cyperaceae, Tribe Schoeneae) species richness in
# quarter degree grid cells in Australia and New Zealand

# Ruan van Mazijk, 2021

# Load packages ----------------------------------------------------------------

library(rgdal)      # Shapefile I/O
library(raster)     # Country borders, other GIS functions
library(tidyverse)  # Data wrangling, figures

# Import data ------------------------------------------------------------------

# Schoenus richness in QDGC
richness <- read_csv("data/Schoenus-Oceania-richness.csv")

# Borders
border_AUS <- getData("GADM", country = "AUS", level = 0)
border_NZL <- getData("GADM", country = "NZL", level = 0)

# Tidy borders -----------------------------------------------------------------

border <- rbind(border_AUS, border_NZL)
# Have a look (WARNING: takes a while)
#plot(border)
# Lots of stray islands...

# Crop to focus on the mainlands
main_extent <- extent(110, 180, -50, -10)
border_cropped <- crop(border, main_extent)

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
# Consequently, will keep the legend for D (shared with and omitted in C)

# Panel D
oceania_plot <- ggplot(richness) +
  # Plot coloured tiles representing QDGC-scale richness
  aes(QDGC_lon, QDGC_lat, fill = richness) +
  geom_tile() +
  # Plot borders of Australia, New Zealand et al.
  geom_polygon(
    data = border_cropped,
    aes(x = long, y = lat, group = group),
    fill = NA, colour = "grey30", size = 0.15
  ) +
  # Panel label
  annotate("text", label = "D", x = 115, y = -13) +
  # Colour scheme customisation
  scale_fill_viridis_c(
    name = "No. species\nper QDGC",
    direction = -1,
    breaks = c(5, 15, 25, 34),
    limits = c(0, 34)
  ) +
  # Lat/lon scale customisation
  coord_equal() +
  scale_x_continuous(
    breaks = c(115, 125, 135, 145, 155, 165, 175),
    labels = scales::label_math(expr = .x*"º")
  ) +
  scale_y_continuous(
    breaks = c(-15, -25, -35, -45),
    labels = scales::label_math(expr = .x*"º")
  ) +
  theme(
    axis.title.x      = element_blank(),
    axis.title.y      = element_blank(),
    legend.text.align = 1
  )
