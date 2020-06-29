# Map Schoenus (Cyperaceae, Tribe Schoeneae) species richness in
# quarter degree grid cells in Australia & New Zealand

# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(rgdal)      # Shapefile I/O
library(raster)     # Country borders, other GIS functions
library(tidyverse)  # Data wrangling

# Import data ------------------------------------------------------------------

# South African Schoenus richness in QDGC (from Tammy Elliott)
richness <- read_csv("Schoenus-occurence-data/SASchoenus_grid_cell_freq.csv")

# Larsen et al. (2009) QDGC for South Africa
QDGC <- readOGR("qdgc_zaf", layer = "qdgc_02_zaf")

# South African border
border_ZAF <- readOGR("border_ZAF", layer = "GID_0")

# Tidy data --------------------------------------------------------------------

# Neaten up data frame
richness <- richness[, -1]
colnames(richness) <- c("QDGC", "richness")

# Make data frame and QDGC shapefile codes comparable
QDGC$qdgc <- QDGC$qdgc %>%
  str_remove("E0") %>%
  str_remove("S")
# Remove spaces
richness$QDGC <- str_remove_all(richness$QDGC, " ")
# Swap lat/lon in South African richness QDGC codes
lat <- substr(richness$QDGC, 1, 2)
lon <- substr(richness$QDGC, 3, 4)
substr(richness$QDGC, 1, 2) <- lon
substr(richness$QDGC, 3, 4) <- lat

# Add lon/lat of QDGC-midpoints to data frame ----------------------------------

richness$QDGC_lon <- NA
richness$QDGC_lat <- NA
for (i in 1:nrow(richness)) {
  richness$QDGC_lon[[i]] <-
    as.numeric(as.character(
      QDGC$lon[
        QDGC$qdgc == richness$QDGC[[i]]
      ]
    ))
  richness$QDGC_lat[[i]] <-
    as.numeric(as.character(
      QDGC$lat[
        QDGC$qdgc == richness$QDGC[[i]]
      ]
    ))
}

# Note:
max(richness$richness)
# 17 is richest square

# Plot maps --------------------------------------------------------------------

# Set ggplot2 theme
theme_set(theme_bw())

ggplot(richness) +
  aes(QDGC_lon, QDGC_lat, fill = richness) +
  geom_tile() +
  geom_polygon(
    data = border_ZAF,
    aes(x = long, y = lat, group = group),
    colour = "black", fill = NA, size = 0.25
  ) +
  coord_equal() +
  scale_fill_viridis_c(
    name = bquote(italic("Schoenus")*" species richness"),
    direction = -1,
    breaks = c(5, 15, 25, 34), limits = c(0, 34)
  ) +
  labs(x = "Longitude (º)", y = "Latitude (º)") +
  theme(legend.position = "top")
