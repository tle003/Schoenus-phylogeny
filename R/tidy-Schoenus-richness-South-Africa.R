# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(rgdal)      # Shapefile I/O
library(raster)     # Country borders, other GIS functions
library(tidyverse)  # Data wrangling

# Import data ------------------------------------------------------------------

# South African Schoenus richness in QDGC (from Tammy Elliott)
richness <- read_csv("data/Schoenus-South-Africa-richness.csv")

# Larsen et al. (2009) QDGC for South Africa
QDGC <- readOGR("data/shapefiles/qdgc_zaf", layer = "qdgc_02_zaf")

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

# Export data ------------------------------------------------------------------

write_csv(richness, "data/Schoenus-South-Africa-richness.csv")
