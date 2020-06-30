# Warning: this script will take a few hours to execute

# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(rgdal)   # Shapefile I/O
library(raster)  # Country borders, other GIS functions
library(rmapshaper)
# (See: <https://cran.r-project.org/web/packages/rmapshaper/vignettes/rmapshaper.html>)

# Import data ------------------------------------------------------------------

# Larsen et al. (2009) QDGC for Australia
QDGC_AUS <- readOGR("qdgc_aus", layer = "qdgc_02_aus")

# Australia border
border_AUS <- getData("GADM", country = "AUS", level = 0)

# Tidy -------------------------------------------------------------------------

# Crop to focus on the mainland
border_AUS <- crop(border_AUS, QDGC_AUS)

# Simplify polygons (don't need that much detail)
border_AUS <- ms_simplify(border_AUS, keep = 0.001)

# Export final border ----------------------------------------------------------

writeOGR(
  border_AUS, "border_AUS",
  layer = names(border_AUS)[[1]],
  driver = "ESRI Shapefile"
)
