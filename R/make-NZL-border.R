# Warning: this script will take an hour or two to execute

# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(rgdal)   # Shapefile I/O
library(raster)  # Country borders, other GIS functions
library(rmapshaper)
# (See: <https://cran.r-project.org/web/packages/rmapshaper/vignettes/rmapshaper.html>)

# Import data ------------------------------------------------------------------

# Larsen et al. (2009) QDGC for New Zealand
QDGC_NZL <- readOGR("data/shapefiles/qdgc_oceania", layer = "qdgc_02_oceania")

# New Zealand border
border_NZL <- getData("GADM", country = "NZL", level = 0)

# Tidy -------------------------------------------------------------------------

# Crop to focus on the mainland
border_NZL <- crop(border_NZL, QDGC_NZL)

# Simplify polygons (don't need that much detail)
border_NZL <- ms_simplify(border_NZL, keep = 0.001)

# Export final border ----------------------------------------------------------

writeOGR(
  border_NZL, "data/shapefiles/border_NZL",
  layer = names(border_NZL)[[1]],
  driver = "ESRI Shapefile"
)
