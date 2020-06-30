# Warning: this script will take a few hours to execute

# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(rgdal)   # Shapefile I/O
library(raster)  # Country borders, other GIS functions
library(rmapshaper)
# (See: <https://cran.r-project.org/web/packages/rmapshaper/vignettes/rmapshaper.html>)

# Import data ------------------------------------------------------------------

# New Caledonia border
border_NCL <- getData("GADM", country = "NCL", level = 0)

# Tidy -------------------------------------------------------------------------

# Simplify polygons (don't need that much detail)
border_NCL <- ms_simplify(border_NCL, keep = 0.001)

# Export final border ----------------------------------------------------------

writeOGR(
  border_NCL, "data/shapefiles/border_NCL",
  layer = names(border_NCL)[[1]],
  driver = "ESRI Shapefile"
)
