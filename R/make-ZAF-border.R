# Warning: this script will take a few hours to execute

# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(rgdal)   # Shapefile I/O
library(raster)  # Country borders, other GIS functions
library(rmapshaper)
# (See: <https://cran.r-project.org/web/packages/rmapshaper/vignettes/rmapshaper.html>)

# Import data ------------------------------------------------------------------

# South Africa border
border_ZAF <- getData("GADM", country = "ZAF", level = 0)

# Tidy -------------------------------------------------------------------------

# Simplify polygons (don't need that much detail)
border_ZAF <- ms_simplify(border_ZAF, keep = 0.001)

# Export final border ----------------------------------------------------------

writeOGR(
  border_ZAF, "border_ZAF",
  layer = names(border_ZAF)[[1]],
  driver = "ESRI Shapefile"
)
