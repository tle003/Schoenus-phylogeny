# Tammy Elliott and Ruan van Mazijk, 2020


# Load the required libraries.
library(sp)
library(rJava)
require(raster)
require(dismo)
library(rgdal)
library(maptools)
library(rgeos)
library(tidyverse)
library(caret)
require(caret)
library(faraway)
library(olsrr)
library(SSDM)
library(gridExtra)
library(RColorBrewer)
library(ggplot2)
library(sf)
library(dplyr)
library(mapproj)
library(rgeos)

# Import data ------------------------------------------------------------------

Schoenus_worldwide <- read.csv("data/occurence-data/Schoenus-worldwide-TDWG.csv")

TDWG <- readOGR("~/wgsrpd/level3", layer = "level3", stringsAsFactors = FALSE)

# Tidy data --------------------------------------------------------------------

# Give species names as rownames
rownames(Schoenus_worldwide) <- Schoenus_worldwide$Species

# Calculate rowsums to get a vector with values per region
Schoenus_worldwide_sums <- as.data.frame(colSums(
  Schoenus_worldwide[, 5:ncol(Schoenus_worldwide)]
))
Schoenus_worldwide_sums <- data.frame(
  LEVEL3_COD = rownames(Schoenus_worldwide_sums),
  Count      = Schoenus_worldwide_sums
)

# Create data frame for TDWG level 3 data
TDWG_level3_df <- data.frame(
  LEVEL3_COD = TDWG@data$LEVEL3_COD,
  Count      = "0"
)

# Merge data with TDWG level 3 names
TDWG_level3_df <- merge(
  TDWG_level3_df, Schoenus_worldwide_sums,
  by  = "LEVEL3_COD",
  all = TRUE
)
TDWG_level3_df[is.na(TDWG_level3_df)] <- 0
colnames(TDWG_level3_df) <- c("id", "Count_x", "Count_y")

# Merge fortified shapefile with count data for each species
TDWG_level3 <- fortify(TDWG, region = "LEVEL3_COD")
TDWG_level3_df <-merge(TDWG_level3, TDWG_level3_df, by = "id", all = TRUE)

# Plot maps --------------------------------------------------------------------

my_colours <- c(
  "#FFFFFF",
  "#FBF3F2",
  "#FAEFEE",
  "#F8EBEA",
  "#F7E7E6",
  "#F6E3E2",
  "#F5E0DE",
  "#F5E0DE",
  "#F0D0CD",
  "#ECC4C1",
  "#E8B9B5",
  "#E6B1AD",
  "#D2736B",
  "#C9584E",
  "#B31205"
)

# Get a vector of counts
count_vector <- sort(unique(TDWG_level3_df$Count_y))

worldwide_plot <- ggplot() +
  geom_polygon(data = TDWG_level3_df,
    aes(
      fill  = factor(Count_y),
      group = group,
      x     = long,
      y     = lat
    ),
    color   = "grey30",
    lwd     = 0.1
  ) +
  scale_fill_manual(name = "Count",
    values = my_colours,
    breaks = count_vector
  ) +
  theme_void()
worldwide_plot

# Format the map
worldwide_plot +
  theme(
    legend.key.width  = unit(0.04, "inch"),
    legend.key.height = unit(0.10, "inch"),
    legend.title      = element_text(size = 7),
    legend.text       = element_text(size = 5),
    legend.position   = c(0.035, 0.5)
  )

# Formate the map (v2)
worldwide_plot +
  theme(
    legend.key.width  = unit(0.1, "inch"),
    legend.key.height = unit(0.2, "inch"),
    legend.title      = element_text(size = 10),
    legend.text       = element_text(size = 8),
    legend.position   = c(0.035, 0.6)
  )
