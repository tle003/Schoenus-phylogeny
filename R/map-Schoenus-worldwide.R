# Tammy Elliott and Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(rgdal)      # Shapefile I/O
library(raster)     # Country borders, other GIS functions
library(tidyverse)  # Data wrangling, figures

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

# Set ggplot2 theme
theme_set(theme_bw())

# Define a red gradient manually for now
red_gradient <- c(
  "#FFFFFF", "#FBF3F2", "#FAEFEE", "#F8EBEA", "#F7E7E6",
  "#F6E3E2", "#F5E0DE", "#F5E0DE", "#F0D0CD", "#ECC4C1",
  "#E8B9B5", "#E6B1AD", "#D2736B", "#C9584E", "#B31205"
)

# Get a vector of counts
count_vector <- sort(unique(TDWG_level3_df$Count_y))

worldwide_plot <- ggplot() +
  geom_polygon(data = TDWG_level3_df,
    aes(x = long, y = lat, group = group, fill = factor(Count_y)),
    colour = "grey30",
    size   = 0.1
  ) +
  coord_equal() +
  scale_fill_manual(name = "No. species",
    values = red_gradient,
    breaks = count_vector
  ) +
  scale_x_continuous(breaks = seq(-180, 180, 60), limits = c(-180, 180)) +
  scale_y_continuous(breaks = seq(-60, 90, 30),   limits = c(-60, 90)) +
  labs(x = "Longitude (º)", y = "Latitude (º)")
