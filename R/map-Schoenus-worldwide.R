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
TDWG_level3_df <- merge(TDWG_level3, TDWG_level3_df, by = "id", all = TRUE)

# Tidy data some more
TDWG_level3_df <- TDWG_level3_df %>%
  as_tibble() %>%
  select(long, lat, group, Count_y) %>%
  rename(richness = Count_y) %>%
  mutate(richness = ifelse(richness == 0, NA, richness)) %>%
  mutate(richness =
    case_when(
      richness ==  1 ~ "1",
      richness <= 11 ~ "2-11",
      richness <= 21 ~ "12-21",
      richness == 45 ~ "45",
      richness == 62 ~ "62"
    ) %>%
    factor(levels = c("1", "2-11", "12-21", "45", "62"))
  )

# Plot maps --------------------------------------------------------------------

# Set ggplot2 theme
theme_set(theme_bw())

worldwide_plot <- ggplot() +
  geom_polygon(data = TDWG_level3_df,
    aes(x = long, y = lat, group = group, fill = richness),
    colour = "grey30",
    size   = 0.1
  ) +
  coord_equal() +
  scale_fill_grey(
    name = "No. species",
    start = 0.9, end = 0.1,
    na.translate = FALSE
  ) +
  scale_x_continuous(breaks = seq(-180, 180, 60), limits = c(-180, 180)) +
  scale_y_continuous(breaks = seq(-60, 90, 30),   limits = c(-60, 90)) +
  labs(x = "Longitude (º)", y = "Latitude (º)")
