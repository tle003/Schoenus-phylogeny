# Map Schoenus (Cyperaceae, Tribe Schoeneae) species richness in
# quarter degree grid cells in Australia, New Caledonia and New Zealand

# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(rgdal)      # Shapefile I/O
library(raster)     # Country borders, other GIS functions
library(tidyverse)  # Data wrangling

# Import data ------------------------------------------------------------------

# Schoenus occurrence records & neat column headings
# (NOTE: not in repository)
occ <- read_csv("Schoenus-occurence-data/Schoenus-Australia-records-2020-06-04/Schoenus-Australia-records-2020-06-04.csv")
headings <- read_csv("Schoenus-occurence-data/Schoenus-Australia-records-2020-06-04/headings.csv")

# Larsen et al. (2009) QDGC for Oceania
QDGC <- readOGR("qdgc_oceania", layer = "qdgc_02_oceania")
# Have a look (WARNING: takes a while)
#plot(QDGC)

# Australia & New Zealand borders
border_AUS <- readOGR("border_AUS", layer = "GID_0")
border_NZL <- readOGR("border_NZL", layer = "GID_0")
border_NCL <- readOGR("border_NCL", layer = "GID_0")
border <- rbind(border_AUS, border_NZL, border_NCL)
# Have a look (WARNING: takes a while)
#plot(border)
# Lots of stray islands...

# Crop to focus on the mainlands (NOTE: takes a while)
main_extent <-extent(110, 180, -50, 0)
border_cropped <- crop(border, main_extent)

# Tidy data --------------------------------------------------------------------

# Check that headings' descriptions are in the same order as the columns
colnames(occ) == headings[[1]]
# Use syntactically valid column names
colnames(occ) <- headings[[4]]

# Check that no super-specific taxa in data
occ %>%
  pull(rank) %>%
  unique()
# There are! Some records are only ID-ed to the genus-level

occ_tidy <- occ %>%
  # Remove genus-only records
  filter(rank != "genus") %>%
  # Remove records without GPS coordinates
  filter(!is.na(latitude) & !is.na(longitude)) %>%
  # Take species-level of all records
  # (i.e. sub-specific taxa are treated as occ. of its species)
  select(species, latitude, longitude)
# Have a look
occ_tidy

# Collate species richness in grid-cells ---------------------------------------

occ_spdf <- SpatialPointsDataFrame(
  coords      = occ_tidy[, c("longitude", "latitude")],
  data        = occ_tidy[, "species"],
  proj4string = QDGC@proj4string
)
# Have a look
#plot(occ_spdf)

occ_spdf$QDGC <- occ_spdf %>%
  over(QDGC) %>%
  pull(qdgc) %>%
  as.character()

# Crop to focus on the mainlands (NOTE: takes a while)
occ_spdf_cropped <- crop(occ_spdf, main_extent)

richness <- occ_spdf_cropped %>%
  as_tibble() %>%
  group_by(QDGC) %>%
  summarise(richness = length(unique(species))) %>%
  filter(!is.na(QDGC))  # NOTE: only 1 species occurs in an NA QDGC

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
# 34 is richest square

# Plot maps --------------------------------------------------------------------

# Set ggplot2 theme
theme_set(theme_bw())

oceania_plot <- ggplot(richness) +
  aes(QDGC_lon, QDGC_lat, fill = richness) +
  geom_tile() +
  geom_polygon(
    data = border_cropped,
    aes(x = long, y = lat, group = group),
    colour = "black", fill = NA, size = 0.25
  ) +
  coord_equal() +
  scale_x_continuous(breaks = c(120, 140, 160, 180)) +
  scale_y_continuous(breaks = c(0, -20, -40)) +
  scale_fill_viridis_c(
    name = bquote(italic("Schoenus")*" species richness"),
    direction = -1,
    breaks = c(5, 15, 25, 34),
    limits = c(0, 34)
  ) +
  xlab("Longitude (º)") +
  theme(axis.title.y = element_blank())
