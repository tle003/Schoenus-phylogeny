# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(rgdal)      # Shapefile I/O
library(tidyverse)  # Data wrangling

# Import data ------------------------------------------------------------------

# Schoenus occurrence records & neat column headings
# (NOTE: not in repository)
occ <- read_csv("data/occurence-data/Schoenus-Australia-records-2020-06-04/Schoenus-Australia-records-2020-06-04.csv")
headings <- read_csv("data/occurence-data/Schoenus-Australia-records-2020-06-04/headings.csv")

# Larsen et al. (2009) QDGC for Oceania
QDGC_AUS <- readOGR("data/shapefiles/qdgc_aus", layer = "qdgc_02_aus")
QDGC_NZL <- readOGR("data/shapefiles/qdgc_nzl", layer = "qdgc_02_nzl")
QDGC <- rbind(QDGC_AUS, QDGC_NZL)
# Have a look (WARNING: takes a while)
#plot(QDGC)

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
#ggplot(occ_tidy) +
#  aes(longitude, latitude) +
#  geom_point(alpha = 0.25) +
#  xlim(110, 180) +
#  coord_equal()

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
main_extent <-extent(110, 180, -50, -10)
occ_spdf_cropped <- crop(occ_spdf, main_extent)

richness <- occ_spdf_cropped %>%
  as_tibble() %>%
  group_by(QDGC) %>%
  summarise(richness = length(unique(species))) %>%
  filter(!is.na(QDGC))  # NOTE: 1x NA QDGC w/ richness 8

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

# NOTE:
max(richness$richness)
# 34 is richest square

# Export data ------------------------------------------------------------------

write_csv(richness, "data/Schoenus-Oceania-richness.csv")
