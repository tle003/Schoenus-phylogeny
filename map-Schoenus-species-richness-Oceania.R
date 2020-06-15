# Map Schoenus (Cyperaceae, Tribe Schoeneae) species richness in
# quarter degree grid cells in Oceania

# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(rgdal)      # Shapefile I/O
library(sp)         # Map projections, other GIS functions
library(raster)     # Country borders, other GIS functions
library(tidyverse)  # Data wrangling

# Import data ------------------------------------------------------------------
# (NOTE: not in repository)

# Schoenus occurrence records & neat column headings
occ <- read_csv("Schoenus-Australia-records-2020-06-04/Schoenus-Australia-records-2020-06-04.csv")
headings <- read_csv("Schoenus-Australia-records-2020-06-04/headings.csv")
# Larsen et al. (2009) QDGC for Oceania
QDGC <- readOGR("qdgc_oceania", layer = "qdgc_02_oceania")

# Oceanian countries' borders
border_AUS <- getData("GADM", country = "AUS", level = 0)  # Australia
border_IDN <- getData("GADM", country = "IDN", level = 0)  # Indonesia
border_MYS <- getData("GADM", country = "MYS", level = 0)  # Malaysia
border_NZL <- getData("GADM", country = "NZL", level = 0)  # New Zealand
border_PNG <- getData("GADM", country = "PNG", level = 0)  # Papua New Guinea
border_PHL <- getData("GADM", country = "PHL", level = 0)  # Phillipines
border <- rbind(
  border_AUS,
  border_IDN,
  border_MYS,
  border_NZL,
  border_PNG,
  border_PHL
)

# Crop to focus on the mainlands (NOTE: takes a while)
border_cropped <- crop(border, QDGC)
# Have a look
#plot(border_cropped)

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

richness <- occ_spdf %>%
  as_tibble() %>%
  group_by(QDGC) %>%
  summarise(richness = length(unique(species))) %>%
  filter(!is.na(QDGC))

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
# FIXME/NOTE: NA QDGC (with richness 18) might be islands?

# Note:
max(richness$richness)
# 38 is richest square

# Plot maps --------------------------------------------------------------------

# Set ggplot2 theme
theme_set(theme_bw())

richness_map <- ggplot(richness) +
  aes(QDGC_lon, QDGC_lat, fill = richness) +
  geom_tile() +
  geom_polygon(
    data = border,
    aes(x = long, y = lat, group = group),
    colour = "black", fill = NA, size = 0.25
  ) +
  coord_equal() +
  scale_x_continuous(breaks = c(110, 120, 130, 140, 150, 160, 170, 180)) +
  labs(x = "Longitude (º)", y = "Latitude (º)") +
  theme(legend.position = "top")

# Make maps in both colour & black-&-white
richness_map_colour <- richness_map +
  scale_fill_viridis_c(
    name = bquote(italic("Schoenus")*" species richness"),
    direction = -1,
    breaks = c(5, 15, 25, 34)
  )
richness_map_bw <- richness_map +
  scale_fill_gradient(
    name = bquote(italic("Schoenus")*" species richness"),
    low = "grey90", high = "black",
    breaks = c(5, 15, 25, 34)
  )

# Save maps --------------------------------------------------------------------
# (as PDFs & PNGs)

ggsave(
  "maps/Schoenus-species-richness-Oceania_colour.pdf",
  richness_map_colour,
  width = 6, height = 4
)
ggsave(
  "maps/Schoenus-species-richness-Oceania_colour.png",
  richness_map_colour,
  width = 6, height = 4,
  dpi = 300
)

ggsave(
  "maps/Schoenus-species-richness-Oceania_bw.pdf",
  richness_map_bw,
  width = 6, height = 4
)
ggsave(
  "maps/Schoenus-species-richness-Oceania_bw.png",
  richness_map_bw,
  width = 6, height = 4,
  dpi = 300
)
