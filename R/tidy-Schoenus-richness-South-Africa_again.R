# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(rgdal)      # Shapefile I/O
library(raster)     # Country borders, other GIS functions
library(tidyverse)  # Data wrangling

# Import data ------------------------------------------------------------------

# South African Schoenus richness in QDGC (from Tammy Elliott)
richness <- read_csv("data/Schoenus-South-Africa-richness.csv")

# Larsen et al. (2009) QDGC for South Africa
QDGC <- readOGR("data/shapefiles/qdgc_zaf", layer = "qdgc_02_zaf")

# South African border
border_ZAF <- getData("GADM", country = "ZAF", level = 0)

# Tidy data --------------------------------------------------------------------

# Make data frame and QDGC shapefile codes comparable
QDGC$qdgc <- QDGC$qdgc %>%
  str_remove("E0") %>%
  str_remove("S")

# Add richness of QDGCs from data frame to QDGC shapefile's data
QDGC$richness <- NA
for (i in 1:nrow(QDGC)) {
  if (QDGC$qdgc[[i]] %in% richness$QDGC) {
    QDGC$richness[[i]] <- richness$richness[
      richness$QDGC == QDGC$qdgc[[i]]
    ]
  }
}

# Figure out which QDGC is outside of South Africa's borders -------------------

# Query the QDGCs for their being within the South African border
QDGC_over_ZAF <- QDGC %over% border_ZAF
# Save that into the QDGC shapefile
QDGC$ZAF <- QDGC_over_ZAF$GID_0
# Ammend it such that if we have richness data for a QDGC,
# but it is outside the border,
# it is distinguished from other QDGC outside the border
QDGC@data <- QDGC@data %>%
  mutate(ZAF = case_when(
    !is.na(richness) & (ZAF == "ZAF") ~ 1,
    !is.na(richness) & is.na(ZAF)     ~ 2,
     is.na(richness) & is.na(ZAF)     ~ 3,
  ))

# Merge fortified shapefile with count data for each species
QDGC_qdgc <- fortify(QDGC, region = "qdgc")
QDGC_over_ZAF2 <- merge(
  QDGC_qdgc,
  QDGC %>%
    as_tibble() %>%
    rename(id = qdgc) %>%
    select(-lon, -lat),
  by = "id",
  all = TRUE
)
# Have a look
ggplot() +
  geom_polygon(
    data = QDGC_over_ZAF2,
    aes(x = long, y = lat, group = group, fill = as.factor(ZAF))
  ) +
  scale_fill_viridis_d(na.value = "white") +
  geom_polygon(
    data = border_ZAF,
    aes(x = long, y = lat, group = group),
    colour = "white", fill = NA
  )
# Works!

# Get the QDGC-code of that "bad" QDGC
QDGC@data %>%
  filter(ZAF == 2) %>%
  select(qdgc, richness)
# Result:
##     qdgc richness
## 1 1834CD        1

# Ammend richness dataset for South Africa accordingly -------------------------

# I.e. removing this erroneous QDGC
richness %>%
  filter(QDGC != "1834CD") %>%
  write_csv("data/Schoenus-South-Africa-richness.csv")
