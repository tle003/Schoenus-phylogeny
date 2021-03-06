# Map Schoenus (Cyperaceae, Tribe Schoeneae) species richness and
# sampling effort (Re: our phylogenetic analyses) in TDWG areas worldwide

# Tammy Elliott and Ruan van Mazijk, 2021

# Load packages ----------------------------------------------------------------

library(rgdal)      # Shapefile I/O
library(raster)     # Country borders, other GIS functions
library(tidyverse)  # Data wrangling, figures

# Import data ------------------------------------------------------------------

Schoenus_worldwide <- read.csv("data/occurence-data/Schoenus-worldwide-TDWG.csv")

TDWG <- readOGR("~/wgsrpd/level3", layer = "level3", stringsAsFactors = FALSE)
# NOTE: this requires an external repo's data
# (namely the TDWG areas from the WGSRPD)

# Tidy data --------------------------------------------------------------------

# Rename some columns
colnames(Schoenus_worldwide)[1:4] <- c(
  "species",
  "authority",
  "source",
  "in_phylogeny"
)

# Give species names as rownames
rownames(Schoenus_worldwide) <- Schoenus_worldwide$species

# Calculate rowsums to get a vector with values per region
Schoenus_worldwide_sums <- as.data.frame(colSums(
  Schoenus_worldwide[, 5:ncol(Schoenus_worldwide)]
))
Schoenus_worldwide_sums <- data.frame(
  LEVEL3_COD = rownames(Schoenus_worldwide_sums),
  Count      = Schoenus_worldwide_sums
)
colnames(Schoenus_worldwide_sums) <- c("LEVEL3_COD", "Count")

# Calculate rowsums but only species in phylogeny
Schoenus_in_phylogeny <- Schoenus_worldwide %>%
  mutate(in_phylogeny = as.logical(in_phylogeny)) %>%
  filter(in_phylogeny)
Schoenus_in_phylogeny_sums <- as.data.frame(colSums(
  Schoenus_in_phylogeny[, 5:ncol(Schoenus_in_phylogeny)]
))
Schoenus_in_phylogeny_sums <- data.frame(
  LEVEL3_COD = rownames(Schoenus_in_phylogeny_sums),
  Count      = Schoenus_in_phylogeny_sums
)
colnames(Schoenus_in_phylogeny_sums) <- c("LEVEL3_COD", "Count_in")

# Calculate rowsums but only species not in phylogeny
Schoenus_not_in_phylogeny <- Schoenus_worldwide %>%
  mutate(in_phylogeny = as.logical(in_phylogeny)) %>%
  filter(!in_phylogeny)
Schoenus_not_in_phylogeny_sums <- as.data.frame(colSums(
  Schoenus_not_in_phylogeny[, 5:ncol(Schoenus_not_in_phylogeny)]
))
Schoenus_not_in_phylogeny_sums <- data.frame(
  LEVEL3_COD = rownames(Schoenus_not_in_phylogeny_sums),
  Count      = Schoenus_not_in_phylogeny_sums
)
colnames(Schoenus_not_in_phylogeny_sums) <- c("LEVEL3_COD", "Count_not_in")

# Merge counts (total, no. in phylogeny, and no. not in phylogeny)
Schoenus_worldwide_sums <-
  full_join(Schoenus_in_phylogeny_sums, Schoenus_not_in_phylogeny_sums) %>%
  full_join(Schoenus_worldwide_sums) %>%
  mutate(
    prop_in_phylogeny     = Count_in     / (Count_in + Count_not_in),
    prop_not_in_phylogeny = Count_not_in / (Count_in + Count_not_in)
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
TDWG_level3_df <- TDWG_level3_df[, -2]
colnames(TDWG_level3_df)[[1]] <- "id"
colnames(TDWG_level3_df)[[4]] <- "Count"

# Merge fortified shapefile with count data for each species
TDWG_level3 <- fortify(TDWG, region = "LEVEL3_COD")
TDWG_level3_df <- merge(TDWG_level3, TDWG_level3_df, by = "id", all = TRUE)

# Tidy data some more for plotting
TDWG_level3_df_tidy <- TDWG_level3_df %>%
  as_tibble() %>%
  dplyr::select(long, lat, group, Count, prop_in_phylogeny) %>%
  rename(richness = Count) %>%
  mutate(
    richness          = ifelse(richness == 0, NA, richness),
    prop_in_phylogeny = ifelse(richness == 0, NA, prop_in_phylogeny)
  ) %>%
  mutate(
    # Discretise richness values
    richness =
      case_when(
        richness ==  1 ~ "1",
        richness <= 10 ~ "2-10",
        richness <= 21 ~ "11-21",
        richness == 46 ~ "46",
        richness == 61 ~ "61"
      ) %>%
      factor(levels = c("61", "46", "11-21", "2-10", "1")),
    # Tidy sampling proportion values
    prop_in_phylogeny = prop_in_phylogeny * 100
  )

# Plot maps --------------------------------------------------------------------

# Set ggplot2 theme
theme_set(theme_bw() + theme(panel.grid = element_blank()))

# Note, the panels of Figure 2 are as follows:
#
#     A--------------------\
#     |                    |
#     | richness           |
#     |                    |
#     \--------------------/
#     B--------------------\
#     |                    |
#     | proportion sampled |
#     |                    |
#     \--------------------/
#     C--------\ D---------\
#     |        | |         |
#     | S Afr  | | W Aus   |
#     |        | |         |
#     \--------/ \---------/
#
# Consequently, will omit y-axis labels on the bottom of panel A

# Panel A
richness_map <- ggplot() +
  # Plot TDWG polygons
  geom_polygon(data = TDWG_level3_df_tidy,
    aes(
      x     = long,
      y     = lat,
      group = group,
      fill  = richness,  # Coloured by richness
      alpha = richness   # (for NAs to be transparent)
    ),
    colour = "grey30",
    size   = 0.1  # (= line weight)
  ) +
  # Panel label
  annotate("text", label = "A", x = -180, y = 85, hjust = 1) +
  # Colour scheme customisation
  scale_fill_viridis_d(na.translate = FALSE) +
  scale_alpha_manual(values = c(1, 1, 1, 0.75, 0.5)) +
  # Lat/lon scale customisation
  coord_equal() +
  scale_x_continuous(
    breaks = seq(-180, 180, 60),
    limits = c(-180, 180)
    # (no labels needed as omitting this axis)
  ) +
  scale_y_continuous(
    breaks = seq(-60, 90, 30),
    limits = c(-60, 90),
    labels = scales::label_math(expr = .x*"º")
  ) +
  # Manually make legend (Re: fiddle with alpha for NAs)
  guides(
    fill = guide_legend(
      title = "No. species",
      override.aes = list(alpha = c(1, 1, 1, 0.75, 0.5))
    ),
    alpha = FALSE
  ) +
  theme(
    # Omit x-axis labels etc.
    axis.ticks.x      = element_blank(),
    axis.text.x       = element_blank(),
    axis.title.x      = element_blank(),
    # But keep y-axis! (except for title)
    axis.title.y      = element_blank(),
    legend.text.align = 1
  )

# Panel B
proportion_sampled_map <- ggplot() +
  geom_polygon(data = TDWG_level3_df_tidy,
    # Plot TDWG polygons
    aes(
      x     = long,
      y     = lat,
      group = group,
      fill  = prop_in_phylogeny  # Coloured by sampling proportion
    ),
    colour = "grey30",
    size = 0.1  # (= line weight)
  ) +
  # Panel label
  annotate("text", label = "B", x = -180, y = 85, hjust = 1) +
  # Colour scheme customisation
  scale_fill_distiller(
    name      = "% species\nin phylogeny",
    palette   = "RdYlGn",
    direction = 1,
    na.value  = "white"
  ) +
  # Lat/lon scale customisation
  coord_equal() +
  scale_x_continuous(
    breaks = seq(-180, 180, 60),
    limits = c(-180, 180),
    labels = scales::label_math(expr = .x*"º")
  ) +
  scale_y_continuous(
    breaks = seq(-60, 90, 30),
    limits = c(-60, 90),
    labels = scales::label_math(expr = .x*"º")
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
  )
