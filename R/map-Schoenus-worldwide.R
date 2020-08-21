# Tammy Elliott and Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(rgdal)      # Shapefile I/O
library(raster)     # Country borders, other GIS functions
library(tidyverse)  # Data wrangling, figures

# Import data ------------------------------------------------------------------

Schoenus_worldwide <- read.csv("data/occurence-data/Schoenus-worldwide-TDWG.csv")

TDWG <- readOGR("~/wgsrpd/level3", layer = "level3", stringsAsFactors = FALSE)

# Tidy data --------------------------------------------------------------------

colnames(Schoenus_worldwide)[1:4] <- c(
  "species",
  "authority",
  "source",
  "in_phylogeny"
)

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

# Merge counts (total, in phylogeny and not in phylogeny)
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

# Tidy data some more
TDWG_level3_df_tidy <- TDWG_level3_df %>%
  as_tibble() %>%
  dplyr::select(long, lat, group, Count) %>%
  rename(richness = Count) %>%
  mutate(richness = ifelse(richness == 0, NA, richness)) %>%
  mutate(richness =
    case_when(
      richness ==  1 ~ "1",
      richness <= 11 ~ "2-11",
      richness <= 21 ~ "12-21",
      richness == 45 ~ "45",
      richness == 62 ~ "62"
    ) %>%
    factor(levels = c("62", "45", "12-21", "2-11", "1"))
  )

# Plot maps --------------------------------------------------------------------

# Set ggplot2 theme
theme_set(theme_bw() + theme(panel.grid = element_blank()))

worldwide_plot <- ggplot() +
  geom_polygon(data = TDWG_level3_df_tidy,
    aes(x = long, y = lat, group = group, fill = richness, alpha = richness),
    colour = "grey30",
    size   = 0.1
  ) +
  annotate("text", label = "A", x = -180, y = 85, hjust = 1) +
  coord_equal() +
  scale_fill_viridis_d(na.translate = FALSE) +
  scale_alpha_manual(values = c(1, 1, 1, 0.75, 0.5)) +
  scale_x_continuous(
    breaks = seq(-180, 180, 60),
    limits = c(-180, 180)
  ) +
  scale_y_continuous(
    breaks = seq(-60, 90, 30),
    limits = c(-60, 90),
    labels = scales::label_math(expr = .x*"º")
  ) +
  guides(
    fill = guide_legend(
      title = "No. species",
      override.aes = list(alpha = c(1, 1, 1, 0.75, 0.5))
    ),
    alpha = FALSE
  ) +
  theme(
    axis.title.y      = element_blank(),
    axis.title.x      = element_blank(),
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    legend.text.align = 1
  )

proportion_sampled_plot <- ggplot() +
  geom_polygon(
    data = TDWG_level3_df %>%
      mutate(prop_in_phylogeny = ifelse(Count == 0, NA, prop_in_phylogeny)) %>%
      mutate(prop_in_phylogeny = prop_in_phylogeny * 100),
    aes(x = long, y = lat, group = group, fill = prop_in_phylogeny),
    colour = "grey30", size = 0.1
  ) +
  annotate("text", label = "B", x = -180, y = 85, hjust = 1) +
  coord_equal() +
  scale_fill_distiller(
    name      = "% species\nin phylogeny",
    palette   = "RdYlGn",
    direction = 1,
    na.value  = "white"
  ) +
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
