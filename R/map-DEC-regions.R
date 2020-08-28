# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

# GIS:
library(rgdal)
library(sp)
library(ggfortify)

# Programming:
library(tidyverse)
library(magrittr)
library(glue)

# Import data ------------------------------------------------------------------

# TDWG shapefile(s):
level_names <- glue("level{1:4}")
TDWG <- map(level_names, ~readOGR(
  glue("~/wgsrpd/{.}"),
  layer = .,
  stringsAsFactors = FALSE
))
names(TDWG) <- level_names

# Tidy data --------------------------------------------------------------------

# Rename TDWG columns:
names(TDWG$level1) <- c("level1_code", "level1_name")
names(TDWG$level2) <- c(
  "level2_code",
  "level1_code",
  "level1_name",
  "level2_name"
)
names(TDWG$level3) <- c(
  "level3_name", "level3_code",
  "level2_code", "level1_code"
)
names(TDWG$level4) <- c(
  "ISO","level4_name", "level4_code", "level4_2",
  "level3_code", "level2_code", "level1_code"
)

# Reorder and/or drop some TDWG columns:
TDWG$level1 <- TDWG$level1[c(
  "level1_name",
  "level1_code"
)]
TDWG$level2 <- TDWG$level2[c(
  "level2_name",
  "level2_code", "level1_code"
)]
TDWG$level3 <- TDWG$level3[c(
  "level3_name",
  "level3_code", "level2_code", "level1_code"
)]
TDWG$level4 <- TDWG$level4[c(
  "level4_name",
  "level4_code", "level3_code", "level2_code", "level1_code"
)]

# FIXME:
#TDWG_fortified <- TDWG %>%
#  map(~fortify(., region = names(.)[[1]]))
## Error in rgeos::gUnaryUnion(spgeom = SpP, id = IDs) :
## TopologyException: Input geom 1 is invalid: Ring Self-intersection
## at or near point 7.5312446199999998 4.6023051600000002 at
## 7.5312446199999998 4.6023051600000002

for (i in seq_along(TDWG$level4$level3_code)) {
  TDWG$level4$level3_name[[i]] <- unique(TDWG$level3$level3_name[
    TDWG$level3$level3_code == TDWG$level4$level3_code[[i]]
  ])
  TDWG$level4$level2_name[[i]] <- unique(TDWG$level2$level2_name[
    TDWG$level2$level2_code == TDWG$level4$level2_code[[i]]
  ])
  TDWG$level4$level1_name[[i]] <- unique(TDWG$level1$level1_name[
    TDWG$level1$level1_code == TDWG$level4$level1_code[[i]]
  ])
}

TDWG$level4$region <- TDWG$level4 %$% {case_when(
    (.$level3_name == "Cape Provinces")     ~ "Cape",

    (.$level1_name == "AFRICA")
  & (.$level2_name != "Northern Africa")
  & (.$level3_name != "Cape Province")      ~ "Africa",

    (.$level4_name == "Western Australia")  ~ "Western Australia",

    (.$level2_name == "Australia")
  & (.$level4_name != "Western Australia")  ~ "Australia",

    (.$level2_name == "New Zealand")
  & (.$level4_name != "Chatham Is.")        ~ "New Zealand",

    (.$level1_name == "SOUTHERN AMERICA")
  | (.$level2_name == "Mexico")             ~ "Neotropics",

    (.$level1_name == "PACIFIC")            ~ "Pacific",

    (.$level1_name == "ASIA-TROPICAL")      ~ "Tropical Asia",

      ((.$level1_name %in% c(
        "NORTHERN AMERICA",
        "EUROPE",
        "ASIA-TEMPERATE"))
    | (.$level2_name == "Northern Africa"))
  & (.$level2_name != "Mexico")             ~ "Holarctic"
)}

TDWG$level4$all_level_names <- paste(sep = "_",
  TDWG$level4$level1_name,
  TDWG$level4$level2_name,
  TDWG$level4$level3_name,
  TDWG$level4$level4_name,
  TDWG$level4$region
)

TDWG_level4_tidy <- fortify(TDWG$level4, region = "all_level_names")
TDWG_level4_DEC_regions <- TDWG_level4_tidy %>%
  as_tibble() %>%
  separate(id, sep = "_",
    into = c("level1_name", "level2_name", "level3_name", "level4_name", "region")
  ) %>%
  mutate(region = factor(region, levels = c(
    "Cape",
    "Africa",
    "Western Australia",
    "Australia",
    "New Zealand",
    "Neotropics",
    "Pacific",
    "Tropical Asia",
    "Holarctic"
  ))) %>%
  filter(!is.na(region))

# Plot -------------------------------------------------------------------------

my_palette <- scales::brewer_pal(palette = "Paired")(
  n = length(c(
    "Cape",
    "Africa",
    "Western Australia",
    "Australia",
    "New Zealand",
    "Neotropics",
    "Pacific",
    "Tropical Asia",
    "Holarctic"
  ))
)
# Darken purple
my_palette[[9]] <- "#AB71C7"

DEC_regions_plot <- ggplot() +
  geom_polygon(
    data = TDWG_level4_DEC_regions,
    aes(
      x = long, y = lat, group = group,
      fill   = region,
      colour = region
    ),
    size = 0.75
  ) +
  scale_fill_manual(values = my_palette) +
  scale_colour_manual(values = my_palette) +
  coord_equal() +
  theme_void() +
  theme(legend.position = "none")

# Save plot --------------------------------------------------------------------

ggsave(
  "figures/DEC_regions_plot.pdf",
  DEC_regions_plot,
  width = 10, height = 5
)
ggsave(
  "figures/DEC_regions_plot.png",
  DEC_regions_plot,
  width = 10, height = 5, dpi = 300
)
