# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(rgdal) # Shapefile I/O
library(sp)

library(tidyverse)
library(magrittr)

library(ggraph)
library(tidygraph)

# Import data ------------------------------------------------------------------

dispersal_means <- read_csv("data/all_dispersals_counts_fromto_means.txt")
dispersal_SDs   <- read_csv("data/all_dispersals_counts_fromto_sds.txt")

df2matrix <- function(df) {
  dispersal_matrix <- as.matrix(df[, -1])
  rownames(dispersal_matrix) <- df[[1]]
  dispersal_matrix
}

dispersal_means %<>% df2matrix()
dispersal_SDs   %<>% df2matrix()

dispersal_means[dispersal_means <= 0.5] <- 0
dispersal_means[(dispersal_means - dispersal_SDs) <= 0.5] <- 0

sort(unique(as.vector(dispersal_means)))

dispersal_means %>%
  as_tbl_graph() %>%
  ggraph(layout = "linear", #tribble(
    #  ~x,   ~y, ~Regions,
    #0.33, 1.00, "Holarctic",
    #0.00, 0.50, "Neotropics",
    #0.75, 0.50, "Tropical Asia",
    #1.00, 0.50, "Pacific",
    #0.50, 0.66, "Africa",
    #0.50, 0.25, "Cape",
    #1.00, 0.15, "New Zealand",
    #0.75, 0.33, "Australia",
    #0.75, 0.25, "Western Australia")
  ) +
    geom_node_point(
      aes(colour = {
        case_when(
          name == "C"      ~ "Cape",
          name == "F"      ~ "Africa",
          name == "W"      ~ "Western Australia",
          name == "A"      ~ "Australia",
          name == "Z"      ~ "New Zealand",
          name == "N"      ~ "Neotropics",
          name == "P"      ~ "Pacific",
          name == "T"      ~ "Tropical Asia",
          name == "H"      ~ "Holarctic"
          ) %>%
        factor(levels = c(
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
      }),
      size = 5
    ) +
    geom_edge_arc2(
      aes(edge_width = weight %>%
        {case_when(
          . >= 16  ~   "16.05",
          . >=  7  ~    "7.30",
          TRUE     ~ "<= 4"
        )} %>%
        factor(levels = c(
          "<= 4",
             "7.30",
            "16.05"
        ))
      ),
      start_cap = circle(10, "points"),
      end_cap = circle(10, "points"),
      sep = unit(10, "points"),
      lineend = "round",
      arrow = grid::arrow(length = unit(0.1, "inches"), type = "open")
    ) +
    scale_edge_width_manual(values = c(0.25, 1, 2)) +
    scale_colour_brewer(palette = "Paired") +
    theme_void()

dispersal <-
  full_join(
    dispersal_means %>%
      as_tibble(rownames = "from") %>%
      gather(to, mean_n_dispersal_events, -from),
    dispersal_SDs %>%
      as_tibble(rownames = "from") %>%
      gather(to, SD_dispersal_events, -from)
  ) %>%
  mutate_at(vars(c("from", "to")), ~factor(., levels = c(
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
  filter(mean_n_dispersal_events > 0)
ggplot(dispersal) +
  aes(x = to, y = from, size = mean_n_dispersal_events) +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggplot(dispersal) +
  aes(x = from, y = mean_n_dispersal_events, colour = to, group = to, ) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(
      ymin = mean_n_dispersal_events - SD_dispersal_events,
      ymax = mean_n_dispersal_events + SD_dispersal_events
    ),
    width = 0,
    position = position_dodge(width = 0.5),
    na.rm = TRUE
  ) +
  scale_colour_brewer(palette = "Paired", drop = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
