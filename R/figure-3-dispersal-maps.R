# ...

# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(rgdal) # Shapefile I/O
library(tidyverse)

# Import data ------------------------------------------------------------------

dispersal_matrix <- read_tsv("data/all_dispersals_counts_fromto_means.txt")
dispersal_matrix_codes <- read_csv("data/all_dispersals_counts_fromto_means_codes.csv")

TDWG <- readOGR("~/wgsrpd/level1", layer = "level1", stringsAsFactors = FALSE)

plot(TDWG$LEVEL1_NAM)

theme_set(theme_classic())
dispersal_matrix %>%
  gather(to, mean_n_dispersal_events, -X1) %>%
  rename(from = X1) %>%
  mutate_at(vars(c("from", "to")), ~map_chr(.,
    ~dispersal_matrix_codes$region[dispersal_matrix_codes$code == .]
  )) %>%
  filter(mean_n_dispersal_events > 0) %>%
  ggplot() +
    aes(x = to, y = from, size = mean_n_dispersal_events) +
    geom_point()

dispersal_matrix %>%
  {
    y <- as.matrix(.[, -1])
    rownames(y) <- .[[1]]
    y
  } %>%
  t() %>%
  round() %>%
  diagram::plotmat()

library(network)
devtools::install_github("briatte/ggnet")
library(ggnet)
dispersal_matrix %>%
  {
    y <- as.matrix(.[, -1])
    rownames(y) <- .[[1]]
    y
  } %>%
  network(directed = TRUE) %>%
  ggnet2(label = "vertex.names")
