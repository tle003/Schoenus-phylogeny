# ...

# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(rgdal) # Shapefile I/O
library(tidyverse)

# Import data ------------------------------------------------------------------

dispersal_means <- read_csv("data/dispersal_means.csv")
dispersal_SDs   <- read_csv("data/dispersal_SDs.csv")

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

library(ggraph)
library(tidygraph)
dispersal_matrix %>%
  {
    y <- as.matrix(.[, -1])
    rownames(y) <- .[[1]]
    y[y == 0] <- NA
    y
  } %>%
  as_tbl_graph() %>%
  ggraph() +
    geom_node_point(size = 15) +
    geom_node_text(aes(label = name), colour = "white") +
    geom_edge_parallel(
      aes(colour = weight),
      sep = unit(0.25, "inches"),
      arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")
    ) +
    scale_edge_colour_viridis(direction = -1, na.value = NA)
