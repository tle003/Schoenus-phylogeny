# ...

# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(rgdal) # Shapefile I/O
library(sp)

library(tidyverse)
library(magrittr)

library(ggraph)
library(tidygraph)

# Import data ------------------------------------------------------------------

dispersal_means <- read_csv("data/dispersal_means.csv")
dispersal_SDs   <- read_csv("data/dispersal_SDs.csv")

df2matrix <- function(df) {
  dispersal_matrix <- as.matrix(df[, -1])
  rownames(dispersal_matrix) <- df[[1]]
  dispersal_matrix
}

dispersal_means %<>% df2matrix()
dispersal_SDs   %<>% df2matrix()

dispersal_means[dispersal_means <= 0.5] <- NA
dispersal_means[(dispersal_means - dispersal_SDs) <= 0] <- NA

dispersal_means %>%
  as_tbl_graph() %>%
  ggraph() +
    #geom_node_point(size = 15) +
    geom_node_text(aes(label = name), colour = "black") +
    geom_edge_parallel(
      aes(colour = weight),
      sep = unit(0.25, "inches"),
      arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")
    ) +
    scale_edge_colour_viridis(direction = -1, na.value = NA) +
    theme_void()

dispersal_means %>%
  as_tbl_graph() %>%
  ggraph(layout = tribble(
      ~x,   ~y, ~Regions,
    0.33, 1.00, "Holarctic",
    0.00, 0.50, "Neotropics",
    0.75, 0.50, "Tropical Asia",
    1.00, 0.50, "Pacific",
    0.50, 0.66, "Africa",
    0.50, 0.25, "Cape",
    1.00, 0.15, "New Zealand",
    0.75, 0.33, "Australia",
    0.75, 0.25, "Western Australia"
  )) +
    #geom_node_point(size = 15) +
    geom_node_text(aes(label = name), colour = "black") +
    geom_edge_parallel(
      aes(edge_width = weight %>%
        {case_when(
          . >= 17  ~ "17.4",
          . >=  8  ~  "8.4",
          . >=  4  ~  "4.1",
          . >=  1  ~  "1-3",
          . <   1  ~  "< 1"
        )} %>%
        factor(levels = c(
           "< 1",
           "1-3",
           "4.1",
           "8.4",
          "17.4"
        ))
      ),
      #sep = unit(0.25, "inches"),
      arrow = grid::arrow(length = unit(0.1, "inches"), type = "closed")
    ) +
    #scale_edge_colour_viridis(direction = -1, na.value = NA) +
    scale_edge_width_manual(values = c(0.25, 0.5, 1, 2.5, 5)) +
    theme_void()

TDWG <- readOGR("~/wgsrpd/level1", layer = "level1", stringsAsFactors = FALSE)

plot(TDWG)

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
