# ...

# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(tidyverse)

# Import data ------------------------------------------------------------------

dispersal_matrix <- read_tsv("data/all_dispersals_counts_fromto_means.txt")
dispersal_matrix_codes <- read_csv("data/all_dispersals_counts_fromto_means_codes.csv")

theme_set(theme_classic())
dispersal_matrix %>%
  gather(to, rate, -X1) %>%
  rename(from = X1) %>%
  mutate_at(vars(c("from", "to")), ~map_chr(.,
    ~dispersal_matrix_codes$region[dispersal_matrix_codes$code == .]
  )) %>%
  filter(rate > 0) %>%
  ggplot() +
    aes(x = to, y = from, size = rate) +
    geom_point()
