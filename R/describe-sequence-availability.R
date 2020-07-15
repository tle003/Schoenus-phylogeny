# First attempt at collating available sequence data
# for Schoenus (Cyperaceae, Tribe Schoeneae)

# Merges Tammy and Ruan's records of available sequence data for Schoeneae
# (files in folder data/sequence-availability)

# Ruan van Mazijk, 2020

# Load packages ----------------------------------------------------------------

library(tidyverse)
library(magrittr)

# Import Tammy and Ruan's data -------------------------------------------------

tle <-
  read_csv("data/sequence-availability/2020-04-03_specimen-tot-nodups_tle.csv")
rvm <-
  read_csv("data/sequence-availability/2020-04-02_sequences-tidy-by-voucher_rvm.csv")

# Tidy both our datasets -------------------------------------------------------
# (make columns names the same, etc.)

tle %<>%
  set_colnames(c(
    "taxon", "specimen", "voucher",
    "ETS", "ITS", "rbcL", "rps16", "trnLF", "n_markers_available"
  )) %>%
  select(-specimen, -n_markers_available)

rvm %<>%
  mutate(taxon = str_replace_all(taxon, "_", " ")) %>%
  select(-n_markers_available)

# Merge datasets ---------------------------------------------------------------

sequences_tidy <-
  full_join(tle, rvm) %>%
  # NOTE: Lots of double rows due to inconsistent voucher name standards
  # Ammended here:
  mutate(voucher = voucher %>%
    str_replace_all("TV", "Verboom ") %>%
    str_replace_all("MM", "Muasya ") %>%
    {ifelse(str_detect(., "^[0-9]{3}$"), paste0("TE2016_", .), .)} %>%
    {ifelse(str_detect(., "^[0-9]{4}$"), paste("Muasya", .), .)}
  ) %>%
  distinct() %>%
  # Score markers available as 1, unavailable as 0
  mutate_at(vars(c("ETS", "ITS", "rbcL", "rps16", "trnLF")),
    ~ifelse(is.na(.), 0, .)
  ) %>%
  group_by(taxon, voucher) %>%
  summarise_at(vars(c("ETS", "ITS", "rbcL", "rps16", "trnLF")),
    ~ifelse(sum(.) > 0, 1, 0)
  ) %>%
  # Sum up number of markers available per voucher
  mutate(n_markers = ETS + ITS + rps16 + trnLF + rbcL) %>%
  arrange(taxon, voucher)
write_csv(sequences_tidy, "data/sequence-availability/sequences-by-voucher.csv")

# Recalculate n_markers by taxa
sequences_tidy_taxa <- sequences_tidy %>%
  select(-n_markers, -voucher) %>%
  group_by(taxon) %>%
  summarise_at(vars(c("ETS", "ITS", "rbcL", "rps16", "trnLF")),
    ~ifelse(sum(.) > 0, 1, 0)
  ) %>%
  arrange(taxon) %>%
  mutate(n_markers = ETS + ITS + rbcL + rps16 + trnLF)
write_csv(sequences_tidy_taxa, "data/sequence-availability/sequences-by-taxon.csv")

# Check for specimens (vouchers) that have differing dets in the 2 matrices
# (to be continued)
sequences_tidy_det_check <- sequences_tidy %>%
  ungroup() %>%
  select(taxon, voucher) %>%
  distinct() %>%
  arrange(voucher, taxon) %>%
  group_by(voucher) %>%
  summarise(n_dets = n(), taxa = paste(unique(taxon), collapse = ", ")) %>%
  arrange(desc(n_dets)) %>%
  filter(n_dets > 1, voucher != "?")
write_csv(sequences_tidy_det_check, "data/sequence-availability/dets-to-check.csv")

# Useful statistics ------------------------------------------------------------

sequences_tidy %>%
  # Recalculate n_markers by voucher specimen **only**
  ungroup() %>%
  select(-taxon) %>%
  distinct() %>%
  select(-n_markers) %>%
  group_by(voucher) %>%
  summarise_at(vars(c("ETS", "ITS", "rbcL", "rps16", "trnLF")),
    ~ifelse(sum(.) > 0, 1, 0)
  ) %>%
  arrange(voucher) %>%
  mutate(n_markers = ETS + ITS + rps16 + trnLF + rbcL) %T>%
  {print(nrow(.))} %>%  # 468 **specimens** sequenced
  filter(n_markers > 1) %>%
  nrow()  # 160 of which have > 1 marker

sequences_tidy_taxa %T>%
  {print(nrow(.))} %>%  # ca. 239 **taxa** sequenced
  filter(n_markers > 1) %>%
  nrow()  # 130 of which have > 1 marker

# Visualise sequence availability by taxon -------------------------------------

theme_set(theme_classic())

marker_names_ordered <- sequences_tidy_taxa %>%
  select(ETS, ITS, rbcL, rps16, trnLF) %>%
  map_dbl(sum) %>%
  sort(decreasing = TRUE) %>%
  names()

sequences_tidy_taxa %<>%
  select(-n_markers) %>%
  gather(marker_name, available, ETS:trnLF) %>%
  mutate(
    marker_name = factor(marker_name, levels = marker_names_ordered),
    available   = (available == 1),
    genus       = str_extract(taxon, "[a-zA-Z]+")
  )

sequences_plot <- ggplot(sequences_tidy_taxa) +
  aes(marker_name, taxon, fill = genus, alpha = available) +
  geom_tile() +
  labs(x = "Marker sequenced?") +
  scale_alpha_manual(values = c(0, 1)) +
  theme(
    axis.text.y     = element_text(hjust = 0),
    axis.title.y    = element_blank(),
    legend.position = "none"
  )

# Save plot --------------------------------------------------------------------

ggsave(
  "figures/sequences-by-taxon.pdf",
  sequences_plot,
  width = 5, height = 30
)
