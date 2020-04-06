library(tidyverse)
library(magrittr)

# Import Tammy and Ruan's data
tle <- read_csv("2020-04-03_specimen-tot-nodups_tle.csv")
rvm <- read_csv("2020-04-02_sequences-tidy-by-voucher_rvm.csv")

# Tidy both our datasets (make columns names the same, etc.)
tle %<>%
  set_colnames(c(
    "taxon", "specimen", "voucher",
    "ETS", "ITS", "rbcL", "rps16", "trnLF", "n_markers_available"
  )) %>%
  select(-specimen, -n_markers_available)
rvm %<>%
  mutate(taxon = str_replace_all(taxon, "_", " ")) %>%
  select(-n_markers_available)

# Merge datasets
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
  group_by(taxon, voucher) %>%
  # Score markers available as 1, unavailable as 0
  mutate_at(vars(c("ETS", "ITS", "rbcL", "rps16", "trnLF")),
    ~ifelse(is.na(.), 0, .)
  ) %>%
  # Sum up number of markers available per voucher
  mutate(n_markers = ETS + ITS + rps16 + trnLF + rbcL) %>%
  arrange(taxon, voucher)
write_csv(sequences_tidy, "sequences-by-voucher.csv")

# Recalculate n_markers by taxa
sequences_tidy_taxa <- sequences_tidy %>%
  select(-n_markers, -voucher) %>%
  group_by(taxon) %>%
  summarise_at(vars(c("ETS", "ITS", "rbcL", "rps16", "trnLF")),
    ~ifelse(sum(.) > 0, 1, 0)
  ) %>%
  arrange(taxon) %>%
  mutate(n_markers = ETS + ITS + rbcL + rps16 + trnLF)
write_csv(sequences_tidy_taxa, "sequences-by-taxon.csv")

# Check for specimens (vouchers) that have differing dets in the 2 matrices
# (to be continued)
sequences_tidy_det_check <- sequences_tidy %>%
  ungroup() %>%
  arrange(voucher) %>%
  group_by(voucher) %>%
  summarise(n_dets = n(), taxa = paste(taxon, collapse = ", ")) %>%
  arrange(desc(n_dets), voucher) %>%
  filter(n_dets > 1, voucher != "?")
write_csv(sequences_tidy_det_check, "dets-to-check.csv")

# Useful statistics:
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
