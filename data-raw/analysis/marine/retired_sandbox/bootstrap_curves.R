library(tidyverse)
library(ggthemes)
source("helper_files.R")

n_iter <- 300
n_units_loop <- seq(1, 10)
n_samples_loop <- seq(1, 2)
n_k_subsamples <- seq(1, 10)
n_iter_loop <- seq(1, n_iter)

dat_may_3_raw <-
  read_csv("ptconception_occupancy_3level_may.csv") |>
  rename(unit = Site,
         sample = FieldReplicate) |>
  mutate(consensus = Consensus,
         y = ifelse(value > 0, 1, 0)) |>
  separate_wider_delim(Consensus, delim = "/",
                       names_sep = "", too_few = "align_start") |>
  select(unit, sample, PCRReplicate,
         species = Consensus1, consensus, y) |>
  ungroup()

curve_data_3_f <- tibble()
system.time({
  for (unit_idx in n_units_loop) {
    for (sample_idx in n_samples_loop) {
      for (k_idx in n_k_subsamples) {
        for (iter in n_iter_loop) {
          curve_data_3_f <-
            bind_rows(curve_data_3_f,
                      species_from_event(n_units = unit_idx,
                                         n_samples = sample_idx,
                                         k_subsamples = k_idx,
                                         data_in = dat_may_3_raw,
                                         replace = FALSE))
        }
      }
    }
  }
}) |>
  print()

write_csv(file = "curve_data_3_rep_f.csv", x = curve_data_3_f)

curve_data_3_t <- tibble()
system.time({
  for (unit_idx in n_units_loop) {
    for (sample_idx in n_samples_loop) {
      for (k_idx in n_k_subsamples) {
        for (iter in n_iter_loop) {
          curve_data_3_t <-
            bind_rows(curve_data_3_t,
                      species_from_event(n_units = unit_idx,
                                         n_samples = sample_idx,
                                         k_subsamples = k_idx,
                                         data_in = dat_may_3_raw,
                                         replace = TRUE))
        }
      }
    }
  }
}) |>
  print()
write_csv(file = "curve_data_3_rep_t.csv", x = curve_data_3_t)

## Run 2 level
dat_may_2_raw <-
  read_csv("ptconception_occupancy_2level_may.csv") |>
  rename(unit = Site,
         sample = FieldReplicate) |>
  mutate(consensus = Consensus,
         y = ifelse(value > 0, 1, 0)) |>
  separate_wider_delim(Consensus, delim = "/",
                       names_sep = "", too_few = "align_start") |>
  select(unit, sample, PCRReplicate,
         species = Consensus1, consensus, y) |>
  ungroup()

n_units_loop <- seq(1, 10)
n_samples_loop <- seq(1, 2)
n_k_subsamples <- seq(1, 2)
n_iter_loop <- seq(1, n_iter)

curve_data_2_f <- tibble()
system.time({
  for (unit_idx in n_units_loop) {
    for (sample_idx in n_samples_loop) {
      for (k_idx in n_k_subsamples) {
        for (iter in n_iter_loop) {
          curve_data_2_f <-
            bind_rows(curve_data_2_f,
                      species_from_event(n_units = unit_idx,
                                         n_samples = sample_idx,
                                         k_subsamples = k_idx,
                                         data_in = dat_may_2_raw,
                                         replace = FALSE))
        }
      }
    }
  }
}) |>
  print()
write_csv(file = "curve_data_2_rep_f.csv", x = curve_data_2_f)
curve_data_2_t <- tibble()
system.time({
  for (unit_idx in n_units_loop) {
    for (sample_idx in n_samples_loop) {
      for (k_idx in n_k_subsamples) {
        for (iter in n_iter_loop) {
          curve_data_2_t <-
            bind_rows(curve_data_2_t,
                      species_from_event(n_units = unit_idx,
                                         n_samples = sample_idx,
                                         k_subsamples = k_idx,
                                         data_in = dat_may_2_raw,
                                         replace = TRUE))
        }
      }
    }
  }
}) |>
  print()

write_csv(file = "curve_data_2_rep_t.csv", x = curve_data_2_t)
