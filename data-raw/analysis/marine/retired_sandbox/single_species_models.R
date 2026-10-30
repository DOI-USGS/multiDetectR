## @knitr load_libraries
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
library(occStan)

## @knitr load_data
dat_may_2 <-
  read_csv("ptconception_occupancy_2level_may.csv") |>
  rename(unit = Site,
         sample = FieldReplicate) |>
  mutate(consensus = Consensus,
         y = ifelse(value > 0, 1, 0)) |>
  separate_wider_delim(Consensus, delim = "/",
                       names_sep = "", too_few = "align_start") |>
  select(unit, sample, PCRReplicate,
         species = Consensus1, consensus, y) |>
  ungroup() |>
  group_by(unit, sample, species, consensus) |>
  summarize(y = sum(y), n_agg = n(), k_size = n(),
            .groups = "drop") |>
  mutate(sample_occ = ifelse(y > 0, 1, 0))

## zero entries have shown up in the past
dat_may_2 |>
  group_by(species) |>
  summarize(total_so = sum(sample_occ)) |>
  filter(total_so == 0) |>
  nrow()

dat_may_2 |>
  group_by(species) |>
  summarize(total_so = sum(sample_occ)) |>
  filter(total_so != 0) |>
  nrow()

dat_may_2 |>
  group_by(species) |>
  summarize(total_so = sum(sample_occ)) |>
  filter(total_so != 0) |>
  nrow()

dat_may_2_unit <-
  dat_may_2 |>
  group_by(unit, species) |>
  summarize(n = n(), z_obs = ifelse(sum(y) > 0, 1, 0),
            .groups = "drop")

dat_may_3 <-
  read_csv("ptconception_occupancy_3level_may.csv") |>
  rename(unit = Site,
         sample = FieldReplicate) |>
  mutate(consensus = Consensus,
         y = ifelse(value > 0, 1, 0)) |>
  separate_wider_delim(Consensus, delim = "/",
                       names_sep = "", too_few = "align_start") |>
  select(unit, sample, PCRReplicate,
         species = Consensus1, consensus, y) |>
  ungroup() |>
  group_by(unit, sample, species, consensus) |>
  summarize(y = sum(y), n_agg = n(), k_size = n(),
            .groups = "drop") |>
  mutate(sample_occ = ifelse(y > 0, 1, 0))

dat_may_3_unit <-
  dat_may_3 |>
  group_by(unit, species) |>
  summarize(n = n(), z_obs = ifelse(sum(y) > 0, 1, 0),
            .groups = "drop")

dat_may_2 |>
  group_by(sample) |>
  summarize(max_k = max(k_size),
            max_y = max(y))

dat_may_3 |>
  group_by(sample) |>
  summarize(max_k = max(k_size),
            max_y = max(y))

## @knitr format_data_2
unit_spp_summary_2_may <-
  dat_may_2_unit |>
  group_by(species) |>
  summarize(n_units_per_site = sum(z_obs))

spp_consensus_key_2_may <-
  dat_may_2 |>
  distinct(species, consensus)

spp_use_2_may <-
  unit_spp_summary_2_may |>
  filter(n_units_per_site > 0) |>
  pull(species)

model_estimates_2_may <- tibble()

system.time({
  for (spp_idx in spp_use_2_may) {
    print(spp_idx)

    sample_dat_use <-
      dat_may_2 |>
      filter(species == spp_idx)

    unit_dat_use <-
      dat_may_2_unit |>
      filter(species == spp_idx)

    stan_data_3 <-
      format_occ_3(unit_data = unit_dat_use,
                   sample_data = sample_dat_use,
                   psi_formula = ~ 1,
                   theta_formula = ~ 1,
                   p_formula = ~ 1)

    out_3 <- occ_3_stan(stan_data_3, iter = 4000)

    out_summary_3 <-
      rstan::summary(out_3, probs = c(0.975, 0.95, 0.5, 0.05, 0.025))$summary

    filter_out <- "logit_|log_|log1m_|target_"
    model_estimates_2_may <-
      bind_rows(
        model_estimates_2_may,
        out_summary_3 |>
          as.data.frame() |>
          rownames_to_column("parameter") |>
          as_tibble() |>
          filter(!grepl(filter_out, parameter)) |>
          mutate(species = spp_idx,
                 prob = plogis(mean),
                 u95 = plogis(`97.5%`),
                 u90 = plogis(`95%`),
                 median = plogis(`50%`),
                 l90 = plogis(`5%`),
                 l95 = plogis(`2.5%`))
      )
  }
})

model_estimates_2_may <-
  model_estimates_2_may |>
  left_join(spp_consensus_key_2_may, by = "species",
            relationship = "many-to-many") |>
  mutate(par_plot = case_when(grepl("psi\\[", parameter) ~ "Ψ",
                              grepl("p\\[", parameter) ~ "p",
                              grepl("theta\\[", parameter) ~ "θ",
                              .default = parameter),
         model = "three-level",
         data = "two-level")

write_csv(file = "model_estimates_2_may.csv", x = model_estimates_2_may)

## @knitr format_data_3
unit_spp_summary_3_may <-
  dat_may_3_unit |>
  group_by(species) |>
  summarize(n_units_per_site = sum(z_obs))

spp_consensus_key_3_may <-
  dat_may_3 |>
  distinct(species, consensus)

spp_use_3_may <-
  unit_spp_summary_3_may |>
  filter(n_units_per_site > 0) |>
  pull(species)

model_estimates_3_may <- tibble()

system.time({
  for (spp_idx in spp_use_3_may) {
    print(spp_idx)

    sample_dat_use <-
      dat_may_3 |>
      filter(species == spp_idx)

    unit_dat_use <-
      dat_may_3_unit |>
      filter(species == spp_idx)

    stan_data_3 <-
      format_occ_3(unit_data = unit_dat_use,
                   sample_data = sample_dat_use,
                   psi_formula = ~ 1,
                   theta_formula = ~ 1,
                   p_formula = ~ 1)

    out_3 <- occ_3_stan(stan_data_3, iter = 4000)

    out_summary_3 <-
      rstan::summary(out_3, probs = c(0.975, 0.95, 0.5, 0.05, 0.025))$summary

    filter_out <- "logit_|log_|log1m_|target_"
    model_estimates_3_may <-
      bind_rows(
        model_estimates_3_may,
        out_summary_3 |>
          as.data.frame() |>
          rownames_to_column("parameter") |>
          as_tibble() |>
          filter(!grepl(filter_out, parameter)) |>
          mutate(species = spp_idx,
                 prob = plogis(mean),
                 u95 = plogis(`97.5%`),
                 u90 = plogis(`95%`),
                 median = plogis(`50%`),
                 l90 = plogis(`5%`),
                 l95 = plogis(`2.5%`))
      )
  }
})

model_estimates_3_may <-
  model_estimates_3_may |>
  left_join(spp_consensus_key_3_may, by = "species",
            relationship = "many-to-many") |>
  mutate(par_plot = case_when(grepl("psi\\[", parameter) ~ "Ψ",
                              grepl("p\\[", parameter) ~ "p",
                              grepl("theta\\[", parameter) ~ "θ",
                              .default = parameter),
         model = "three-level",
         data = "three-level")

write_csv(file = "model_estimates_3_may.csv", x = model_estimates_3_may)
