library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
library(occStan)

## two-level model with two-level data
## Important assumption to ask about
## merge two PCR replicates into one
dat_2_may <-
  read_csv("ptconception_occupancy_2level_may.csv",
           show_col_types = FALSE) |>
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
  summarize(y = max(y), n_agg = n(), k_size = 1,
            .groups = "drop")

dat_2_may_unit <-
  dat_2_may |>
  group_by(unit, species) |>
  summarize(n = n(), z_obs = ifelse(sum(y) > 0, 1, 0))

unit_spp_summary <-
  dat_2_may_unit |>
  group_by(species) |>
  summarize(n_units_per_site = sum(z_obs))

spp_consensus_key <-
  dat_2_may |>
  distinct(species, consensus)

spp_use <-
  unit_spp_summary |>
  filter(n_units_per_site > 0) |>
  pull(species)

model_estimates_2 <- tibble()

system.time({
  for (spp_idx in spp_use) {
    print(spp_idx)
  
    sample_dat_use <-
      dat_2_may |>
      filter(species == spp_idx)
    
    unit_dat_use <-
      dat_2_may_unit |>
      filter(species == spp_idx)
    
    stan_data_2 <-
      format_occ_2(unit_data = unit_dat_use,
                   sample_data = sample_dat_use)
    
    out_2 <- occ_2_stan(stan_data_2, iter = 4000)
  
    out_summary_2 <-
      rstan::summary(out_2, probs = c(0.975, 0.95, 0.5, 0.05, 0.025))$summary
    
    model_estimates_2 <-
      bind_rows(model_estimates_2,
                out_summary_2 |>
                  as.data.frame() |>
                  rownames_to_column("parameter") |>
                  as_tibble() |>
                  filter(!grepl("logit_p|logit_psi", parameter)) |>
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

model_estimates_2 <-
  model_estimates_2 |>
  left_join(spp_consensus_key, by = "species",
            relationship = "many-to-many") |>
  mutate(par_plot = case_when(grepl("psi\\[", parameter) ~ "Ψ",
                              grepl("p\\[", parameter) ~ "p",
                              grepl("theta\\[", parameter) ~ "θ",
                              .default = parameter),
         model = "two-level",
         data = "two-level")

write_csv(file = "model_estimates_2_data_2_model.csv", x = model_estimates_2)



## three-level model with two-level data
dat_3_may <-
  read_csv("ptconception_occupancy_2level_may.csv", show_col_types = FALSE) |>
  rename(unit = Site,
         sample = FieldReplicate) |>
  mutate(consensus = Consensus,
         y = ifelse(value > 0, 1, 0)) |>
  separate_wider_delim(Consensus, delim = "/",
                       names_sep = "", too_few = "align_start") |>
  select(unit, sample, PCRReplicate,
         species = Consensus1, consensus, y) |>
  group_by(unit, sample, species, consensus) |>
  summarize(y = sum(y), k_size = n(), .groups = "drop") |>
  mutate(sample_occ = ifelse(y > 0, 1, 0)) |>
  group_by(unit, species) |>
  mutate(unit_occ = ifelse(sum(y) > 0, 1, 0),
         unit_factor = factor(unit)) |>
  ungroup()

dat_3_may_unit <-
  dat_3_may |>
  group_by(unit, species) |>
  summarize(unit_occ = mean(unit_occ), .groups = "drop") |>
  mutate(unit_factor = factor(unit))


unit_summary <-
  dat_3_may_unit |>
  group_by(species) |>
  summarize(rough_psi = mean(unit_occ),
            total_occ = sum(unit_occ),
            .groups = "drop")


## loop over species
all_spp <-
  unit_summary |>
  filter(total_occ > 0) |>
  pull(species)

spp_conesnsus_key <-
  dat_3_may |>
  distinct(species, consensus)

model_estimates_3 <- tibble()

system.time({
  for (spp_idx in all_spp){
    print(spp_idx)
    unit_data <-
      dat_3_may_unit |>
      filter(species == spp_idx)
    
    sample_data <-
      dat_3_may |>
      filter(species == spp_idx)
    
    stan_data <-
      format_occ_3(unit_data = unit_data,
                   sample_data = sample_data,
                   psi_formula = ~1,
                   theta_formula = ~1,
                   p_formula = ~1)
    
    occ_3_out <- occ_3_stan(stan_data,
                            iter = 4000)
    
    out_summary <-
      rstan::summary(occ_3_out, probs = c(0.975, 0.95, 0.5, 0.05, 0.025))$summary
    
    model_estimates_3 <-
      bind_rows(model_estimates_3,
                out_summary |>
                  as.data.frame() |>
                  rownames_to_column("parameter") |>
                  as_tibble() |>
                  filter(!grepl("logit_p|logit_theta|logit_psi", parameter)) |>
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

model_estimates_3 <-
  model_estimates_3 |>
  full_join(spp_conesnsus_key, by = "species",
            relationship = "many-to-many") |>
  mutate(par_plot = case_when(grepl("psi\\[", parameter) ~ "Ψ",
                              grepl("p\\[", parameter) ~ "p",
                              grepl("theta\\[", parameter) ~ "θ",
                              .default = parameter),
         model = "three-level",
         data = "two-level")


write_csv(file = "model_estimates_2_data_3_model.csv", x = model_estimates_3)


## three-level model with three-level data
dat_3_may <-
  read_csv("ptconception_occupancy_3level_may.csv", show_col_types = FALSE) |>
  rename(unit = Site,
         sample = FieldReplicate) |>
  mutate(consensus = Consensus,
         y = ifelse(value > 0, 1, 0)) |>
  separate_wider_delim(Consensus, delim = "/",
                       names_sep = "", too_few = "align_start") |>
  select(unit, sample, PCRReplicate,
         species = Consensus1, consensus, y) |>
  group_by(unit, sample, species, consensus) |>
  summarize(y = sum(y), k_size = n(), .groups = "drop") |>
  mutate(sample_occ = ifelse(y > 0, 1, 0)) |>
  group_by(unit, species) |>
  mutate(unit_occ = ifelse(sum(y) > 0, 1, 0),
         unit_factor = factor(unit)) |>
  ungroup()

dat_3_may_unit <-
  dat_3_may |>
  group_by(unit, species) |>
  summarize(unit_occ = mean(unit_occ), .groups = "drop") |>
  mutate(unit_factor = factor(unit))


unit_summary <-
  dat_3_may_unit |>
  group_by(species) |>
  summarize(rough_psi = mean(unit_occ),
            total_occ = sum(unit_occ),
            .groups = "drop")


## loop over species
all_spp <-
  unit_summary |>
  filter(total_occ > 0) |>
  pull(species)

spp_conesnsus_key <-
  dat_3_may |>
  distinct(species, consensus)

model_estimates_3 <- tibble()

system.time({
  for (spp_idx in all_spp){
    print(spp_idx)
    unit_data <-
      dat_3_may_unit |>
      filter(species == spp_idx)
  
    sample_data <-
      dat_3_may |>
      filter(species == spp_idx)
  
    stan_data <-
      format_occ_3(unit_data = unit_data,
                   sample_data = sample_data,
                   psi_formula = ~1,
                   theta_formula = ~1,
                   p_formula = ~1)
  
    occ_3_out <- occ_3_stan(stan_data,
                            iter = 4000)
  
    out_summary <-
      rstan::summary(occ_3_out, probs = c(0.975, 0.95, 0.5, 0.05, 0.025))$summary
  
    model_estimates_3 <-
      bind_rows(model_estimates_3,
                out_summary |>
                  as.data.frame() |>
                  rownames_to_column("parameter") |>
                  as_tibble() |>
                  filter(!grepl("logit_p|logit_theta|logit_psi", parameter)) |>
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

model_estimates_3 <-
  model_estimates_3 |>
  full_join(spp_conesnsus_key, by = "species",
            relationship = "many-to-many") |>
  mutate(par_plot = case_when(grepl("psi\\[", parameter) ~ "Ψ",
                              grepl("p\\[", parameter) ~ "p",
                              grepl("theta\\[", parameter) ~ "θ",
                              .default = parameter),
         model = "three-level",
         data = "three-level")


write_csv(file = "model_estimates_3_data_3_model.csv", x = model_estimates_3)
