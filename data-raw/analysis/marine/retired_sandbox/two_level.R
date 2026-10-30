library(tidyverse)
library(occstanhm)
library(ggthemes)

dat_3_may <-
  read_csv("ptconception_occupancy_3level_may.csv", show_col_types = FALSE) |>
  rename(
    sample = PCRReplicate,
    revisit = FieldReplicate,
    unit = Site
  ) |>
  mutate(consensus = Consensus,
         y = ifelse(value > 0, 1, 0)) |>
  separate_wider_delim(Consensus, delim = "/",
                       names_sep = "", too_few = "align_start") |>
  select(unit, revisit, sample, species = Consensus1, consensus, y) |>
  mutate(unit_id = as.integer(factor(unit)),
         species_id = as.integer(factor(species)),
         revisit_id = as.integer(revisit)) |>
  group_by(unit, revisit, species) |>
  mutate(z_obs = ifelse(sum(y) > 0, 1, 0)) |>
  ungroup()

dat_3_may

dat_3_may |>
  filter(z_obs != y)

## Remove species
spp_unit_summary <-
  dat_3_may |>
  group_by(unit, species) |>
  summarize(y_sum = sum(y),
            z_sum = sum(z_obs),
            .groups = "drop") |>
  group_by(species) |>
  summarize(sites_with_spp = sum(y_sum), .groups = "drop")

spp_unit_summary |>
  arrange(sites_with_spp)


spp_unit_summary |>
  arrange(-sites_with_spp) |>
  print(n = 40)

spp_remove <-
  spp_unit_summary |>
  filter(sites_with_spp < 10) |>
  pull(species)

dat_3_may_2c <-
  dat_3_may |>
  filter(!species %in% spp_remove) |>
  arrange(unit, species, revisit) |>
  mutate(unit = unit,
         unit_id = as.integer(unit_id),
         species_id = as.integer(factor(species))) |>
  rowid_to_column("index")

dat_3_may_2c |>
  distinct(species) |>
  pull(species)

many_plot <-
  dat_3_may_2c |>
  group_by(unit, species, revisit) |>
  mutate(revisit = factor(revisit)) |>
  summarize(total_y = sum(y),
            .groups = "drop") |>
  ggplot(aes(x = revisit, y = species, fill = total_y)) +
  geom_tile() +
  scale_fill_gradient("Samples\nwith\nPCR\ndetections",
                      low = "skyblue", high = "black") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.spacing = unit(2, "lines")) +
  facet_wrap(vars(unit),
             nrow = 2) +
  ylab("Species") +
  xlab("Field Replicate (Sample)")
many_plot

## Continue formatting for Stan
revisit_species_data <-
  dat_3_may_2c |>
  ungroup() |>
  group_by(unit, species) |>
  mutate(min_idx = min(index) - 1) |>
  ungroup() |>
  group_by(unit, unit_id, species, species_id, revisit) |>
  summarize(z_obs = ifelse(sum(y) > 0, 1, 0),
            n_samples = n(),
            revisit_dat_start = min(index - min_idx),
            revisit_dat_stop = max(index - min_idx),
            .groups = "drop") |>
  rowid_to_column("index") |>
  ungroup()

revisit_species_data |>
  print(n = 20)

revisit_species_at_unit_data  <-
  revisit_species_data |>
  group_by(unit, species) |>
  summarize(n_revisits = n(),
            revisit_unit_start = min(index),
            revisit_unit_stop = max(index),
            .groups = "drop") |>
  ungroup()

revisit_species_at_unit_data


unit_species_data <-
  dat_3_may_2c |>
  group_by(unit, unit_id, species, species_id) |>
  summarize(psi_obs = mean(z_obs),
            n_samples = n(),
            unit_spp_dat_start = min(index),
            unit_spp_dat_stop = max(index),
            .groups = "drop")

unit_species_data |>
  print(n = 20)

x_psi <- model.matrix(~ species - 1, revisit_species_data)
colnames(x_psi) <- gsub("species", "", colnames(x_psi))
x_psi |> dim()

x_psi_star <-
  dat_3_may_2c |>
  distinct(unit) |>
  model.matrix(object = ~ 1)

colnames(x_psi_star) <- gsub("species", "", colnames(x_psi_star))
x_psi_star |> head()
x_psi_star |> dim()

x_psi_hyper <-
  dat_3_may_2c |>
  group_by(unit, species) |>
  summarize(n = n(), .groups = "drop") |>
  mutate(one = 1)

## check math for beta coefficient
beta_psi_check <- matrix(1,
                         nrow = dat_3_may_2c |> distinct(unit) |> nrow(),
                         ncol = x_psi |> ncol())
beta_psi_check |> dim()

logit_psi_check <- matrix(1,
                          nrow = revisit_species_data |> nrow(),
                          ncol = 1)

logit_psi_check |> dim()
rowSums(x_psi, beta_psi_check[revisit_species_data$unit_id, ]) |> length()

beta_psi_star_prior <-
  matrix(1.0,
         nrow = x_psi_star |> ncol(),
         ncol = x_psi |> ncol())

## check math for prior
x_psi_star_check <-
  matrix(1,
         nrow = x_psi_star |> nrow(),
         ncol = x_psi_star |> ncol())

x_psi_star_check %*% beta_psi_star_prior  |> dim()
beta_psi_check |> dim()


beta_psi_star_check <-
  matrix(1,
         nrow = 1,
         ncol = dat_3_may_2c |> distinct(species) |> nrow())

# lhs
beta_psi_check |> dim()

# rhs elements
x_psi_star |> dim()
beta_psi_star_check |> dim()

# rhs product
x_psi_star %*%
  beta_psi_star_check |>
  dim()

v_p <- model.matrix(~ species - 1, dat_3_may_2c)
colnames(v_p) <- gsub("species", "", colnames(v_p))
v_p |> head()
v_p |> dim()


delta_p_check <- matrix(1,
                        nrow = dat_3_may_2c |> distinct(unit) |> nrow(),
                        ncol = v_p |> ncol())

delta_p_check |> dim()
logit_p_check <- matrix(1,
                        nrow = dat_3_may_2c |> nrow(),
                        ncol = 1)
logit_p_check |> dim()

rowSums(v_p, delta_p_check[dat_3_may_2c$unit_id, ]) |> length()

## hyper parameters
v_p_star <-
  dat_3_may_2c |>
  distinct(unit) |>
  model.matrix(object = ~ 1)

v_p_star |> dim()
v_p_star |> head()

delta_p_star_prior <-
  matrix(1,
         nrow = v_p_star |> ncol(),
         ncol = v_p |> ncol())

v_p_star_check <-
  matrix(1,
         nrow = v_p_star |> nrow(),
         ncol = v_p_star |> ncol())

delta_p_star_check <-
  matrix(1,
         nrow = 1,
         ncol = dat_3_may_2c |> distinct(species) |> nrow())
# lhs
delta_p_check |> dim()

# rhs elements
v_p_star |> dim()
delta_p_check |> dim()

# rhs product
v_p_star %*%
  beta_psi_star_check |>
  dim()


stan_data <-
  list(y = dat_3_may_2c |> pull(y),
       z_obs = revisit_species_data |> pull(z_obs),
       # next are the dimensions for data (including vectors)
       n_unit_species = unit_species_data |> nrow(),
       n_total_revisits = revisit_species_data |> nrow(),
       n_total_samples = dat_3_may_2c |> nrow(),
       n_revisits = revisit_species_at_unit_data |> pull(n_revisits),
       n_samples = revisit_species_data |> pull(n_samples),
       m_beta = x_psi |> ncol(),
       n_beta = dat_3_may_2c |> distinct(unit) |> nrow(),
       m_beta_star = x_psi_star |> ncol(),
       n_beta_star = x_psi_star |> nrow(),
       m_delta = v_p |> ncol(),
       n_delta = dat_3_may_2c |> distinct(unit_id) |> nrow(),
       m_delta_star = v_p_star |> ncol(),
       n_delta_star = v_p_star |> nrow(),
       # indexing
       unit_spp_dat_start = unit_species_data |> pull(unit_spp_dat_start),
       unit_spp_dat_stop = unit_species_data |> pull(unit_spp_dat_stop),
       revisit_dat_start = revisit_species_data |> pull(revisit_dat_start),
       revisit_dat_stop = revisit_species_data |> pull(revisit_dat_stop),
       revisit_unit_start = revisit_species_at_unit_data |>
         pull(revisit_unit_start),
       revisit_unit_stop = revisit_species_at_unit_data |>
         pull(revisit_unit_stop),
       # psi-level inputs
       x_psi = x_psi,
       x_psi_star = x_psi_star,
       # p-level inputs
       v_p = v_p,
       v_p_star = v_p_star,
       # psi-level priors
       beta_psi_star_prior = beta_psi_star_prior,
       beta_psi_star_sd_prior = 0.1,
       eta_psi = 1.0,
       jj_psi = revisit_species_data |> dplyr::pull(unit_id),
       # p-level priors
       delta_p_star_prior = delta_p_star_prior,
       delta_p_star_sd_prior = 0.1,
       eta_p = 1.0,
       jj_p = dat_3_may_2c |> dplyr::pull(unit_id),
       # priors for sigmas
       sigma_psi_prior = 2.5,
       sigma_psi_prior_sd = 0.05,
       sigma_p_prior = 1.75,
       sigma_p_prior_sd = 0.05,
       # setting for reduce_sum
       grainsize = 1,
       ## Store data (not used in Stan)
       dat = dat_3_may_2c,
       unit_species_data = unit_species_data)


n_chains <- 3
n_parallel_chains <- 3
n_threads_per_chain <- 5
n_refresh <- 100

run_stan <- TRUE
if (run_stan) {
  # for sampling numbers
  n_warmup <- 500 # 2000
  n_sample <- 200 # 2000

  fit <- occstanhm_2(stan_data,
                     n_chains = n_chains,
                     n_parallel_chains = n_parallel_chains,
                     n_threads_per_chain = n_threads_per_chain,
                     n_refresh = n_refresh,
                     n_warmup = n_warmup,
                     n_sample = n_sample)

  fit_summary <-
    fit$summary(.cores = 15)
  saveRDS(fit, "fit_2_level.RDS")
  saveRDS(fit_summary, "fit_summary_2_level.RDS")
} else {
  fit <- readRDS("fit_2_level.RDS")
  fit_summary <- readRDS("fit_summary_2_level.RDS")
}

## Look at results
beta_detail <-
  fit_summary |>
  filter(grepl("beta_psi\\[", variable)) |>
  select(variable, q5, median, q95) |>
  mutate(l95 = plogis(q5),
         u95 = plogis(q95),
         prob = plogis(median)) |>
  select(-q5, -q95, -median) |>
  mutate(unit_id =
         as.integer(gsub("beta_psi\\[(\\d+),(\\d+)\\]", "\\1", variable)),
         species_id = as.integer(gsub("beta_psi\\[(\\d+),(\\d+)\\]",
                                      "\\2", variable))) |>
  full_join(stan_data$unit_species_data,
            by = c("species_id", "unit_id")) |>
  select(variable, unit, species, l95, prob, u95, psi_obs) |>
  mutate(unit = gsub("Unit Unit", "Unit", unit))

beta_detail |>
  ggplot(aes(x = `species`,
             y = prob,
             group = `unit`,
             ymin = l95, ymax = u95)) +
  geom_point() +
  geom_linerange() +
  coord_flip() +
  facet_wrap(facets = vars(`unit`), nrow = 2) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  ylim(c(0, 1)) +
  scale_size_continuous("Proportion\nof\nrevisits\nwith\ndetections") +
  ylab(expression("estimated" ~ psi))


dat_rough_p <-
  stan_data$dat |>
  filter(z_obs > 0) |>
  group_by(unit, unit_id, species, species_id) |>
  summarize(p_rough = mean(y), .groups = "drop")

dat_p_all <-
  stan_data$dat |>
  distinct(unit, unit_id, species, species_id) |>
  full_join(dat_rough_p, by = c("unit", "unit_id",
                                "species",
                                "species_id"))
p_details <-
  fit_summary |>
  filter(grepl("delta_p\\[", variable)) |>
  mutate(unit_id = as.integer(gsub("delta_p\\[(\\d+),(\\d+)\\]", "\\1",
                                   variable)),
         species_id = as.integer(gsub("delta_p\\[(\\d+),(\\d+)\\]",
                                      "\\2", variable))) |>
  full_join(dat_p_all, by = c("species_id", "unit_id")) |>
  mutate(prob = plogis(median),
         l95 = plogis(q5),
         u95 = plogis(q95)) |>
  select(variable, unit, species, l95, prob, u95, p_rough) |>
  mutate(diff = p_rough - prob)

p_details |> print(n = 30)

p_details |>
  ggplot(aes(x = `species`,
             y = prob,
             group = `unit`,
             ymin = l95, ymax = u95)) +
  geom_point() +
  geom_linerange() +
  coord_flip() +
  facet_wrap(facets = vars(`unit`), nrow = 2) +
  theme_bw() +
  theme(strip.background = element_blank()) +
  ylim(c(0, 1)) +
  ylab("p")

fit_summary |>
  filter(grepl("^Omega_psi\\[", variable)) |>
  select(variable, q5, median, q95)


fit_summary |>
  filter(grepl("^Omega_p\\[", variable)) |>
  select(variable, q5, median, q95)


gsub_pattern <- "^Omega_(\\w+)\\[(\\d+),(\\d+)\\]"

spp_key <-
  stan_data$unit_species_data |>
  distinct(species, species_id) |>
  mutate(species_row = species,
         species_col = species) |>
  select(-species)

corr_summary <-
  fit_summary |>
  filter(grepl(gsub_pattern, variable)) |>
  mutate(species_row_id = as.integer(gsub(gsub_pattern, "\\2", variable)),
         species_col_id = as.integer(gsub(gsub_pattern, "\\3", variable)),
         level = gsub(gsub_pattern, "\\1", variable)) |>
  full_join(spp_key |> select(species_id, species_row),
            by = c("species_row_id" = "species_id")) |>
  full_join(spp_key |> select(species_id, species_col),
            by = c("species_col_id" = "species_id"))

corr_summary |>
  select(-species_row_id, -species_col_id, -variable) |>
  filter(mean != 1) |>
  arrange(-abs(median)) |>
  print(n = 20) |>
  filter(level == "psi")

corr_summary |>
  select(-species_row_id, -species_col_id, -variable) |>
  filter(mean != 1) |>
  arrange(-abs(median)) |>
  print(n = 20) |>
  filter(level == "p")

fit_summary |>
  filter(grepl(gsub_pattern, variable)) |>
  mutate(species_row_id = as.integer(gsub(gsub_pattern, "\\2", variable)),
         species_col_id = as.integer(gsub(gsub_pattern, "\\3", variable)),
         level = gsub(gsub_pattern, "\\1", variable)) |>
  full_join(spp_key |> select(species_id, species_row),
            by = c("species_row_id" = "species_id")) |>
  full_join(spp_key |> select(species_id, species_col),
            by = c("species_col_id" = "species_id")) |>
  ggplot(aes(x = species_row, y = median, ymin = q5, ymax = q95,
             color = level)) +
  geom_point(position = position_dodge(width = 1.01)) +
  geom_linerange(position = position_dodge(width = 1.01)) +
  coord_flip() +
  facet_grid(cols = vars(species_col)) +
  geom_hline(yintercept = 0, color = "red") +
  theme_bw() +
  theme(strip.background = element_blank()) +
  scale_color_colorblind()

corr_summary |>
  filter(species_row != species_col) |>
  ggplot(aes(x = species_row, y = species_col, fill = median)) +
  geom_tile() +
  scale_fill_gradient2("correlation",
                       low = "red", high = "blue", mid = "white",
                       breaks = seq(-1, 1, by = 0.25),
                       limits = c(-1, 1)) +
  facet_grid(cols = vars(level)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -65, vjust = 0.5, hjust = 0),
        strip.background = element_blank())

corr_summary |>
  filter(species_row != species_col) |>
  select(median, level, species_row, species_col) |>
  pivot_wider(names_from = level, values_from = median) |>
  ggplot(aes(x = psi, y = p)) +
  geom_point()
