## assumption 1, there are zeros in data, not detection only
## assumption 2, ESVId can be ignored for this analysis

library(tidyverse)
library(occstanhm)
library(ggthemes)
library(corrplot)

dat_3_may <-
  read_csv("ptconception_occupancy_3level_may.csv", show_col_types = FALSE)

dat_3_may |>
  distinct(Site)

dat_3_all <-
  dat_3_may |>
  rename(unit = Site,
         sample = FieldReplicate) |>
  mutate(consensus = Consensus) |>
  separate_wider_delim(Consensus, delim = "/",
                       names_sep = "", too_few = "align_start") |>
  mutate(y = ifelse(value > 0, 1, 0),
         species = Consensus1) |>
  select(unit, sample, PCRReplicate,
         species, consensus, y)

dat_3_all |>
  distinct(unit, sample, PCRReplicate, species, consensus, y)

## Site is unit
## FieldReplicate is sample
## Consensus is species
## Each row is a molecular replicate that needs to be summed.
## There is only one revisit


dat_3_all_a <-
  dat_3_all |>
  group_by(unit, sample, species, consensus) |>
  summarize(k_subsample = n(), y = sum(y), .groups = "drop") |>
  group_by(unit, sample, species) |>
  mutate(a_obs = ifelse(sum(y) > 0, 1, 0),
         revisit = 1,
         unit_id = factor(unit)) |>
  ungroup() |>
  group_by(unit, species) |>
  mutate(z_obs = ifelse(sum(a_obs) > 0, 1, 0)) |>
  ungroup() |>
  mutate(unit_id = as.numeric(unit_id))

dat_3_all_a |>
  distinct(unit, unit_id)

## Remove species with few detections
spp_unit_summary <-
  dat_3_all_a |>
  group_by(unit, species) |>
  summarize(y_sum = sum(y),
            a_sum = sum(a_obs),
            .groups = "drop")

spp_unit_summary |>
  arrange(a_sum)

spp_unit_summary |>
  arrange(-a_sum) |>
  print(n = 40)

species_at_sites <-
  spp_unit_summary |>
  group_by(species) |>
  summarize(sites_with_spp = sum(a_sum), n = n())

species_at_sites |>
  arrange(sites_with_spp)

species_at_sites |>
  arrange(-sites_with_spp) |>
  print(n = 20)

species_at_sites |>
  ggplot(aes(x = sites_with_spp)) +
  geom_histogram()

spp_remove <-
  species_at_sites |>
  filter(sites_with_spp < 10) |>
  pull(species)

dat_3_all_a |>
  filter(!species %in% spp_remove)

dat_3_all_b <-
  dat_3_all_a |>
  filter(!species %in% spp_remove) |>
  arrange(unit, species, revisit) |>
  mutate(unit = unit,
         unit_id = as.integer(unit_id),
         species_id = as.integer(factor(species))) |>
  rowid_to_column("index")

dat_3_all_b |>
  distinct(species) |>
  pull(species)

many_plot <-
  dat_3_all_b |>
  group_by(unit, species, sample) |>
  summarize(total_y = sum(y), .groups = "drop") |>
  mutate(sample = factor(sample)) |>
  ggplot(aes(x = sample, y = species, fill = total_y)) +
  geom_tile() +
  scale_fill_gradient("Samples\nwith\ndetections",
                      low = "white", high = "navyblue",
                      breaks = seq(0, 10, by = 2)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.spacing = unit(2, "lines")) +
  facet_wrap(vars(unit),
             nrow = 2) +
  ylab("Species") +
  xlab("Field Replicate (Sample)")
many_plot

dat_3_all_corr <-
  dat_3_all_b |>
  group_by(unit, sample, species) |>
  summarize(a_obs = mean(a_obs, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = species, values_from = a_obs) |>
  select(-unit, -sample) |>
  cor()

corrplot(dat_3_all_corr, tl.col = "blue", bg = "White", tl.srt = 35,
         addCoef.col = "black", type = "lower")

## Aggregate to revisit species data
revisit_species_data <-
  dat_3_all_b |>
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
  ungroup() |>
  mutate(jj_p = unit_id) ## Need to figure out

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
  dat_3_all_b |>
  group_by(unit, unit_id, species, species_id) |>
  summarize(psi_obs = mean(z_obs),
            n_samples = n(),
            unit_spp_dat_start = min(index),
            unit_spp_dat_stop = max(index),
            .groups = "drop")
## species - 1
x_psi <- model.matrix(~ species - 1, revisit_species_data)
colnames(x_psi) <- gsub("species", "", colnames(x_psi))
x_psi |> dim()
x_psi |> head()
x_psi_star <-
  dat_3_all_b |>
  distinct(unit) |>
  model.matrix(object = ~ 1)
colnames(x_psi_star) <- gsub("species", "", colnames(x_psi_star))
beta_psi_star_prior <-
  matrix(1.0,
         nrow = x_psi_star |> ncol(),
         ncol = x_psi |> ncol())
# theta-level
w_theta <- model.matrix(~ species - 1, dat_3_all_b)
colnames(w_theta) <- gsub("species", "", colnames(w_theta))
## hyper parameters
w_theta_star <-
  dat_3_all_b |>
  distinct(unit) |>
  model.matrix(object = ~ 1)

alpha_theta_star_prior <-
  matrix(1,
         nrow = w_theta_star |> ncol(),
         ncol = w_theta |> ncol())
# p-level
v_p <- model.matrix(~ species - 1, dat_3_all_b)
colnames(v_p) <- gsub("species", "", colnames(v_p))
## check math for delta coefficient
delta_p_check <- matrix(1,
                        nrow = dat_3_all_b |> distinct(unit) |> nrow(),
                        ncol = v_p |> ncol())
## hyper parameters
v_p_star <-
  dat_3_all_b |>
  distinct(unit) |>
  model.matrix(object = ~ 1)

delta_p_star_prior <-
  matrix(1,
         nrow = v_p_star |> ncol(),
         ncol = v_p |> ncol())
### Stan data
stan_data <-
  # first are observed data
  list(
    y = dat_3_all_b |> pull(y),
    a_obs = dat_3_all_b |> pull(a_obs),
    k_subsamples = dat_3_all_b |> pull(k_subsample),
    z_obs = revisit_species_data |> pull(z_obs),

    # next are the dimensions for data (including vectors)
    n_unit_species = unit_species_data |> nrow(),
    n_total_revisits = revisit_species_data |> nrow(),
    n_total_samples = dat_3_all_b |> nrow(),

    n_revisits = revisit_species_at_unit_data |> pull(n_revisits),
    n_samples = revisit_species_data |> pull(n_samples),

    m_beta = x_psi |> ncol(),
    n_beta = dat_3_all_b |> distinct(unit) |> nrow(),

    m_beta_star = x_psi_star |> ncol(),
    n_beta_star = x_psi_star |> nrow(),

    m_alpha = w_theta |> ncol(),
    n_alpha = dat_3_all_b |> distinct(unit_id) |> nrow(),

    m_alpha_star = w_theta_star |> ncol(),
    n_alpha_star = w_theta_star |> nrow(),

    m_delta = v_p |> ncol(),
    n_delta = dat_3_all_b |> distinct(unit_id) |> nrow(),

    m_delta_star = v_p_star |> ncol(),
    n_delta_star = v_p_star |> nrow(),

    # indexing
    unit_spp_dat_start = unit_species_data |> pull(unit_spp_dat_start),
    unit_spp_dat_stop = unit_species_data |> pull(unit_spp_dat_stop),

    revisit_dat_start = revisit_species_data |> pull(revisit_dat_start),
    revisit_dat_stop = revisit_species_data |> pull(revisit_dat_stop),

    revisit_unit_start =
    revisit_species_at_unit_data |> pull(revisit_unit_start),
    revisit_unit_stop =
    revisit_species_at_unit_data |> pull(revisit_unit_stop),

    # psi-level inputs
    x_psi = x_psi,
    x_psi_star = x_psi_star,

    # psi-level inputs
    w_theta = w_theta,
    w_theta_star = w_theta_star,

    # p-level inputs
    v_p = v_p,
    v_p_star = v_p_star,

    # psi-level priors
    beta_psi_star_prior = beta_psi_star_prior,
    beta_psi_star_sd_prior = 0.1,
    eta_psi = 1.0,
    jj_psi = revisit_species_data |> dplyr::pull(unit_id),

    # theta-level priors
    alpha_theta_star_prior = alpha_theta_star_prior,
    alpha_theta_star_sd_prior = 0.1,
    eta_theta = 1.0,
    jj_theta = dat_3_all_b |> dplyr::pull(unit_id),

    # p-level priors
    delta_p_star_prior = delta_p_star_prior,
    delta_p_star_sd_prior = 0.1,
    eta_p = 1.0,
    jj_p = dat_3_all_b |> dplyr::pull(unit_id),

    # priors for sigma values
    sigma_psi_prior = 2.5,
    sigma_psi_prior_sd = 0.05,

    sigma_theta_prior = 1.75,
    sigma_theta_prior_sd = 0.05,

    sigma_p_prior = 1.75,
    sigma_p_prior_sd = 0.05,

    # setting for reduce_sum
    grainsize = 1,

    ## Store data (not used in Stan)
    dat = dat_3_all_b,
    unit_species_data = unit_species_data
  )

run_stan <- TRUE
if (run_stan) {
  n_chains <- 3
  n_parallel_chains <- 3
  n_threads_per_chain <- 5
  n_refresh <- 100

  # for sampling numbers
  n_warmup <- 500 # 10000
  n_sample <- 200 # 10000

  fit <- occstanhm_3(stan_data,
                     n_chains = n_chains,
                     n_parallel_chains = n_parallel_chains,
                     n_threads_per_chain = n_threads_per_chain,
                     n_refresh = n_refresh,
                     n_warmup = n_warmup,
                     n_sample = n_sample)
  fit_summary <-
    fit$summary(.cores = 15)
  fit_summary

  saveRDS(fit, "fit_3_level.RDS")
  saveRDS(fit_summary, "fit_summary_3_level.RDS")
} else {
  fit <- readRDS("fit_3_level.RDS")
  fit_summary <- readRDS("fit_summary_3_level.RDS")
}

## look at outputs
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
  select(variable, unit, unit_id, species, species_id, l95, prob, u95, psi_obs)

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


## theta details
theta_detail <-
  fit_summary |>
  filter(grepl("alpha_theta\\[", variable)) |>
  select(variable, q5, median, q95) |>
  mutate(l95 = plogis(q5),
         u95 = plogis(q95),
         prob = plogis(median)) |>
  select(-q5, -q95, -median) |>
  mutate(unit_id =
         as.integer(gsub("alpha_theta\\[(\\d+),(\\d+)\\]", "\\1", variable)),
         species_id = as.integer(gsub("alpha_theta\\[(\\d+),(\\d+)\\]",
                                      "\\2", variable))) |>
  full_join(stan_data$unit_species_data,
            by = c("species_id", "unit_id")) |>
  select(variable, unit, unit_id, species, species_id, l95, prob, u95, psi_obs)

theta_detail |>
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
  ylab(expression("estimated" ~ theta))

dat_rough_p <-
  stan_data$dat |>
  filter(z_obs > 0) |>
  group_by(unit, unit_id, species, species_id) |>
  summarize(p_rough = mean(y / k_subsample), .groups = "drop")

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
  select(variable, unit, unit_id, species, l95, prob, u95, p_rough) |>
  mutate(p_rough = ifelse(is.na(p_rough), 0, p_rough),
         diff = p_rough - prob)

p_details

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

#####
fit_summary |>
  filter(grepl("^Omega_psi\\[", variable)) |>
  select(variable, q5, median, q95)

fit_summary |>
  filter(grepl("^Omega_theta\\[", variable)) |>
  select(variable, q5, median, q95)

fit_summary |>
  filter(grepl("^Omega_p\\[", variable)) |>
  select(variable, q5, median, q95)

###
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
            by = c("species_col_id" = "species_id")) |>
  mutate(level_plot = factor(level,
                             levels = c("psi", "theta", "p"),
                             labels = c("\u03A8",
                                        "\u0398",
                                        "p")))
corr_plot_ci_2 <-
  corr_summary |>
  ggplot(aes(x = species_row, y = median, ymin = q5, ymax = q95,
             color = level_plot, group = level_plot)) +
  geom_point(position = position_dodge(width = 1.01)) +
  geom_linerange(position = position_dodge(width = 1.01)) +
  coord_flip() +
  facet_wrap(vars(species_col)) +
  geom_hline(yintercept = 0, color = "red") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.spacing = unit(1.15, "lines")) +
  xlab("Species") +
  ylab("Estimated correlation") +
  scale_color_colorblind("Model\nlevel")
print(corr_plot_ci_2)

corr_summary |>
  mutate(median_plot = round(median, 2)) |>
  filter(species_row != species_col) |>
  ggplot(aes(x = species_row, y = species_col, fill = median,
             label = median_plot)) +
  geom_tile() +
  scale_fill_gradient2("correlation",
                       low = "red", high = "blue", mid = "white",
                       breaks = seq(-1, 1, by = 0.25),
                       limits = c(-1, 1)) +
  geom_text() +
  facet_grid(cols = vars(level)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = -65, vjust = 0.5, hjust = 0),
        strip.background = element_blank())


corr_summary |>
  filter(species_row != species_col,
         level == "psi") |>
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
  filter(species_row != species_col,
         level == "psi") |>
  arrange(-abs(median))

bayesplot::mcmc_trace(fit$draws(),
                      pars = c("lp__"))
