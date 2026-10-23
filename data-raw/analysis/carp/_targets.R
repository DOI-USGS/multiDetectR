library(targets)
library(tarchetypes)
library(here)

# Input data lives in the separate data-release repository, not in this package
# repo (the raw USFWS eDNA data is sensitive and is removed from the public
# package per the DGEC review). By default we read it from the `data-release`
# repo checked out alongside the package (a sibling of the package root that
# here::here() anchors to). Override with the MULTIDETECTR_DATA environment
# variable to point at the data wherever it is obtained (e.g. a released DOI
# archive). See data-raw/analysis/carp/data/README.md.
release_dir <- Sys.getenv(
  "MULTIDETECTR_DATA",
  unset = file.path(dirname(here::here()), "data-release")
)
data_dir <- file.path(release_dir, "input", "carp")
r_dir <- here::here("data-raw", "analysis", "carp", "R")

all_species_path <- file.path(data_dir, "all_species.csv")

tar_source(files = r_dir)
tar_option_set(packages = c("tibble", "readr", "dplyr", "ggplot2", "here",
  "tidyverse", "rstan", "occStan", "tibble",
  "ggthemes", "scales", "multiDetectR", "tidyr", "forcats"
))

cost_in_2 <- c(200, 20)
cost_in_3 <- c(200, 20, 2)
d_in <- 0.95

list(
  # format input data
  tar_target(name = file, command = all_species_path, format = "file"),
  tar_target(name = all_species_data, command = get_data(file)),
  tar_target(name = plot_input, command = plot_inputs(all_species_data)),
  # fit 2-level model, format and extract outputs
  tar_target(name = all_coef_data_2,
             command = fit_species_stan_2(all_species_data)),
  tar_target(name = plot_coef_2, command = plot_outputs_2(all_coef_data_2)),
  tar_target(name = coef_plot_2, command = plot_theta_hist_2(all_coef_data_2)),
  tar_target(name = all_coef_wide_2, command = format_wide_2(all_coef_data_2)),
  tar_target(name = coef_corr_plot_2,
             command = plot_coef_corr_2(all_coef_wide_2)),
  # allocate resources for 2-level model
  tar_target(name = allocate_rep_2,
             command = represent_theta_2(cost_in_2, d_in)),
  tar_target(name = allocate_emp_2,
             command = empirical_theta_2(all_coef_wide_2, cost_in_2, d_in)),
  # fit 3-level model, format and extract outputs
  tar_target(name = all_coef_data_3,
             command = fit_species_stan_3(all_species_data)),
  tar_target(name = plot_coef_3, command = plot_outputs_3(all_coef_data_3)),
  tar_target(name = coef_plot_3, command = plot_theta_hist_3(all_coef_data_3)),
  tar_target(name = all_coef_wide_3, command = format_wide_3(all_coef_data_3)),
  tar_target(name = coef_corr_plot_3,
             command = plot_coef_corr_3(all_coef_wide_3)),
  # allocate resources for 3-level model
  tar_target(name = allocate_rep_3,
             command = represent_theta_3(cost_in_3, d_in)),
  tar_target(name = allocate_emp_3,
             command = empirical_theta_3(all_coef_wide_3, cost_in_3, d_in)),
  # Larson's method: species accumulation and turnover curves
  # (freshwater carp data cannot support the Pareto front or power-law
  # richness model, which are marine-only per outline_DL5.qmd)
  tar_target(name = accum_curve_plot,
             command = accum_curve(all_species_data)),
  tar_target(name = accum_curve_year_plot,
             command = accum_curve_year(all_species_data)),
  tar_target(name = accum_reps_plot,
             command = accum_reps(all_species_data)),
  tar_target(name = turnover_curve_plot,
             command = turnover_curve(all_species_data))
)