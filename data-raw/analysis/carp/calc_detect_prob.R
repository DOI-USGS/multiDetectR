
data_dir_out <- here::here("data-raw", "analysis", "carp", "data")
all_coef_data <- readr::read_csv(file = file.path(data_dir_out, "all_summary.csv"))



all_coef_wide <-
  all_coef_data |>
  dplyr::filter(parameter != "lp__") |>
  dplyr::mutate(par_plot = gsub("beta_|delta_", "", parameter)) |>
  dplyr::mutate(par_plot = ifelse(par_plot == "psi", "θ[1]", "θ[2]")) |>
  dplyr::select(species, par_plot, med_prob) |>
  tidyr::pivot_wider(names_from = par_plot, values_from = med_prob)


