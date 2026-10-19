format_wide_2 <- function(all_coef_data) {
  parameter <- med_prob <- species <- par_plot <- NULL
  all_coef_wide <-
    all_coef_data |>
    dplyr::filter(parameter != "lp__") |>
    dplyr::mutate(par_plot = gsub("beta_|delta_", "", parameter)) |>
    dplyr::mutate(par_plot = ifelse(par_plot == "psi", "θ[1]", "θ[2]")) |>
    dplyr::select(species, par_plot, med_prob) |>
    tidyr::pivot_wider(names_from = par_plot, values_from = med_prob)
  all_coef_wide
}