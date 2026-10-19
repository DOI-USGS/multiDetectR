format_wide_3 <- function(all_coef_data) {
  parameter <- med_prob <- species <- par_plot <- #nolint
    `θ[1]` <- `θ[2]` <- `θ[3]` <- NULL #nolint

  # Map the three-level parameters to the nesting levels used by
  # multiDetectR::allocate_assemblage (outermost to innermost):
  #   beta_psi    -> theta[1] (site occupancy, psi)
  #   alpha_theta -> theta[2] (sample availability, theta)
  #   delta_p     -> theta[3] (sub-replicate detection, p)
  all_coef_wide <-
    all_coef_data |>
    dplyr::filter(parameter != "lp__") |>
    dplyr::mutate(par_plot = dplyr::case_when(
      grepl("beta_psi", parameter) ~ "θ[1]",
      grepl("alpha_theta", parameter) ~ "θ[2]",
      grepl("delta_p", parameter) ~ "θ[3]"
    )) |>
    dplyr::select(species, par_plot, med_prob) |>
    tidyr::pivot_wider(names_from = par_plot, values_from = med_prob) |>
    dplyr::select(species, `θ[1]`, `θ[2]`, `θ[3]`)
  all_coef_wide
}
