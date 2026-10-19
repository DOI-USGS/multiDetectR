empirical_theta_3 <- function(all_coef_wide, cost_in, d_in) {

  txt_dir <- here::here("data-raw", "analysis", "carp", "outputs")

  species <- NULL
  theta <-
    all_coef_wide |>
    dplyr::select(-species) |>
    as.matrix()

  rownames(theta) <-
    all_coef_wide |>
    dplyr::pull(species)

  colnames(theta) <- NULL
  allocate_out <-
    multiDetectR::allocate_assemblage(theta = theta,
                                      cost = cost_in, D_target = d_in)

  sink(here::here(txt_dir, "allocate_empirical_3.txt"))
  print(allocate_out)
  sink()
  allocate_out
}
