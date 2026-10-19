represent_theta_3 <- function(cost_in, d_in) {
  txt_dir <- here::here("data-raw", "analysis", "carp", "outputs")

  # Representative three-level detectabilities (site, sample, replicate)
  # spanning few/medium/many occupied-site scenarios.
  theta <-
    rbind(
      few_site_1 = c(0.15, 0.5, 0.5),
      few_site_2 = c(0.15, 0.15, 0.5),
      medium_site_1 = c(0.4, 0.15, 0.5),
      medium_site_2 = c(0.4, 0.2, 0.5),
      many_site_1 = c(0.75, 0.15, 0.5),
      many_site_2 = c(0.75, 0.75, 0.5)
    )

  allocate_out <-
    multiDetectR::allocate_assemblage(theta = theta,
                                      cost = cost_in, D_target = d_in)

  sink(here::here(txt_dir, "allocate_represent_3.txt"))
  print(allocate_out)
  sink()
  allocate_out
}
