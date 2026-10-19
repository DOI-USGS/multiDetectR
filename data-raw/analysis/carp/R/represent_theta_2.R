represent_theta_2 <- function(cost_in, d_in) {
  txt_dir <- here::here("data-raw", "analysis", "carp", "outputs")

  theta <-
    rbind(
      few_site_1 = c(0.15, 0.5),
      few_site_2 = c(0.15, 0.15),
      medium_site_1 = c(0.4, 0.15),
      medium_site_1 = c(0.4, 0.2),
      many_site_1 = c(0.75, 0.15),
      many_site_1 = c(0.75, 0.75)
    )

  allocate_out <-
    multiDetectR::allocate_assemblage(theta = theta,
                                      cost = cost_in, D_target = d_in)

  sink(here::here(txt_dir, "allocate_represent_2.txt"))
  print(allocate_out)
  sink()
  allocate_out
}