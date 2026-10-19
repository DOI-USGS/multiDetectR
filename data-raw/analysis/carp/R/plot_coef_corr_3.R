plot_coef_corr_3 <- function(all_coef_wide) {

  fig_dir <- here::here("data-raw", "analysis", "carp", "figures")

  `θ[1]` <- `θ[2]` <- `θ[3]` <- NULL #nolint

  # Three levels give three pairwise relationships. Show them side by side:
  #   psi vs theta, psi vs p, and theta vs p.
  pair_data <-
    dplyr::bind_rows(
      all_coef_wide |>
        dplyr::transmute(x = `θ[1]`, y = `θ[2]`,
                         pair = "θ[1] vs θ[2]"),
      all_coef_wide |>
        dplyr::transmute(x = `θ[1]`, y = `θ[3]`,
                         pair = "θ[1] vs θ[3]"),
      all_coef_wide |>
        dplyr::transmute(x = `θ[2]`, y = `θ[3]`,
                         pair = "θ[2] vs θ[3]")
    )

  x <- y <- pair <- NULL
  coef_corr_plot_3 <-
    pair_data |>
    ggplot2::ggplot(ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(ggplot2::vars(pair)) +
    ggplot2::xlab("Probability (outer level)") +
    ggplot2::ylab("Probability (inner level)") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_blank()) +
    ggplot2::xlim(c(0, 1)) +
    ggplot2::ylim(c(0, 1))

  ggplot2::ggsave(filename = here::here(fig_dir, "coef_corr_plot_3.jpg"),
                  plot = coef_corr_plot_3,
                  dpi = 300,
                  width = 12, height = 6)
  coef_corr_plot_3
}
