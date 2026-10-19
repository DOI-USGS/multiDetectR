plot_coef_corr_2 <- function(all_coef_wide) {

  fig_dir <- here::here("data-raw", "analysis", "carp", "figures")

  `θ[1]` <- `θ[2]` <- NULL #nolint
  coef_corr_plot_2 <- 
    all_coef_wide |>
    ggplot2::ggplot(ggplot2::aes(x = `θ[1]`, y = `θ[2]`)) +
    ggplot2::geom_point() +
    ggplot2::xlab(parse(text = "θ[1]")) +
    ggplot2::ylab(parse(text = "θ[2]")) +
    ggplot2::theme_bw() +
    ggplot2::xlim(c(0, 1)) +
    ggplot2::ylim(c(0, 1))

  ggplot2::ggsave(filename = here::here(fig_dir, "coef_corr_plot_2.jpg"),
                  plot = coef_corr_plot_2,
                  dpi = 300,
                  width = 8, height = 12)
  coef_corr_plot_2
}