plot_theta_hist_3 <- function(all_coef_data){
  fig_dir <- here::here("data-raw", "analysis", "carp", "figures")

  med_prob <- par_plot <- parameter <- NULL

  coef_three_plot <-
    all_coef_data |>
    dplyr::filter(parameter != "lp__") |>
    dplyr::mutate(par_plot = dplyr::case_when(
      grepl("beta_psi", parameter) ~ "θ[1]",
      grepl("alpha_theta", parameter) ~ "θ[2]",
      grepl("delta_p", parameter) ~ "θ[3]"
    )) |>
    ggplot2::ggplot(ggplot2::aes(x = med_prob)) +
    ggplot2::facet_grid(cols = ggplot2::vars(par_plot),
                        labeller = ggplot2::label_parsed) +
    ggplot2::theme_bw() +
    ggplot2::geom_histogram() +
    ggplot2::theme(strip.background = ggplot2::element_blank())

  ggplot2::ggsave(filename = here::here(fig_dir, "coef_three_plot.jpg"),
                  plot = coef_three_plot,
                  dpi = 300,
                  width = 8, height = 12)
  coef_three_plot
}
