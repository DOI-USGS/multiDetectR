plot_theta_hist_2 <- function(all_coef_data){
  fig_dir <- here::here("data-raw", "analysis", "carp", "figures")

  med_prob <- par_plot <- NULL

  coef_two_plot <-
    all_coef_data |>
    dplyr::filter(parameter != "lp__") |>
    dplyr::mutate(par_plot = gsub("beta_|delta_", "", parameter)) |>
    dplyr::mutate(par_plot = ifelse(par_plot == "psi", "θ[1]", "θ[2]")) |>
    ggplot2::ggplot(ggplot2::aes(x = med_prob)) +
    ggplot2::facet_grid(cols = ggplot2::vars(par_plot),
                        labeller = ggplot2::label_parsed) +
    ggplot2::theme_bw() +
    ggplot2::geom_histogram() +
    ggplot2::theme(strip.background = ggplot2::element_blank())

  ggplot2::ggsave(filename = here::here(fig_dir, "coef_two_plot.jpg"),
                  plot = coef_two_plot,
                  dpi = 300,
                  width = 8, height = 12)
  coef_two_plot
}