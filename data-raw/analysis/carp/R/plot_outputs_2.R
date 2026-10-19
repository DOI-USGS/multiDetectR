plot_outputs_2 <- function(all_coef_data) {

  fig_dir <- here::here("data-raw", "analysis", "carp", "figures")

  parameter <- par_plot <- med_prob <- l95 <- u95 <- l90 <- u90 <-
  species <- NULL
  pd <- position_dodge(0.4)

  coef_plot <-
    all_coef_data |>
    dplyr::mutate(par_plot = gsub("beta_|delta_", "", parameter)) |>
    dplyr::mutate(par_plot = ifelse(par_plot == "psi", "θ[1]", "θ[2]")) |>
    dplyr::filter(parameter != "lp__") |>
    ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(species),
                                 y = med_prob, color = par_plot)) +
    ggplot2::geom_point(size = 2, position = pd) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = l95, ymax = u95),
                            linewidth = 1.5,
                            position = pd) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = l90, ymax = u90),
                            linewidth = 1.0,
                            position = pd) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::xlab("Probability") +
    ggplot2::ylab("Species") +
    ggplot2::theme(strip.background = ggplot2::element_blank()) +
    ggthemes::scale_color_colorblind(labels = scales::label_parse()) +
    ggplot2::labs(color = "Parameter")
  ggplot2::ggsave(filename = here::here(fig_dir, "coef_plot.jpg"), dpi = 300,
                  width = 8, height = 12)
  coef_plot
}
