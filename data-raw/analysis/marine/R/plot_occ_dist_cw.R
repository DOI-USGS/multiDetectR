plot_occ_dist_cw <- function(summary_3, label) {
  fig_dir <- here::here("data-raw", "analysis", "marine", "figures")

  parameter <- med_prob <- par_plot <- NULL #nolint

  # Distribution of the per-species median occupancy estimates for one marker,
  # faceted by level: site occupancy (psi), sample availability (theta), and
  # sub-replicate detection (p). A histogram scales to the hundreds of species
  # in the coastwide markers, where a per-species point-range plot would be
  # unreadable. Uses the three-level Government Point summary.
  occ_dist_plot <-
    summary_3 |>
    dplyr::filter(parameter != "lp__") |>
    dplyr::mutate(par_plot = dplyr::case_when(
      grepl("beta_psi", parameter) ~ "psi",
      grepl("alpha_theta", parameter) ~ "theta",
      grepl("delta_p", parameter) ~ "italic(p)"),
      par_plot = factor(par_plot, levels = c("psi", "theta", "italic(p)"))) |>
    ggplot2::ggplot(ggplot2::aes(x = med_prob)) +
    ggplot2::geom_histogram(bins = 20) +
    ggplot2::facet_grid(cols = ggplot2::vars(par_plot),
                        labeller = ggplot2::label_parsed) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_blank()) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::ggtitle(paste("Marker:", label)) +
    ggplot2::xlab("Median posterior probability") +
    ggplot2::ylab("Number of species")

  ggplot2::ggsave(
    filename = here::here(fig_dir, paste0("occ_dist_", label, ".jpg")),
    plot = occ_dist_plot,
    dpi = 300,
    width = 9, height = 4)
  occ_dist_plot
}
