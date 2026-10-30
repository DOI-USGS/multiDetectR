compare_occ_levels_cw <- function(summary_2, summary_3, label) {
  fig_dir <- here::here("data-raw", "analysis", "marine", "figures")
  txt_dir <- here::here("data-raw", "analysis", "marine", "outputs")

  parameter <- species <- se_mean <- se_two <- se_three <- #nolint
    pct_change_se <- parameter_plot <- NULL #nolint

  # Precision gain of the three-level over the two-level occupancy model at
  # Government Point, for one marker. Compares the Monte Carlo standard error of
  # the shared parameters, site occupancy (beta_psi) and detection (delta_p),
  # per species; a negative percent change means the three-level model, which
  # uses the molecular sub-samples, gave a smaller SE (more precise). The old
  # single-event marine Supplement Figure 0 analysis, generalized to any marker.
  shared <- c("beta_psi", "delta_p")

  two_lvl <-
    summary_2 |>
    dplyr::filter(parameter %in% shared) |>
    dplyr::select(species, parameter, se_two = se_mean)

  three_lvl <-
    summary_3 |>
    dplyr::filter(parameter %in% shared) |>
    dplyr::select(species, parameter, se_three = se_mean)

  cmp <-
    dplyr::inner_join(two_lvl, three_lvl,
                      by = c("species", "parameter")) |>
    dplyr::mutate(pct_change_se = 100 * (se_three - se_two) / abs(se_two),
                  parameter_plot = dplyr::recode(parameter,
                                                 "beta_psi" = "psi",
                                                 "delta_p" = "p"))

  compare_plot <-
    cmp |>
    ggplot2::ggplot(ggplot2::aes(x = forcats::fct_rev(species),
                                 y = pct_change_se, fill = parameter_plot)) +
    ggplot2::geom_col(alpha = 0.85) +
    ggplot2::facet_wrap(ggplot2::vars(parameter_plot),
                        labeller = ggplot2::label_parsed) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        color = "gray40") +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      name = "Parameter",
      values = ggthemes::colorblind_pal()(8),
      labels = c(psi = expression(psi), p = expression(p))) +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank()) +
    ggplot2::ggtitle(paste("Marker:", label)) +
    ggplot2::xlab("Marine species") +
    ggplot2::ylab("Percent change in standard error (%)")

  ggplot2::ggsave(
    filename = here::here(fig_dir, paste0("compare_occ_levels_", label, ".jpg")),
    plot = compare_plot,
    dpi = 300,
    width = 8, height = 11)

  readr::write_csv(
    x = cmp,
    file = file.path(txt_dir,
                     paste0("compare_occ_levels_", label, ".csv")))
  compare_plot
}
