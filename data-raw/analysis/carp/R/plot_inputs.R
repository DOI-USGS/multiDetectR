plot_inputs <- function(all_species_data) {
  Scientific_Name <- Common_Name <- year <- month <- revisit <- #nolint
  sample_detection <- total <- NULL #nolint

  fig_dir <- here::here("data-raw", "analysis", "carp", "figures")

  detect_summary <-
    all_species_data |>
    dplyr::group_by(Scientific_Name, Common_Name, year, month, revisit) |>
    dplyr::summarise(n = dplyr::n(), total = sum(sample_detection),
                     .groups = "drop")

  detect_plot <-
    detect_summary |>
    dplyr::filter(total != 0) |>
    ggplot2::ggplot(ggplot2::aes(x = factor(month),
                                 y = Scientific_Name, fill = total)) +
    ggplot2::geom_tile() +
    ggplot2::facet_grid(cols = ggplot2::vars(year),
                        scales = "free_x", space = "free") +
    ggplot2::scale_fill_continuous("Total samples\nwith detection",
                                   transform = "pseudo_log",
                                   low = "skyblue", high = "navyblue") +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_blank()) +
    ggplot2::xlab("Month") +
    ggplot2::ylab("Scientific name") +
    ggplot2::scale_y_discrete(limits = rev)
  ggplot2::ggsave(filename = here::here(fig_dir, "raw_data_plot.jpg"),
                  plot = detect_plot,
                  dpi = 300,
                  width = 8, height = 12)
  detect_plot
}