accum_curve_year <- function(all_species_data, n_boot = 300, seed = 20260101) {
  fig_dir <- here::here("data-raw", "analysis", "carp", "figures")
  txt_dir <- here::here("data-raw", "analysis", "carp", "outputs")

  # Fix the random seed so the bootstrap summary is reproducible run-to-run.
  set.seed(seed)

  SiteID <- Scientific_Name <- sample_detection <- year <- #nolint
    sites <- richness <- lower <- upper <- NULL #nolint

  out_list <- list()

  # Bootstrap a species accumulation curve (sampling sites with replacement)
  # separately for each sampling year.
  for (yr in sort(unique(all_species_data$year))) {
    df_sub <-
      all_species_data |>
      dplyr::filter(year == yr)

    pa <-
      df_sub |>
      dplyr::select(SiteID, Scientific_Name, sample_detection) |>
      tidyr::pivot_wider(names_from = Scientific_Name,
                         values_from = sample_detection,
                         values_fill = 0,
                         values_fn = max)
    pa$SiteID <- NULL
    pa_mat <- as.matrix(pa)
    n_sites <- nrow(pa_mat)
    if (n_sites < 1) next

    accum_mat <- matrix(NA, nrow = n_boot, ncol = n_sites)
    for (b in seq_len(n_boot)) {
      sampled <- sample(seq_len(n_sites), n_sites, replace = TRUE)
      species_seen <- rep(FALSE, ncol(pa_mat))
      rich_vec <- numeric(n_sites)
      for (i in seq_len(n_sites)) {
        species_seen <- species_seen | (pa_mat[sampled[i], ] == 1)
        rich_vec[i] <- sum(species_seen)
      }
      accum_mat[b, ] <- rich_vec
    }

    out_list[[as.character(yr)]] <-
      data.frame(
        year = yr,
        sites = seq_len(n_sites),
        richness = colMeans(accum_mat),
        lower = apply(accum_mat, 2, stats::quantile, probs = 0.025),
        upper = apply(accum_mat, 2, stats::quantile, probs = 0.975)
      )
  }

  accum_year_summary <- dplyr::bind_rows(out_list)

  accum_year_plot <-
    accum_year_summary |>
    ggplot2::ggplot(ggplot2::aes(x = sites, y = richness,
                                 color = factor(year))) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper,
                                      fill = factor(year)),
                         alpha = 0.2, color = NA) +
    ggplot2::geom_line() +
    ggplot2::theme_bw() +
    ggplot2::theme(strip.background = ggplot2::element_blank()) +
    ggplot2::scale_color_manual(name = "Year",
                                values = ggthemes::colorblind_pal()(8)) +
    ggplot2::scale_fill_manual(name = "Year",
                               values = ggthemes::colorblind_pal()(8)) +
    ggplot2::xlab("Number of sites (X)") +
    ggplot2::ylab("Accumulated species richness")

  ggplot2::ggsave(filename = here::here(fig_dir, "accum_curve_year.jpg"),
                  plot = accum_year_plot,
                  dpi = 300,
                  width = 8, height = 6)

  readr::write_csv(x = accum_year_summary,
                   file = here::here(txt_dir, "accum_curve_year_summary.csv"))
  accum_year_plot
}
