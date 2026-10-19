accum_curve <- function(all_species_data, n_perm = 300, seed = 20260101) {
  fig_dir <- here::here("data-raw", "analysis", "carp", "figures")
  txt_dir <- here::here("data-raw", "analysis", "carp", "outputs")

  # Fix the random seed so the bootstrap summary is reproducible run-to-run.
  set.seed(seed)

  SiteID <- Scientific_Name <- sample_detection <- #nolint
    sites <- richness <- lower <- upper <- NULL #nolint

  # Build a site-by-species presence/absence matrix.
  pa <-
    all_species_data |>
    dplyr::select(SiteID, Scientific_Name, sample_detection) |>
    tidyr::pivot_wider(names_from = Scientific_Name,
                       values_from = sample_detection,
                       values_fill = 0,
                       values_fn = max)
  pa$SiteID <- NULL
  pa_mat <- as.matrix(pa)
  n_sites <- nrow(pa_mat)

  # Permutation-based species accumulation (sampling sites without replacement).
  accum_mat <- matrix(NA, nrow = n_perm, ncol = n_sites)
  for (perm in seq_len(n_perm)) {
    site_order <- sample(seq_len(n_sites), n_sites, replace = FALSE)
    species_seen <- rep(FALSE, ncol(pa_mat))
    richness_vec <- numeric(n_sites)
    for (i in seq_len(n_sites)) {
      species_seen <- species_seen | (pa_mat[site_order[i], ] == 1)
      richness_vec[i] <- sum(species_seen)
    }
    accum_mat[perm, ] <- richness_vec
  }

  accum_summary <-
    data.frame(
      sites = seq_len(n_sites),
      richness = colMeans(accum_mat),
      lower = apply(accum_mat, 2, stats::quantile, probs = 0.025),
      upper = apply(accum_mat, 2, stats::quantile, probs = 0.975)
    )

  accum_plot <-
    accum_summary |>
    ggplot2::ggplot(ggplot2::aes(x = sites, y = richness)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = "lightblue", alpha = 0.4) +
    ggplot2::geom_line(color = "blue") +
    ggplot2::theme_bw() +
    ggplot2::xlab("Number of sites (X)") +
    ggplot2::ylab("Accumulated species richness")

  ggplot2::ggsave(filename = here::here(fig_dir, "accum_curve.jpg"),
                  plot = accum_plot,
                  dpi = 300,
                  width = 8, height = 6)

  readr::write_csv(x = accum_summary,
                   file = here::here(txt_dir, "accum_curve_summary.csv"))
  accum_plot
}
