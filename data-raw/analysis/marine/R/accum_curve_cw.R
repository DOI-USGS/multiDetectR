accum_curve_cw <- function(marine_data, label, n_perm = 300,
                           seed = 20260101) {
  fig_dir <- here::here("data-raw", "analysis", "marine", "figures")
  txt_dir <- here::here("data-raw", "analysis", "marine", "outputs")

  # Fix the random seed so the bootstrap summary is reproducible run-to-run.
  set.seed(seed)

  unit <- species <- present <- #nolint
    sites <- richness <- lower <- upper <- NULL #nolint

  # Coastwide species accumulation across sites for one marker. The reader
  # output lists detections, so a distinct (unit, species) pair is a presence;
  # pivot to a site-by-species presence/absence matrix (every sampled site has
  # detections, so all appear).
  pa <-
    marine_data |>
    dplyr::distinct(unit, species) |>
    dplyr::mutate(present = 1L) |>
    tidyr::pivot_wider(names_from = species, values_from = present,
                       values_fill = 0)
  pa$unit <- NULL
  pa_mat <- as.matrix(pa)
  n_sites <- nrow(pa_mat)

  # Permutation-based accumulation (sampling sites without replacement).
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
      marker = label,
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
    ggplot2::ggtitle(paste("Marker:", label)) +
    ggplot2::xlab("Number of sites (X)") +
    ggplot2::ylab("Accumulated species richness")

  ggplot2::ggsave(
    filename = here::here(fig_dir, paste0("accum_curve_cw_", label, ".jpg")),
    plot = accum_plot,
    dpi = 300,
    width = 8, height = 6)

  readr::write_csv(
    x = accum_summary,
    file = file.path(txt_dir, paste0("accum_curve_cw_", label, ".csv")))
  accum_plot
}
