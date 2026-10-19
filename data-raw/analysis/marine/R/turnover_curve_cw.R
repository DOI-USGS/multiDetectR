turnover_curve_cw <- function(marine_data, label, n_boot = 300,
                              seed = 20260101) {
  fig_dir <- here::here("data-raw", "analysis", "marine", "figures")
  txt_dir <- here::here("data-raw", "analysis", "marine", "outputs")

  # Fix the random seed so the bootstrap summary is reproducible run-to-run.
  set.seed(seed)

  unit <- species <- present <- #nolint
    sites <- turnover <- lower <- upper <- NULL #nolint

  # Site-by-species presence/absence matrix for one marker (see accum_curve_cw).
  pa <-
    marine_data |>
    dplyr::distinct(unit, species) |>
    dplyr::mutate(present = 1L) |>
    tidyr::pivot_wider(names_from = species, values_from = present,
                       values_fill = 0)
  pa$unit <- NULL
  pa_mat <- as.matrix(pa)
  n_sites <- nrow(pa_mat)

  # Jaccard dissimilarity between the running community and the next site.
  jaccard <- function(a, b) {
    shared <- sum(a == 1 & b == 1)
    only_a <- sum(a == 1 & b == 0)
    only_b <- sum(a == 0 & b == 1)
    denom <- shared + only_a + only_b
    if (denom == 0) return(0)
    (only_a + only_b) / denom
  }

  # Bootstrap the site order (with replacement) and track how much the
  # cumulative community changes as each site is added.
  diss_mat <- matrix(NA, nrow = n_boot, ncol = n_sites)
  for (b in seq_len(n_boot)) {
    sampled <- sample(seq_len(n_sites), n_sites, replace = TRUE)
    cum_comm <- rep(0, ncol(pa_mat))
    dvec <- numeric(n_sites)
    for (i in seq_len(n_sites)) {
      prev <- cum_comm
      cum_comm <- as.numeric((prev + pa_mat[sampled[i], ]) > 0)
      dvec[i] <- if (i == 1) NA else jaccard(prev, cum_comm)
    }
    diss_mat[b, ] <- dvec
  }

  # Turnover is undefined at the first site; drop that NA row.
  turnover_summary <-
    data.frame(
      marker = label,
      sites = seq_len(n_sites),
      turnover = apply(diss_mat, 2, mean, na.rm = TRUE),
      lower = apply(diss_mat, 2, stats::quantile, probs = 0.025,
                    na.rm = TRUE),
      upper = apply(diss_mat, 2, stats::quantile, probs = 0.975,
                    na.rm = TRUE)
    ) |>
    dplyr::filter(sites > 1)

  turnover_plot <-
    turnover_summary |>
    ggplot2::ggplot(ggplot2::aes(x = sites, y = turnover)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = "pink", alpha = 0.3) +
    ggplot2::geom_line(color = "darkred") +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste("Marker:", label)) +
    ggplot2::xlab("Number of sites (X)") +
    ggplot2::ylab("Mean Jaccard turnover")

  ggplot2::ggsave(
    filename = here::here(fig_dir, paste0("turnover_curve_cw_", label, ".jpg")),
    plot = turnover_plot,
    dpi = 300,
    width = 8, height = 6)

  readr::write_csv(
    x = turnover_summary,
    file = file.path(txt_dir, paste0("turnover_curve_cw_", label, ".csv")))
  turnover_plot
}
