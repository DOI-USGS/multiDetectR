accum_reps <- function(all_species_data, n_boot = 300, seed = 20260101) {
  fig_dir <- here::here("data-raw", "analysis", "carp", "figures")
  txt_dir <- here::here("data-raw", "analysis", "carp", "outputs")

  # Fix the random seed so the bootstrap summary is reproducible run-to-run.
  set.seed(seed)

  SiteID <- Sample_Barcode_ID <- Scientific_Name <- sample_detection <- #nolint
    n_reps <- replicates <- richness <- lower <- upper <- NULL #nolint

  # Keep only sites that have all three PCR technical replicates (K = 3),
  # so the richness-versus-K curve is built on complete sub-sampling.
  complete_sites <-
    all_species_data |>
    dplyr::distinct(SiteID, Sample_Barcode_ID) |>
    dplyr::count(SiteID, name = "n_reps") |>
    dplyr::filter(n_reps == 3) |>
    dplyr::pull(SiteID)

  df2 <-
    all_species_data |>
    dplyr::filter(SiteID %in% complete_sites)

  pa <-
    df2 |>
    dplyr::select(SiteID, Sample_Barcode_ID, Scientific_Name,
                  sample_detection) |>
    dplyr::group_by(SiteID, Sample_Barcode_ID, Scientific_Name) |>
    dplyr::summarize(sample_detection = max(sample_detection),
                     .groups = "drop") |>
    tidyr::pivot_wider(names_from = Scientific_Name,
                       values_from = sample_detection,
                       values_fill = 0)

  sites <- unique(pa$SiteID)
  boot_all <- list()

  # Bootstrap pooled richness from 1, 2, or 3 replicates within each site.
  for (s in sites) {
    pa_s <- pa |> dplyr::filter(SiteID == s)
    pa_s$SiteID <- NULL
    mat <- as.matrix(pa_s[order(pa_s$Sample_Barcode_ID), -1, drop = FALSE])
    if (!(nrow(mat) == 3)) next

    boot_mat <- matrix(NA, nrow = n_boot, ncol = 3)
    for (b in seq_len(n_boot)) {
      for (k in 1:3) {
        sampled <- sample(1:3, k, replace = TRUE)
        pooled <- colSums(mat[sampled, , drop = FALSE]) > 0
        boot_mat[b, k] <- sum(pooled)
      }
    }
    boot_all[[s]] <- boot_mat
  }

  combined <- do.call(rbind, boot_all)
  reps_summary <-
    data.frame(
      replicates = c(1, 2, 3),
      richness = colMeans(combined),
      lower = apply(combined, 2, stats::quantile, probs = 0.025),
      upper = apply(combined, 2, stats::quantile, probs = 0.975)
    )

  reps_plot <-
    reps_summary |>
    ggplot2::ggplot(ggplot2::aes(x = replicates, y = richness)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                         fill = "gold", alpha = 0.25) +
    ggplot2::geom_line(color = "#B8860B", linewidth = 2) +
    ggplot2::geom_point(color = "#B8860B", size = 5) +
    ggplot2::scale_x_continuous(breaks = c(1, 2, 3)) +
    ggplot2::theme_bw() +
    ggplot2::xlab("Number of qPCR sub samples (K)") +
    ggplot2::ylab("Accumulated species richness")

  ggplot2::ggsave(filename = here::here(fig_dir, "accum_reps.jpg"),
                  plot = reps_plot,
                  dpi = 300,
                  width = 6, height = 6)

  reps_summary$n_sites_used <- length(boot_all)
  readr::write_csv(x = reps_summary,
                   file = here::here(txt_dir, "accum_reps_summary.csv"))
  reps_plot
}
