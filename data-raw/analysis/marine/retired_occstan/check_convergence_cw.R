check_convergence_cw <- function(summaries, rhat_max = 1.01,
                                 n_eff_min = 400) {
  txt_dir <- here::here("data-raw", "analysis", "marine", "outputs")

  parameter <- Rhat <- n_eff <- species <- marker <- subset <- NULL #nolint

  # Pool a named list of per-marker/subset occupancy summaries and check the
  # standard Stan convergence diagnostics across every per-species fit:
  # Rhat < rhat_max and n_eff > n_eff_min. Each list element is tagged with its
  # name as `subset` (e.g. "gp_3_12S", "cw_2_18S"). Returns the flagged table.
  diag <-
    dplyr::bind_rows(summaries, .id = "subset") |>
    dplyr::filter(parameter != "lp__")

  flagged <-
    diag |>
    dplyr::filter(Rhat > rhat_max | n_eff < n_eff_min) |>
    dplyr::select(subset, species, parameter, Rhat, n_eff) |>
    dplyr::arrange(dplyr::desc(Rhat))

  summary_tbl <-
    diag |>
    dplyr::group_by(subset) |>
    dplyr::summarize(n_fits = dplyr::n(),
                     max_Rhat = max(Rhat, na.rm = TRUE),
                     min_n_eff = min(n_eff, na.rm = TRUE),
                     n_flagged = sum(Rhat > rhat_max | n_eff < n_eff_min),
                     .groups = "drop")

  converged <- nrow(flagged) == 0

  sink(file.path(txt_dir, "convergence_report_cw.txt"))
  cat("Coastwide marine occupancy MCMC convergence check\n")
  cat("Thresholds: Rhat <", rhat_max, "and n_eff >", n_eff_min, "\n\n")
  print(as.data.frame(summary_tbl))
  cat("\nOverall converged:", converged, "\n")
  if (!converged) {
    cat("\nFlagged (", nrow(flagged), "):\n", sep = "")
    print(as.data.frame(flagged))
  } else {
    cat("\nAll species-parameters pass both thresholds.\n")
  }
  sink()

  readr::write_csv(x = flagged,
                   file = file.path(txt_dir, "convergence_flagged_cw.csv"))

  list(converged = converged, summary = summary_tbl, flagged = flagged)
}
