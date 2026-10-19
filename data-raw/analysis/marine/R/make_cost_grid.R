make_cost_grid <- function(x_sites_vec = 1:10,
                            n_samples_vec = 1:2,
                            k_subsamples_vec = 1:10,
                            site_cost = 100,
                            sample_cost = 100,
                            event_intercept = 500,
                            subsample_rate = 5,
                            n_events = 1) {
  txt_dir <- here::here("data-raw", "analysis", "marine", "outputs")

  X_sites <- N_samples <- K_subsamples <- NULL #nolint

  # Build the full design grid over sites (X), samples (N), and sub-samples (K),
  # and attach per-level costs plus the per-event travel intercept. The default
  # ranges give the 200 Point Conception scenarios (10 x 2 x 10). Ported from
  # Larson's generate_cost_grid in costcalculations_DL.R.
  grid <- expand.grid(X_sites = x_sites_vec,
                      N_samples = n_samples_vec,
                      K_subsamples = k_subsamples_vec,
                      KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)

  total_samples <- grid$X_sites * grid$N_samples
  total_pcr_reps <- total_samples * grid$K_subsamples

  cost_events <- event_intercept * n_events
  cost_sites <- site_cost * grid$X_sites
  cost_samples <- sample_cost * total_samples
  cost_subsamples <- subsample_rate * total_pcr_reps

  cost_grid <-
    grid |>
    dplyr::mutate(total_samples = total_samples,
                  total_pcr_reps = total_pcr_reps,
                  cost_events = cost_events,
                  cost_sites = cost_sites,
                  cost_samples = cost_samples,
                  cost_subsamples = cost_subsamples,
                  total_cost = cost_events + cost_sites +
                    cost_samples + cost_subsamples) |>
    dplyr::arrange(X_sites, N_samples, K_subsamples)

  readr::write_csv(x = cost_grid,
                   file = file.path(txt_dir, "edna_cost_scenarios.csv"))
  cost_grid
}
