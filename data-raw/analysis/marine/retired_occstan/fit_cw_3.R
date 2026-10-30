fit_cw_3 <- function(marine_data, label) {
  options(mc.cores = parallel::detectCores())

  unit <- sample <- PCRReplicate <- species <- consensus <- y <- #nolint
    k_size <- sample_occ <- n_units_per_site <- z_obs <- parameter <- #nolint
    `50%` <- `2.5%` <- `5%` <- `95%` <- `97.5%` <- NULL #nolint

  rstan::rstan_options(threads_per_chain = 1)

  # Write regenerated summaries to the data-release repo (output/), not the
  # package repo. Mirrors the MULTIDETECTR_DATA convention used by _targets.R.
  release_dir <- Sys.getenv(
    "MULTIDETECTR_DATA",
    unset = file.path(dirname(here::here()), "data-release")
  )
  data_dir_out <- file.path(release_dir, "output", "marine")
  dir.create(data_dir_out, recursive = TRUE, showWarnings = FALSE)

  # Three-level occupancy for one marker/subset of the Dryad release.
  # Levels: unit (Site) -> sample (Barcode filter) -> PCR library (K).
  # The reader lists only detections, but every sequenced library has at least
  # one in-scope detection, so distinct (sample, PCRReplicate) over all species
  # is the full roster of PCR libraries attempted per sample.

  # Roster: number of PCR libraries attempted per sample (k_size).
  roster <-
    marine_data |>
    dplyr::distinct(unit, sample, PCRReplicate) |>
    dplyr::group_by(unit, sample) |>
    dplyr::summarize(k_size = dplyr::n(), .groups = "drop")

  # Per sample x species: y = number of PCR libraries that detected the species.
  detect <-
    marine_data |>
    dplyr::distinct(unit, sample, species, consensus, PCRReplicate) |>
    dplyr::group_by(unit, sample, species, consensus) |>
    dplyr::summarize(y = dplyr::n(), .groups = "drop")

  # Build the full sample x species frame: every species gets a row for every
  # sampled unit/sample (y = 0 where undetected), with k_size from the roster.
  all_species <- sort(unique(marine_data$species))
  spp_key <-
    marine_data |>
    dplyr::distinct(species, consensus) |>
    dplyr::group_by(species) |>
    dplyr::summarize(consensus = dplyr::first(consensus), .groups = "drop")

  all_coef <- tibble::tibble()

  for (spp_idx in all_species) {
    sample_data_use <-
      roster |>
      dplyr::left_join(
        detect |>
          dplyr::filter(species == spp_idx) |>
          dplyr::select(unit, sample, y),
        by = c("unit", "sample")) |>
      dplyr::mutate(y = ifelse(is.na(y), 0L, as.integer(y)),
                    sample_occ = as.integer(y > 0),
                    species = spp_idx)

    unit_data_use <-
      sample_data_use |>
      dplyr::group_by(unit, species) |>
      dplyr::summarize(z_obs = ifelse(sum(y) > 0, 1, 0), .groups = "drop")

    # Skip species never detected in this subset.
    if (sum(unit_data_use$z_obs) == 0) next
    print(spp_idx)

    stan_data <-
      occStan::format_occ_3(unit_data = unit_data_use,
                            sample_data = sample_data_use,
                            psi_formula = ~1,
                            theta_formula = ~1,
                            p_formula = ~1,
                            unit_name = "unit",
                            k_size_name = "k_size",
                            sample_seen_name = "sample_occ",
                            y_name = "y")

    out <-
      occStan::occ_3_stan(stan_data = stan_data, iter = 2000)
    pars_save <- c("beta_psi", "alpha_theta", "delta_p", "lp__")

    all_coef <-
      dplyr::bind_rows(
        all_coef,
        rstan::summary(out, pars = pars_save,
                       probs = c(0.025, 0.05, 0.5, 0.95, 0.975))$summary |>
          as.data.frame() |>
          tibble::rownames_to_column("parameter") |>
          dplyr::mutate(parameter = gsub("\\[1\\]", "", parameter),
                        mean_prob = plogis(mean),
                        med_prob = plogis(`50%`),
                        l95 = plogis(`2.5%`),
                        l90 = plogis(`5%`),
                        u90 = plogis(`95%`),
                        u95 = plogis(`97.5%`),
                        species = spp_idx)
      )
  }

  all_coef <-
    all_coef |>
    dplyr::left_join(spp_key, by = "species") |>
    dplyr::mutate(marker = label)

  readr::write_csv(x = all_coef,
                   file = file.path(data_dir_out,
                                    paste0("cw_summary_3_", label, ".csv")))
  all_coef
}
