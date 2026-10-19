fit_cw_2 <- function(marine_data, label) {
  options(mc.cores = parallel::detectCores())

  unit <- sample <- species <- consensus <- y <- k_size <- #nolint
    n_units_per_site <- z_obs <- parameter <- #nolint
    `50%` <- `2.5%` <- `5%` <- `95%` <- `97.5%` <- NULL #nolint

  rstan::rstan_options(threads_per_chain = 1)

  data_dir_out <- here::here("data-raw", "analysis", "marine", "data")

  # Two-level occupancy for one marker/subset: unit (Site) -> sample (Barcode
  # filter), collapsing PCR libraries to a single sample-level detection
  # (k_size = 1). Estimates site occupancy (psi) and sample-level detection (p).

  # Roster of samples attempted per unit (every sample appears in the reader
  # output because every library has at least one in-scope detection).
  roster <-
    marine_data |>
    dplyr::distinct(unit, sample)

  # Per sample x species: detected at the sample level or not.
  detect <-
    marine_data |>
    dplyr::distinct(unit, sample, species, consensus)

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
      dplyr::mutate(
        y = as.integer(paste(unit, sample) %in%
                         with(dplyr::filter(detect, species == spp_idx),
                              paste(unit, sample))),
        k_size = 1L,
        species = spp_idx)

    unit_data_use <-
      sample_data_use |>
      dplyr::group_by(unit, species) |>
      dplyr::summarize(z_obs = ifelse(sum(y) > 0, 1, 0), .groups = "drop")

    if (sum(unit_data_use$z_obs) == 0) next
    print(spp_idx)

    stan_data <-
      occStan::format_occ_2(unit_data = unit_data_use,
                            sample_data = sample_data_use)

    out <-
      occStan::occ_2_stan(stan_data = stan_data, iter = 2000)
    pars_save <- c("beta_psi", "delta_p", "lp__")

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
                                    paste0("cw_summary_2_", label, ".csv")))
  all_coef
}
