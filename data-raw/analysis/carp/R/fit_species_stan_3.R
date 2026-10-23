fit_species_stan_3 <- function(all_species_data) {
  options(mc.cores = parallel::detectCores())

  revisit <- year <- month <- SiteID <- Scientific_Name <- Common_Name <- #nolint
    sample_detection <- y <- unit <- probs <- parameter <- sample_occ <- #nolint
    unit_factor <- #nolint
    `50%` <- `2.5%` <- `5%` <- `95%` <- `97.5%` <- k_size <- NULL #nolint

  rstan::rstan_options(threads_per_chain = 1)

  # Write regenerated summaries to the data-release repo (output/), not the
  # package repo. Mirrors the MULTIDETECTR_DATA convention used by _targets.R.
  release_dir <- Sys.getenv(
    "MULTIDETECTR_DATA",
    unset = file.path(dirname(here::here()), "data-release")
  )
  data_dir_out <- file.path(release_dir, "output", "carp")
  dir.create(data_dir_out, recursive = TRUE, showWarnings = FALSE)

  # Format data for occStan three-level model.
  #
  # The data nest into three levels:
  #   - unit (occupancy, psi):      revisit (a site-visit)
  #   - sample (availability, theta): SiteID within a revisit (a water sample)
  #   - sub-replicate (detection, p): Sample_Barcode_ID (PCR replicate)
  #
  # One row per sample: k_size counts the sub-replicates (PCR replicates),
  # y counts the positive sub-replicates, and sample_occ flags any detection.
  sample_data <-
    all_species_data |>
    dplyr::group_by(revisit, year, month,
                    SiteID, Scientific_Name, Common_Name) |>
    dplyr::summarize(k_size = dplyr::n(),
                     y = sum(sample_detection), .groups = "drop") |>
    dplyr::mutate(sample_occ = as.integer(y > 0),
                  unit = as.integer(factor(revisit)),
                  unit_factor = factor(unit))

  sample_group_index <-
    sample_data |>
    dplyr::distinct(revisit, SiteID) |>
    dplyr::group_by(revisit) |>
    dplyr::mutate(sample = dplyr::row_number())

  sample_data <-
    sample_data |>
    dplyr::full_join(sample_group_index, by = c("revisit", "SiteID"),
                     relationship = "many-to-many")

  unit_data <-
    sample_data |>
    dplyr::group_by(Scientific_Name, unit) |>
    dplyr::summarize(unit_occ = ifelse(sum(y > 0), 1, 0),
                     .groups = "drop")

  # 3. loop over all species
  all_species_names <-
    sample_data |>
    dplyr::distinct(Scientific_Name) |>
    dplyr::pull(Scientific_Name)
  all_coef <- tibble::tibble()

  for (spp_idx in all_species_names){
    print(spp_idx)
    sample_data_use <-
      sample_data |>
      dplyr::filter(Scientific_Name == spp_idx) |>
      dplyr::arrange(unit)

    unit_data_use <-
      unit_data |>
      dplyr::filter(Scientific_Name == spp_idx) |>
      dplyr::arrange(unit)

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
      occStan::occ_3_stan(stan_data = stan_data,
                          iter = 2000)
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

  readr::write_csv(x = all_coef,
                   file = file.path(data_dir_out, "all_summary_3.csv"))
  all_coef
}
