fit_species_stan_2 <- function(all_species_data) {
  options(mc.cores = parallel::detectCores())

  revisit <- year <- month <- SiteID <- Scientific_Name <- Common_Name <- #nolint
    sample_detection <- y <- unit <- probs <- parameter <- #nolint
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

  # Format data for occStan
  ## take product to keep everything between 0 and 1
  sample_data <-
    all_species_data |>
    dplyr::group_by(revisit, year, month,
                    SiteID, Scientific_Name, Common_Name) |>
    dplyr::summarize(k_size = 1,
                     y = prod(sample_detection), .groups = "drop") |>
    dplyr::mutate(unit = as.integer(factor(revisit)))

  sample_group_index <-
    sample_data |>
    dplyr::distinct(revisit, SiteID) |>
    dplyr::group_by(revisit) |>
    dplyr::mutate(sample = dplyr::row_number())

  sample_data <-
    sample_data |>
    dplyr::full_join(sample_group_index, by = c("revisit", "SiteID"),
                     relationship = "many-to-many")

  sample_data  |>
    dplyr::select(-Common_Name)

  sample_data |>
    dplyr::select(-Common_Name) |>
    dplyr::filter((k_size > 1))

  sample_data |>
    dplyr::select(-Common_Name) |>
    dplyr::filter((k_size > 1) & (y > 0))

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
      dplyr::filter(Scientific_Name == spp_idx)

    unit_data_use <-
      unit_data |>
      dplyr::filter(Scientific_Name == spp_idx)

    stan_data <-
      occStan::format_occ_2(unit_data = unit_data_use,
                            sample_data = sample_data_use)

    out <-
      occStan::occ_2_stan(stan_data = stan_data,
                          iter = 2000)
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

  readr::write_csv(x = all_coef,
                   file = file.path(data_dir_out, "all_summary.csv"))
  all_coef
}