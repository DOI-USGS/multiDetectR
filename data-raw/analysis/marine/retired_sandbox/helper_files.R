#' Bootstrap species accumulation curve from empirical sampling events
#'
#' This function summarizes the number of distinct species detected based on
#' specified units, samples, and PCR replicates from the input data.
#'
#' @param n_units An integer specifying the number of units to sample.
#'   Default is 10.
#' @param n_samples An integer specifying the number of samples to sample.
#'   Default is 2.
#' @param k_subsamples An integer specifying the number of PCR subsamples to
#'   sample. Default is 10.
#' @param data_in A data frame containing the sampling data. Default is
#'   `dat_may_3_raw`.
#' @param replace A logical indicating whether sampling should be with
#'   replacement. Default is TRUE.
#' @return A tibble summarizing the count of distinct species detected
#'   along with the parameters used for sampling.
#' @examples
#' species_summary <- species_from_event(n_units = 5, n_samples = 3, k_subsamples = 8,
#'                                       data_in = some_data_frame)
#' species_summary
#' @export
dat_may_3_raw <- NULL

species_from_event  <-
  function(n_units = 10,
           n_samples = 2,
           k_subsamples = 10,
           data_in = dat_may_3_raw,
           replace = TRUE) {
    unit <- NULL
    sample <- NULL
    PCRReplicate <- NULL
    species <- NULL
    y <- NULL

    units_use <-
      data_in |>
      dplyr::distinct(unit) |>
      dplyr::pull(unit) |>
      base::sample(size = n_units, replace = replace) |>
      base::unique()

    samples_use <-
      data_in |>
      dplyr::distinct(sample) |>
      dplyr::pull(sample) |>
      base::sample(size = n_samples, replace = replace) |>
      base::unique()

    pcr_replicate <-
      data_in |>
      dplyr::distinct(PCRReplicate) |>
      dplyr::pull(PCRReplicate) |>
      base::sample(size = k_subsamples, replace = replace) |>
      base::unique()

    data_sample <-
      data_in |>
      dplyr::filter(unit %in% units_use,
                    sample %in% samples_use,
                    PCRReplicate %in% pcr_replicate,
                    y > 0)

    data_summary <-
      tibble::tibble(
        n_units = n_units,
        n_samples = n_samples,
        k_subsamples = k_subsamples,
        n_species = data_sample |>
          dplyr::distinct(species) |>
          nrow()
      )
    return(data_summary)
  }