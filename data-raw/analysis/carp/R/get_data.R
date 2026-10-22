#' Read and reshape the Burns Harbor freshwater metabarcoding data
#'
#' The raw data release (`Burns_Harbor_Metabarcoding_2023-2025.csv`) lists only
#' detected species (one row per observed barcode x species), with no explicit
#' zero-detection rows, and names its columns differently from what the rest of
#' the carp pipeline expects. This reader reshapes that detections-only file
#' into the complete presence/absence grid the pipeline consumes:
#'
#'   Sample_Barcode_ID, Scientific_Name, SiteID, Read_Count, Common_Name,
#'   year, month, day, revisit, sample_detection
#'
#' Reshaping steps:
#'   1. Rename columns (SiteID_ID -> SiteID, Year -> year, etc.).
#'   2. Derive year/month/day from the timestamp encoded in the site ID
#'      (a 3-letter prefix followed by YYYYMMDD and a time of day), and
#'      revisit = "year_month" (the occupancy unit, a site-visit).
#'   3. Expand to the full species x barcode grid so every barcode carries a
#'      row for every species in the pooled species list, with
#'      sample_detection = 1 where the species was observed for that barcode
#'      and 0 otherwise (Read_Count carried through, 0 for non-detections).
#'
#' The species pool is every Scientific_Name observed across all years in the
#' file (2023-2025), matching the pipeline's universal-pool convention.
#'
#' If the file is already in the reshaped presence/absence format (it contains a
#' `sample_detection` column), it is returned unchanged, so the reader is
#' backwards compatible with a pre-reshaped `all_species.csv`.
get_data <- function(file) {
  raw <- readr::read_csv(file, col_types = readr::cols())

  # Backwards compatibility: a file already in the pipeline format is returned
  # as-is (it already has explicit zeros and the derived columns).
  if ("sample_detection" %in% names(raw)) {
    return(raw)
  }

  SiteID_ID <- SiteID <- Sample_Barcode_ID <- Scientific_Name <- #nolint
    Common_Name <- Read_Count <- Year <- year <- month <- day <- #nolint
    revisit <- sample_detection <- NULL #nolint

  # --- 1. rename to the pipeline's column names ---
  detections <-
    raw |>
    dplyr::rename(SiteID = SiteID_ID) |>
    dplyr::mutate(year = as.integer(Year))

  # --- 2. derive month/day from the timestamp in the site ID ---
  # Site IDs look like "GRB202307301516" or "FWS20250707142828":
  # a 3+ letter prefix, then YYYYMMDD, then a time of day (HHMM or HHMMSS).
  # Take the 8 digits immediately after the leading alphabetic prefix as the
  # date; ignore the trailing time (the pipeline only needs y/m/d and revisit).
  date_digits <- sub("^[A-Za-z]+", "", detections$SiteID)
  detections <-
    detections |>
    dplyr::mutate(
      month = as.integer(substr(date_digits, 5, 6)),
      day   = as.integer(substr(date_digits, 7, 8)),
      revisit = paste(year, month, sep = "_")
    )

  # Guard: the date parse must yield valid calendar months and days.
  bad <- which(is.na(detections$month) | detections$month < 1 |
                 detections$month > 12 | is.na(detections$day) |
                 detections$day < 1 | detections$day > 31)
  if (length(bad) > 0) {
    stop("get_data(): could not parse a valid date from ",
         length(bad), " site ID(s), e.g. '",
         detections$SiteID[bad[1]], "'. Check the site-ID timestamp format.")
  }

  # --- 3. expand detections to the full species x barcode grid ---
  # One Common_Name per Scientific_Name (first non-missing), so the expanded
  # grid can carry the common name on every row.
  species_pool <-
    detections |>
    dplyr::group_by(Scientific_Name) |>
    dplyr::summarize(Common_Name = dplyr::first(Common_Name),
                     .groups = "drop")

  # The identity of each barcode: its site, revisit, and date. A barcode belongs
  # to exactly one site-visit, so these are constant within a barcode.
  barcode_info <-
    detections |>
    dplyr::distinct(Sample_Barcode_ID, SiteID, year, month, day, revisit)

  # The observed detections, keyed by barcode x species, with their read counts.
  observed <-
    detections |>
    dplyr::group_by(Sample_Barcode_ID, Scientific_Name) |>
    dplyr::summarize(Read_Count = sum(Read_Count), .groups = "drop") |>
    dplyr::mutate(sample_detection = 1L)

  # Full grid: every barcode crossed with every species in the pool.
  grid <-
    tidyr::expand_grid(
      Sample_Barcode_ID = barcode_info$Sample_Barcode_ID,
      Scientific_Name = species_pool$Scientific_Name
    ) |>
    dplyr::left_join(barcode_info, by = "Sample_Barcode_ID") |>
    dplyr::left_join(species_pool, by = "Scientific_Name") |>
    dplyr::left_join(observed, by = c("Sample_Barcode_ID", "Scientific_Name")) |>
    dplyr::mutate(
      Read_Count = dplyr::coalesce(Read_Count, 0),
      sample_detection = dplyr::coalesce(sample_detection, 0L)
    )

  # Return the columns the pipeline expects, in the original order.
  grid |>
    dplyr::select(Sample_Barcode_ID, Scientific_Name, SiteID, Read_Count,
                  Common_Name, year, month, day, revisit, sample_detection) |>
    dplyr::arrange(Sample_Barcode_ID, Scientific_Name)
}
