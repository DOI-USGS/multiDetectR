get_dryad_marine <- function(reads_file, samples_file,
                              platform = c("illumina", "nanopore", "all"),
                              area = NULL,
                              site = NULL,
                              rank_keep = c("species", "genus")) {
  platform <- match.arg(platform)

  Barcode <- Site <- Area <- batch <- esv_id <- library_id <- #nolint
    reads <- taxon <- rank <- unit <- sample <- PCRReplicate <- #nolint
    species <- y <- platform_group <- NULL #nolint

  # Reader for the coastwide marine Dryad release. Produces the pipeline's
  # unit / sample / PCRReplicate / species / y shape for one marker, where:
  #   unit        = Site                (occupancy level, X / psi)
  #   sample      = Barcode (filter)    (field-sample level, N / theta)
  #   PCRReplicate = sequencing library (sub-replicate level, K / p)
  #   species     = the per-site taxonomic call carried in the reads file
  #   y           = 1 (a reads > 0 row is a detection)
  #
  # The reads file lists only reads > 0, so each returned row is a detection;
  # non-detection is absence from the frame, which the occupancy and
  # accumulation functions reconstruct from the sampled units. Analyses are
  # restricted to a platform-consistent subset (decision: restrict, not model).

  # Map the Jonah Ventures batch to sequencing platform. Illumina: the MiSeq /
  # NovaSeq runs through May 2024 (Government Point / Dangermond). Nanopore: the
  # Summer 2024 coastwide expedition (JVB3735) and all later JVB5xxx / JVB6xxx.
  illumina_batches <- c("JVB2844", "JVB3105", "JVB3506")

  samples <-
    readr::read_csv(samples_file, col_types = readr::cols()) |>
    dplyr::select(Barcode, Site, Area, bag, SubSite, batch = batch_12S) |>
    dplyr::mutate(platform_group = ifelse(batch %in% illumina_batches,
                                          "illumina", "nanopore"))

  reads <-
    readr::read_csv(reads_file, col_types = readr::cols()) |>
    dplyr::filter(reads > 0, rank %in% rank_keep) |>
    dplyr::select(esv_id, library_id, Barcode, Site, taxon, rank)

  out <-
    reads |>
    dplyr::inner_join(samples, by = c("Barcode", "Site"))

  if (platform != "all") {
    out <- dplyr::filter(out, platform_group == platform)
  }
  if (!is.null(area)) {
    out <- dplyr::filter(out, Area %in% area)
  }
  if (!is.null(site)) {
    out <- dplyr::filter(out, Site %in% site)
  }

  out |>
    dplyr::transmute(unit = Site,
                     sample = Barcode,
                     PCRReplicate = library_id,
                     species = taxon,
                     consensus = taxon,
                     y = 1L)
}
