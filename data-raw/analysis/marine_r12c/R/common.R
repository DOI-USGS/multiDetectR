# Frame 2 (Kevin Lafferty, 2026-10-10): hierarchy marker < PCR < filter < sample (bag) < site < location < visit < coastline.
# Intertidal species ("rocky intertidal" + "partly"). 12S and COI (18S dropped). Prediction for Nanopore only.
# Pt. Conception sites (Government Point, Cojo, Jalama) are BRANCHES (fixed, ranked habitat quality GP > CJ > JB).
# The 20 other locations are EXCHANGEABLE (one site, one visit each) -> branch "COAST" with its own presence probability.
# PCRs per filter are a design choice; pooled reactions share one extract (finite-extract correction, fraction a per PCR).
# Reads ONLY the deposit (or archive) CSVs; writes only under this folder; never overwrites.
suppressPackageStartupMessages(library(dplyr))
# --- PORTED to the multiDetectR targets workflow (2026-10-26) ---
# ROOT was Kevin's Mac path. It now comes from the MULTIDETECTR_R12C_ROOT
# environment variable, which the _targets.R driver sets to this analysis
# directory (the one holding costs.csv, out/, and R/). As a convenience for
# running a script by hand, it falls back to the current working directory.
# ROOT keeps its trailing slash and its original roles: outputs under
# ROOT/out_<label>/, costs at ROOT/../costs.csv (the driver also sets COSTS to
# the right path), and model.R sourced from ROOT (model.R is placed at ROOT by
# the driver; see _targets.R). Scripts' own `source(... common.R)` line is
# rewritten to a relative path by the port (see README_port.md).
ROOT <- paste0(sub("/*$", "", Sys.getenv("MULTIDETECTR_R12C_ROOT", getwd())), "/")
# Default since 2026-10-13 (Kevin: "the code was run on whatever is in Dryad"): the Dryad deposit, files <name>.csv, every read
# checked against the deposit's MANIFEST.csv md5. r2 (dryad_2026-10-13) until 2026-10-25; since then r12c (dryad_2026-10-24, the
# deposit the paper cites; Kevin "yes, do it" 2026-10-25). The older internal archive layout (<name>_<date>_<label>.csv, e.g.
# R1009f) is still read when CI_ARCHIVE points at it and CI_DEPOSIT is set empty; r2 with CI_DEPOSIT=<its path> CI_LABEL=r2 SCREEN=0.
# SCREEN (default 1, since r12c): use the calls as the report products show them. Each <marker>_reads row is joined to the deposit's
# screening_log.csv on (esv_id, Site, marker): taxon/rank become taxon_as_shown/rank_as_shown (expert move-ups, renames) and rows
# with status "dropped" (expert exclusions, blank flags, the excluded tube) go. Tubes with tube_excluded (RSTJ54ZI, study lead
# 2026-10-15) are also removed from samples and libraries, so they count as neither detections nor effort. SCREEN=0 = the calls
# as the pipeline made them (outputs then go to out_<label>_unscreened/).
# PORTED: the Dryad r12c deposit path. CI_DEPOSIT (Kevin's original override)
# still wins if set. Otherwise default to the r12c release located via the
# repo's MULTIDETECTR_DATA convention (the data-release checkout, where the
# marine input lives), falling back to the dryad_2026-10-24 folder at the
# repository root (a sibling of the multiDetectR package).
.default_deposit <- {
  md <- Sys.getenv("MULTIDETECTR_DATA", "")
  cand <- if (nzchar(md)) c(file.path(md, "input", "marine", "dryad_2026-10-24"),
                            file.path(md, "dryad_2026-10-24")) else character(0)
  # repo-root fallback: .../<repo>/dryad_2026-10-24 (two levels above the package's data-raw)
  cand <- c(cand, file.path(getwd(), "dryad_2026-10-24"))
  hit <- cand[file.exists(file.path(cand, "MANIFEST.csv"))]
  if (length(hit)) hit[1] else cand[length(cand)]
}
DEPOSIT <- Sys.getenv("CI_DEPOSIT", .default_deposit)
SCREEN  <- Sys.getenv("SCREEN", "1") == "1"
if (nzchar(DEPOSIT)) {
  ARCHIVE <- DEPOSIT; ARCHIVE_DATE <- sub("^dryad_", "", basename(DEPOSIT))
  RUN_LABEL <- paste0(Sys.getenv("CI_LABEL", "r12c"), if (!SCREEN) "_unscreened" else "")
  MANIFEST_MD5 <- unname(tools::md5sum(file.path(DEPOSIT, "MANIFEST.csv")))
  .man <- utils::read.csv(file.path(DEPOSIT, "MANIFEST.csv"), stringsAsFactors = FALSE)
  .arc_raw <- function(name) { p <- file.path(DEPOSIT, paste0(name, ".csv")); f <- basename(p)
    stopifnot(f %in% .man$file, unname(tools::md5sum(p)) == .man$md5[.man$file == f])          # the file is the deposited one
    utils::read.csv(p, stringsAsFactors = FALSE) }
  .excl <- NULL; .slog <- NULL
  .excluded_tubes <- function() { if (is.null(.excl)) { s <- .arc_raw("samples")
    .excl <<- if ("tube_excluded" %in% names(s)) s$Barcode[s$tube_excluded %in% TRUE] else character(0) }; .excl }
  arc <- function(name) { x <- .arc_raw(name); if (!SCREEN) return(x)
    if (name == "samples") return(x[!x$Barcode %in% .excluded_tubes(), , drop = FALSE])
    if (name == "libraries") return(x[!x$Barcode %in% .excluded_tubes(), , drop = FALSE])
    if (grepl("^(12S|18S|COI)_reads$", name)) { m <- sub("_reads$", "", name)
      if (is.null(.slog)) .slog <<- .arc_raw("screening_log")
      sl <- .slog[.slog$marker == m, c("esv_id", "Site", "taxon_as_shown", "rank_as_shown", "status")]
      n0 <- nrow(x); x <- dplyr::left_join(x, sl, by = c("esv_id", "Site")); stopifnot(nrow(x) == n0, !anyNA(x$status))   # one log row per ESV x site
      x <- x[x$status == "kept" & !x$Barcode %in% .excluded_tubes(), , drop = FALSE]
      x$taxon <- x$taxon_as_shown; x$rank <- x$rank_as_shown; return(x[, setdiff(names(x), c("taxon_as_shown", "rank_as_shown", "status"))]) }
    x }
} else {
  ARCHIVE <- Sys.getenv("CI_ARCHIVE", paste0(ROOT, "archive"))  # PORTED: old-archive fallback, only used when CI_DEPOSIT is empty (not in the targets workflow)
  .f <- basename(Sys.glob(file.path(ARCHIVE, "samples_*.csv"))); stopifnot(length(.f) == 1)
  ARCHIVE_DATE <- sub("^samples_(\\d{4}-\\d{2}-\\d{2})_.*$", "\\1", .f)
  RUN_LABEL    <- sub("^samples_\\d{4}-\\d{2}-\\d{2}_(.*)\\.csv$", "\\1", .f); MANIFEST_MD5 <- NA_character_
  arc <- function(name) utils::read.csv(file.path(ARCHIVE, sprintf("%s_%s_%s.csv", name, ARCHIVE_DATE, RUN_LABEL)), stringsAsFactors = FALSE)
}
OUT <- paste0(ROOT, "out_", RUN_LABEL, "/"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
prov <- function() data.frame(data = RUN_LABEL, source = ARCHIVE, manifest_md5 = MANIFEST_MD5,
  calls = if (SCREEN) "screened (screening_log.csv; excluded tubes removed)" else "as made by the pipeline", written = format(Sys.time(), "%Y-%m-%d %H:%M %Z"))
SITES    <- c("GP", "CJ", "JB", "COAST")              # model branches; COAST = any of the 20 exchangeable locations
PC_SITES <- c(`Government Point` = "GP", `Cojo` = "CJ", `Jalama Beach` = "JB")
MARKERS  <- c("12S", "COI")
HABITAT_KEEP <- c("rocky intertidal", "partly")
PLATFORM <- c(JVB2844 = "illumina", JVB3105 = "illumina", JVB3506 = "illumina")   # every other batch is Nanopore (workflow PLATFORM_OF_BATCH)
REF_BATCH <- "JVB6097"                                   # current protocol (Sep 2025 - Jan 2026); design PCRs are in its units
A_GRID <- c(0.01, 0.03, 0.05)                            # fraction of the extract used per PCR (peer bracket; elution volume unknown)
# fits and build always come from FITDIR. Env COSTS = another cost table, COSTTAG = its tag: design outputs then go to OUT/costs_<tag>/.
# Default since 2026-10-11: ../costs.csv = technician-hour vector, outputs in OUT/costs_hours_2026-10-11/. The earlier $400 vector is
# ../costs_usd_2026-10-10.csv; its outputs are the top level of OUT (reproduce with COSTS=<that file> COSTTAG= set but empty).
FITDIR <- OUT; COSTS_FILE <- Sys.getenv("COSTS", paste0(ROOT, "../costs.csv")); COSTTAG <- Sys.getenv("COSTTAG", "hours_2026-10-11")
if (nzchar(COSTTAG)) { OUT <- paste0(OUT, "costs_", COSTTAG, "/"); dir.create(OUT, showWarnings = FALSE) }
# costs: editable table in ../costs.csv (default)
read_costs <- function(p = COSTS_FILE) { x <- utils::read.csv(p, stringsAsFactors = FALSE); setNames(x$cost_usd, x$item) }
lg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), ..., "\n", file = stderr())
wcsv <- function(x, name) { p <- paste0(OUT, name, ".csv"); if (file.exists(p)) stop("exists, not overwritten: ", p); utils::write.csv(x, p, row.names = FALSE, na = ""); invisible(p) }
# effective single-PCR count of a pool of k reactions, each drawing a fraction a of the extract: (1 - (1-a)^k) / a  (-> k as a -> 0)
m_eff <- function(k, a) (1 - (1 - a)^k) / a
