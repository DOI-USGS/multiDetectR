# r12c gap-1 rebuild wrapper (sampling-opt, 2026-10-25): the deposited build_reads_long_hopflagged.R with PATHS ONLY changed:
# code from code_2026-10-24_r12c.zip, data from the Dryad deposit raw_inputs.zip (both unzipped to scratch); one marker per run (env MK).
# =============================================================================
# INGEST-LEVEL index-hop flagging, run ONCE per marker and shared by every consumer.
#
# WHY A SEPARATE SCRIPT. make_esv_site_detections.R, make_esv_sample_detections.R and the
# workflow each call pc_build_inputs() independently. Adding the hop step inside each would run it
# three times and let the three copies drift. This builds the flagged reads_long once and writes it.
#
# NOT WIRED IN (corrected 2026-10-08). The plan above was that the others would read this file. They
# never did: the workflow and both make_esv_*_detections.R still call pc_build_inputs() directly, and
# DROP_FLAGS below is empty. So the hop flag is DIAGNOSTIC ONLY and has never filtered a detection or
# a call in any marker. The only readers are studies_2026-10-07/s1_sampling*/01_build.R. Rebuild this
# whenever the reads change (e.g. a library reissue): the 18S file was stale for JVB6164/JVB6334 until
# 2026-10-08.
#
# WHAT IT IS FOR, stated plainly so nobody over-reads the output: this dataset has essentially NO
# index hopping. Measured across 35 run x marker combinations, most runs are bound_only with an
# upper bound near 1e-5; the 12S ESVs on JVB3105 with a median 1.9M reads appear in 0 of 9 blanks,
# where 0.1% hopping predicts about 4 reads per blank. The value here is the RECORD -- an explicit,
# measured statement that hopping was assessed and bounded -- not the rows it removes.
#
# FILTER POLICY = (c) KEEP EVERYTHING, CARRY validity_flag AS PROVENANCE.
# The only substantial flags are ~1.34% of 12S detections, all on JVB3735, and the shape of those
# reads says clumped well-to-well or handling contamination, NOT diffuse hopping: the ten biggest
# ESVs occupy a median 5% of that run's 29 blanks (three of them 0%), where diffuse hopping at the
# fitted 0.15% would put reads of each into essentially every blank. Removing data on a mechanism
# the evidence contradicts is worse than carrying a column. Switch to (a) or (b) by editing
# DROP_FLAGS below -- and note NA means NOT ASSESSED (refused implausible_rate runs, no evidence,
# or control rows) and must never be filtered as bad. Hence !(flag %in% DROP) and never == "valid".
# =============================================================================
if (!exists(".ci_paths_loaded")) local({ d <- getwd(); while (!file.exists(file.path(d, "CaliforniaIntertidal", "paths.R")) && dirname(d) != d) d <- dirname(d); source(Sys.getenv("CI_PATHS", file.path(d, "CaliforniaIntertidal", "paths.R")), local = globalenv()) })   # machine paths: CaliforniaIntertidal/paths.R
DROP_FLAGS <- character(0)   # (a) "invalid_index_hop"; (b) add "questionable_index_hop"

setwd("/private/tmp/claude-501/-Users-lafferty/3861b05b-4942-4975-afc5-5207b58e09cb/scratchpad/r12c_code/code_2026-10-24/code/PtConception/Multi_site_multi_marker")
mem.maxVSize(4000)
suppressPackageStartupMessages({library(dplyr)})
ARCHIVE_ROOT <- "/private/tmp/claude-501/-Users-lafferty/3861b05b-4942-4975-afc5-5207b58e09cb/scratchpad/r12c_raw/raw_inputs"; CROSSWALK <- file.path(ARCHIVE_ROOT, "derived", "ptcon_site_crosswalk.csv"); CORRECTIONS <- file.path(ARCHIVE_ROOT, "derived", "ptcon_field_corrections.csv")
suppressMessages({source("ptcon_ingest.R")})
CI <- CI_DIR
Sys.setenv(TAXAID_LIBRARY_DROPS = "/private/tmp/claude-501/-Users-lafferty/3861b05b-4942-4975-afc5-5207b58e09cb/scratchpad/r12c_raw/raw_inputs/decisions/CalIntertidal_library_drops_2026-10-09.csv"); source(file.path(CI, "library_drops.R"))   # Kevin 2026-10-09: drop failed/re-run libraries + 12S JVB3506 block A
OUT_HOP <- Sys.getenv("TAXAID_HOP_OUT_DIR", CI)   # override = write elsewhere (e.g. scratch) instead of the live files
MARKERS <- c(MiFishU = "12S", `18Sv9_89` = "18S", UniCOI = "COI")[Sys.getenv("MK")]; stopifnot(!anyNA(MARKERS))

if (!requireNamespace("TaxaFlag", quietly = TRUE) ||
    !exists("flag_hopped_detections", asNamespace("TaxaFlag"))) {
  stop("flag_hopped_detections() is not in the installed TaxaFlag (built ",
       packageDescription("TaxaFlag")$Built, ").\n",
       "It is merged to main as 69feacd but NOT yet installed; the package-review session owns\n",
       "the reinstall gate. Re-run this after they confirm. Nothing has been written.")
}

# CONTROLS ARE CHOSEN EXPLICITLY. flag_hopped_detections() is hops-only by decision; it does not
# drop contaminated blanks, so the control set is ours to curate. S067800 is excluded because
# validate_controls() returns RESEMBLES_SAMPLE on all three of its replicates in BOTH 18S and COI
# -- two independent assays of one physical tube -- and its contents are a marine plankton
# community (Sinocalanus, Tripos, Lingulodinium polyedra, Akashiwo sanguinea, tintinnids) rather
# than the lab/human profile of the clean blank S067808 beside it on the same run.
# THIS LIST IS A FLOOR, NOT A CEILING: on 18S JVB3735/JVB3506/JVB2844 and COI JVB3735 the
# blank-versus-field contrast is collapsed, so "consistent_with_control" there means untested.
EXCLUDE_CONTROLS <- c("S067800")
# 2026-10-08 (Kevin "B yes"): the hand list is replaced by the identity decisions -- every tube with disposition
# "exclude_tube" (13, S067800 among them) is a contaminated control and is not used as one. That is change 1 of the
# two below; change 2 (ubiquitous PositiveControl features) is still open.
IDD <- "/private/tmp/claude-501/-Users-lafferty/3861b05b-4942-4975-afc5-5207b58e09cb/scratchpad/r12c_raw/raw_inputs/decisions/CalIntertidal_identity_decisions_2026-09-28.csv"
.idd <- read.csv(IDD, stringsAsFactors = FALSE)
EXCLUDE_CONTROLS <- sort(unique(.idd$sample[.idd$disposition %in% "exclude_tube"]))
stopifnot("identity decisions list no exclude_tube" = length(EXCLUDE_CONTROLS) > 0, "S067800" %in% EXCLUDE_CONTROLS)
#
# TWO CHANGES REQUIRED ONCE THE SAMPLE-IDENTITY GATE LANDS (branch sample-identity-gate):
# 1. Replace EXCLUDE_CONTROLS entirely with the gate's admit_as == "control" set. The hand list
#    above is already known to be insufficient -- 30 of 30 blanks on 18S JVB3735 hold protist
#    communities richer than that run's field samples, with Gambierdiscus at 2,162 reads in blanks
#    against 269 across the whole field set. Two competing definitions of "which blanks are usable"
#    is the failure shape this project keeps reproducing, so there must be exactly one.
# 2. Set aside the gate's runs$ubiquitous_taxa features BEFORE hop detection on the runs that
#    carry the ubiquitous PositiveControl ESV. MEASURED: a single 180 bp ESV, byte-identical
#    across runs and labelled PositiveControl in read-data (it is ESV_075708 on JVB6097), is
#    present in 100% of field libraries and 97-100% of blanks on JVB5058/5059/5060/5061, 6164,
#    6334 and 6097, and entirely absent (0%) from JVB2844, JVB3105, JVB3506 and JVB3735 --
#    verified across all 12 MiFishU runs. So it appears from Event 5 onward and not in Events 1-4.
#    ITS CAUSE IS UNKNOWN AND IS A QUESTION FOR JV. The user did NOT spike these samples, so this
#    is most likely JV's own positive control reaching every library, either added by them or
#    leaking from their positive-control well. Do not describe it as a spike-in by us.
#    Why it matters here regardless of cause: a feature present in every field library AND every
#    blank at high abundance is the textbook hopping signature and would read as catastrophic
#    carryover. Nothing reported so far is affected -- every measured 12S hop flag is on JVB3735,
#    which does not carry it.

rates <- list()
for (mk in names(MARKERS)) {
  short <- MARKERS[[mk]]
  message(sprintf("\n--- %s ---", short))
  inp <- pc_build_inputs(ARCHIVE_ROOT, marker = mk, crosswalk = CROSSWALK, corrections = CORRECTIONS)
  rl <- inp$reads_long |>
    mutate(is_blank     = dplyr::coalesce(is_blank, FALSE),
           is_lab_blank = dplyr::coalesce(is_lab_blank, FALSE),
           observation_id = paste0(short, "_", ESVId))

  # EMPTY BLANKS MUST BE REPRESENTED. pc_read_reads() keeps only n_reads > 0, so a blank that
  # sequenced nothing has NO rows at all -- yet those are the best evidence a run is clean, since a
  # control with zero reads of an abundant ESV is precisely what "no hopping" looks like. Without a
  # row the run looks merely unmeasured instead of clean. 7 of 59 blanks have no 12S reads and 2
  # have no COI reads.
  #
  # The batch for such a blank cannot come from reads_long (it has no rows) and the manifest maps
  # FILE -> batch, not barcode -> batch. The authoritative source is the read file's own header: a
  # zero-read blank still exists there as an all-zero column. So read header lines only (cheap) and
  # recover barcode -> batch from the column names of each marker's read-data files.
  mfiles <- inp$manifest |> filter(marker == mk, grepl("read-data", path)) |> distinct(path, batch)
  hdr <- lapply(seq_len(nrow(mfiles)), function(i) {
    h <- try(names(read.csv(mfiles$path[i], nrows = 1, check.names = FALSE)), silent = TRUE)
    if (inherits(h, "try-error")) return(NULL)
    j <- which(h == "# species"); if (!length(j) || j[1] >= length(h)) return(NULL)
    data.frame(event_id = h[(j[1]+1L):length(h)], batch = mfiles$batch[i], stringsAsFactors = FALSE)
  })
  hdr <- dplyr::bind_rows(hdr)
  if (nrow(hdr)) {
    hdr$Barcode <- sub("\\..*$", "", hdr$event_id)
    absent <- hdr |> filter(Barcode %in% inp$blank_barcodes, !(event_id %in% rl$event_id))
    if (nrow(absent)) {
      # one zero-count row each, on an ESV that exists for this marker so the join keys are valid
      anchor <- rl$observation_id[1]
      add <- absent |> transmute(event_id, Barcode, batch, observation_id = anchor,
                                 ESVId = sub(paste0("^", short, "_"), "", anchor),
                                 n_reads = 0, is_blank = TRUE, is_lab_blank = FALSE)
      for (cc in setdiff(names(rl), names(add))) add[[cc]] <- NA
      rl <- dplyr::bind_rows(rl, add[, names(rl), drop = FALSE])
      message(sprintf("  %d blank(s) had NO reads -- added zero rows so they count as evidence: %s",
                      nrow(absent), paste(utils::head(absent$event_id, 6), collapse = ", ")))
    } else message("  every blank has at least one read row")
  } else message("  NOTE: could not read file headers -- absent blanks NOT represented")

  rl <- drop_listed_libraries(rl, short)   # after the zero-read blanks are added, so listed empty blanks go too
  ctl <- unique(rl$event_id[rl$is_blank | rl$is_lab_blank])
  drop_ctl <- grep(paste0("^X?(", paste(EXCLUDE_CONTROLS, collapse="|"), ")\\."), ctl, value = TRUE)
  ctl <- setdiff(ctl, drop_ctl)
  message(sprintf("  controls: %d used, %d excluded (%s)", length(ctl), length(drop_ctl),
                  if (length(drop_ctl)) paste(drop_ctl, collapse=", ") else "none"))

  hop <- TaxaFlag::flag_hopped_detections(
    rl, event_col = "event_id", taxon_col = "observation_id",
    count_col = "n_reads", run_col = "batch",          # batch, NOT event: JVB6097 spans Events 10/11
    control_samples = ctl)

  rr <- attr(hop, "run_rates"); rates[[short]] <- rr
  if (!is.null(rr)) { message("  run rates:"); print(as.data.frame(rr), row.names = FALSE) }
  ca <- attr(hop, "controls_absent")
  if (!is.null(ca) && length(ca)) message(sprintf("  controls absent (no reads): %s", paste(ca, collapse=", ")))
  vf <- if ("validity_flag" %in% names(hop)) hop$validity_flag else rep(NA_character_, nrow(hop))
  message(sprintf("  flags: %s", paste(names(table(vf, useNA="ifany")), table(vf, useNA="ifany"), sep="=", collapse=" ")))

  if (length(DROP_FLAGS)) {
    n0 <- nrow(hop); hop <- hop[!(vf %in% DROP_FLAGS), , drop = FALSE]
    message(sprintf("  DROP_FLAGS applied: %s row(s) removed", format(n0 - nrow(hop), big.mark=",")))
  } else message("  DROP_FLAGS empty -- nothing removed, validity_flag carried as provenance")

  saveRDS(hop, file.path(OUT_HOP, sprintf("CalIntertidal_reads_long_hopflagged_%s.rds", short)))
  message(sprintf("  wrote CalIntertidal_reads_long_hopflagged_%s.rds (%s rows)", short,
                  format(nrow(hop), big.mark = ",")))
}
saveRDS(list(rates = rates, drop_flags = DROP_FLAGS, excluded_controls = EXCLUDE_CONTROLS,
             built = Sys.time(), taxaflag = packageDescription("TaxaFlag")$Built),
        file.path(OUT_HOP, sprintf("CalIntertidal_hop_run_rates_%s.rds", Sys.getenv("MK"))))
message("\nwrote CalIntertidal_hop_run_rates.rds -- this is the provenance line for the richness and gradient tables")
