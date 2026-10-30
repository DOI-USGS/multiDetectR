#!/usr/bin/env Rscript
# verify_against_handoff.R -- acceptance test for the r12c port.
#
# Kevin's analysis is deterministic (no seeds), so a correct port must regenerate
# his result CSVs byte-for-byte in value. This script compares every CSV this
# pipeline produced against the matching file in handoff/ and reports, per file,
# whether the numbers agree.
#
# Two output locations are checked (see README_port.md):
#   out_r12c/        -- the t_*.csv summary/design tables and the fit .rds files.
#   handoff_out/     -- the handoff CSVs written by scripts 12/30/31/32, which
#                       use Kevin's exact filenames.
#
# Matching strategy (no hardcoded name map, so it can't silently rot):
#   1. handoff_out/<f>.csv            <-> handoff/<f>.csv            (exact name)
#   2. out_r12c/t_<f>.csv             <-> handoff/<f>.csv            (drop the t_ prefix)
#   3. out_r12c/t_<f>.csv             <-> handoff/cmp_<f>.csv        (t_ -> cmp_)
# A handoff file matched by more than one rule is reported once (first match wins,
# in the order above). Files on only one side are listed so nothing is missed.
#
# USAGE (from marine_r12c/, after tar_make() finishes):
#   Rscript verify_against_handoff.R
#   Rscript verify_against_handoff.R --tol 1e-6 --out verify_report.csv
#
# Exit status: 0 if every matched file agrees within tolerance, 1 otherwise.

args <- commandArgs(trailingOnly = TRUE)
get_opt <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
TOL        <- as.numeric(get_opt("--tol", "1e-6"))    # absolute tolerance
RTOL       <- as.numeric(get_opt("--rtol", "1e-4"))   # relative tolerance; a value
# passes if it is within TOL (absolute) OR RTOL (relative). RTOL covers large
# optimizer-derived outputs (e.g. Lc_km ~ 8500) whose last digits drift across
# machines; deterministic design tables agree far inside both.
REPORT_OUT <- get_opt("--out", "verify_report.csv")

ANALYSIS_DIR <- normalizePath(".")
OUT_R12C     <- file.path(ANALYSIS_DIR, "out_r12c")
HANDOFF_OUT  <- Sys.getenv("R12C_HANDOFF_OUT",
                           file.path(ANALYSIS_DIR, "handoff_out"))
# handoff/ lives at the repo root; find it by walking up.
find_handoff <- function(start) {
  d <- start
  for (i in 1:8) {
    cand <- file.path(d, "handoff")
    if (dir.exists(cand) && file.exists(file.path(cand, "community_summary.csv")))
      return(normalizePath(cand))
    d <- dirname(d)
  }
  stop("could not locate the repo-root handoff/ directory above ", start)
}
HANDOFF <- find_handoff(ANALYSIS_DIR)

cat("Analysis dir :", ANALYSIS_DIR, "\n")
cat("handoff/     :", HANDOFF, "\n")
cat("out_r12c/    :", if (dir.exists(OUT_R12C)) OUT_R12C else "(missing)", "\n")
cat("handoff_out/ :", if (dir.exists(HANDOFF_OUT)) HANDOFF_OUT else "(missing)", "\n")
cat("tolerance    : abs", TOL, " rel", RTOL, "\n\n")

# ---- build the candidate (produced_path, reference_path, rule) list ----------
ref_files <- list.files(HANDOFF, pattern = "\\.csv$", full.names = FALSE)
prod_out  <- if (dir.exists(OUT_R12C))
  list.files(OUT_R12C, pattern = "\\.csv$", full.names = FALSE) else character(0)
prod_hand <- if (dir.exists(HANDOFF_OUT))
  list.files(HANDOFF_OUT, pattern = "\\.csv$", full.names = FALSE) else character(0)

pairs <- list()
used_ref <- character(0)
add_pair <- function(prod_dir, prod_name, ref_name, rule) {
  if (ref_name %in% used_ref) return(invisible())     # one rule per ref file
  if (!(ref_name %in% ref_files)) return(invisible())
  pairs[[length(pairs) + 1L]] <<- data.frame(
    reference = ref_name,
    produced  = file.path(basename(prod_dir), prod_name),
    rule      = rule,
    prod_path = file.path(prod_dir, prod_name),
    ref_path  = file.path(HANDOFF, ref_name),
    stringsAsFactors = FALSE
  )
  used_ref <<- c(used_ref, ref_name)
}

# rule 1: handoff_out/<f> == handoff/<f>
for (f in prod_hand) add_pair(HANDOFF_OUT, f, f, "exact(handoff_out)")
# rule 2: out_r12c/t_<f> -> handoff/<f>
for (f in prod_out) if (startsWith(f, "t_"))
  add_pair(OUT_R12C, f, sub("^t_", "", f), "t_->name")
# rule 3: out_r12c/t_<f> -> handoff/cmp_<f>
for (f in prod_out) if (startsWith(f, "t_"))
  add_pair(OUT_R12C, f, paste0("cmp_", sub("^t_", "", f)), "t_->cmp_")

pair_df <- if (length(pairs)) do.call(rbind, pairs) else
  data.frame(reference = character(0), produced = character(0),
             rule = character(0), prod_path = character(0),
             ref_path = character(0), stringsAsFactors = FALSE)

# ---- compare one pair --------------------------------------------------------
compare_csv <- function(prod_path, ref_path, tol) {
  a <- tryCatch(read.csv(prod_path, stringsAsFactors = FALSE, check.names = FALSE),
                error = function(e) e)
  b <- tryCatch(read.csv(ref_path, stringsAsFactors = FALSE, check.names = FALSE),
                error = function(e) e)
  if (inherits(a, "error")) return(list(ok = FALSE, note = paste("read produced:", conditionMessage(a))))
  if (inherits(b, "error")) return(list(ok = FALSE, note = paste("read reference:", conditionMessage(b))))
  if (!identical(dim(a), dim(b)))
    return(list(ok = FALSE, note = sprintf("dim %s vs %s",
                                           paste(dim(a), collapse = "x"),
                                           paste(dim(b), collapse = "x"))))
  if (!identical(sort(names(a)), sort(names(b))))
    return(list(ok = FALSE, note = "column names differ"))
  b <- b[, names(a), drop = FALSE]                     # align column order
  # A value passes if it agrees to the ABSOLUTE tol OR the RELATIVE tol. The
  # relative test matters for optimizer-derived quantities (e.g. the spatial
  # decay length Lc_km ~ 8500): a pure 1e-6 absolute tol would flag a
  # last-digit, cross-machine BLAS/optimizer difference as a DIFF even though
  # the reported result (integer crossings, designs) is identical. Deterministic
  # tables (designs, rule checks) still agree far inside both tols.
  max_abs <- 0; bad_cols <- character(0)
  for (nm in names(a)) {
    x <- a[[nm]]; y <- b[[nm]]
    nx <- suppressWarnings(as.numeric(x)); ny <- suppressWarnings(as.numeric(y))
    numeric_col <- !all(is.na(nx)) && !all(is.na(ny))
    if (numeric_col) {
      # Non-finite values (NA, Inf, -Inf) are compared by identity, not
      # arithmetic: abs(Inf - Inf) is NaN, which would otherwise flag an
      # IDENTICAL column as a DIFF. A pair passes the non-finite test when both
      # are the same non-finite value (or both finite and handled below).
      same_nonfinite <- (is.na(nx) & is.na(ny)) |
                        (!is.finite(nx) & !is.finite(ny) & (nx == ny) %in% TRUE) |
                        (is.infinite(nx) & is.infinite(ny) & sign(nx) == sign(ny))
      both_finite <- is.finite(nx) & is.finite(ny)
      d <- abs(nx - ny)
      rel <- d / pmax(abs(nx), abs(ny), .Machine$double.eps)
      # a value FAILS only if both finite AND outside both tolerances, OR the
      # two are non-finite-but-not-the-same, OR one finite and one not.
      fails <- (both_finite & (d > tol) & (rel > RTOL)) |
               (!both_finite & !same_nonfinite)
      if (any(fails, na.rm = TRUE)) {
        bad_cols <- c(bad_cols, nm)
        fd <- d[both_finite]
        if (length(fd)) max_abs <- max(max_abs, max(fd, na.rm = TRUE))
      }
    } else {
      if (!identical(as.character(x), as.character(y)))
        bad_cols <- c(bad_cols, nm)
    }
  }
  if (length(bad_cols))
    return(list(ok = FALSE,
                note = sprintf("cols differ: %s (max abs num diff %.3g)",
                               paste(bad_cols, collapse = ","), max_abs)))
  list(ok = TRUE, note = sprintf("match (%d x %d)", nrow(a), ncol(a)))
}

# ---- run ---------------------------------------------------------------------
results <- data.frame(reference = character(0), produced = character(0),
                      rule = character(0), status = character(0),
                      note = character(0), stringsAsFactors = FALSE)
for (i in seq_len(nrow(pair_df))) {
  r <- compare_csv(pair_df$prod_path[i], pair_df$ref_path[i], TOL)
  results <- rbind(results, data.frame(
    reference = pair_df$reference[i], produced = pair_df$produced[i],
    rule = pair_df$rule[i], status = if (r$ok) "OK" else "DIFF",
    note = r$note, stringsAsFactors = FALSE))
}
results <- results[order(results$status, results$reference), ]

n_ok   <- sum(results$status == "OK")
n_diff <- sum(results$status == "DIFF")
cat(sprintf("Compared %d file pair(s): %d OK, %d DIFF\n\n",
            nrow(results), n_ok, n_diff))
if (n_diff) {
  cat("DIFFERENCES:\n")
  d <- results[results$status == "DIFF", ]
  for (i in seq_len(nrow(d)))
    cat(sprintf("  [DIFF] %-44s  %s\n", d$reference[i], d$note[i]))
  cat("\n")
}

# unmatched produced files (informational)
matched_prod <- basename(pair_df$produced)
all_prod <- c(prod_out, prod_hand)
unmatched_prod <- setdiff(all_prod, matched_prod)
if (length(unmatched_prod)) {
  cat("Produced CSVs with no handoff counterpart (not checked):\n")
  for (f in sort(unmatched_prod)) cat("  -", f, "\n")
  cat("\n")
}
# reference files never matched (informational: expected until full run / scripts 12,30-32)
unmatched_ref <- setdiff(ref_files, used_ref)
if (length(unmatched_ref)) {
  cat(sprintf("handoff/ CSVs not yet produced (%d; expected until full run):\n",
              length(unmatched_ref)))
  for (f in sort(unmatched_ref)) cat("  -", f, "\n")
  cat("\n")
}

write.csv(results, REPORT_OUT, row.names = FALSE)
cat("Wrote", REPORT_OUT, "\n")

quit(status = if (n_diff == 0) 0L else 1L)
