# _targets.R -- marine r12c sampling-optimization pipeline (ported from Kevin
# Lafferty's handoff/scripts_2026-10-25 / run_r12c.sh, release r12c).
#
# DESIGN OF THE PORT
# Kevin's scripts are run UNCHANGED (only paths were de-hardcoded; see
# README_port.md). They communicate through .rds files written into out_<label>/
# and write t_*.csv result tables, and they are configured entirely by
# environment variables. So each target here runs one script via run_step(),
# which sets the working directory to R/ and the step's env vars, then executes
# the script with Rscript in a fresh process (callr). Dependencies are expressed
# by passing upstream target results into downstream commands (ordering), not by
# targets' file tracking -- the scripts own their file I/O. Each run_step()
# returns a small manifest (script, tag, timestamp, and the files it produced)
# so targets has a real value to track and so a rebuild is visible.
#
# This mirrors run_r12c.sh exactly: same scripts, same env tags, same order.
# The shell driver's machine-specific MEMORY GATE and two-lane scheduling are
# dropped; `targets` handles ordering and (optionally) parallelism.
#
# RUN:  cd multiDetectR/data-raw/analysis/marine_r12c && Rscript -e 'targets::tar_make()'
# The r12c deposit is found via MULTIDETECTR_DATA or the repo-root
# dryad_2026-10-24/ (see common.R). Deterministic (no seeds): the regenerated
# t_*.csv must match handoff/*.csv exactly. ~6 h end to end; longest is
# rule_visit_bands (script 29).

library(targets)

# ---------------------------------------------------------------------------
# Clean the environment the step processes inherit.
# callr::rscript merges its `env` with the PARENT environment, so any of Kevin's
# knobs left exported in the shell that launched tar_make() would leak into every
# child -- and because the scripts read them as Sys.getenv(x, <default>), a stale
# value (even "") silently overrides the intended default. (This caused two
# failures: a leaked COSTTAG="" sent designs to the wrong dir; and SCREEN=""
# flipped the run to "unscreened" -> out_r12c_unscreened.) Mirror
# run_r12c_original.sh's top-of-file `unset ...`: remove these from the parent
# env ONCE here, so callr inherits a clean slate and each script gets its proper
# default. run_step() then sets per-step knobs (incl. COSTTAG) explicitly.
# NOTE: these are UNSET (removed), never set to "" -- "" is not the same as unset
# for Sys.getenv(x, default). CI_DEPOSIT / MULTIDETECTR_DATA are left intact (they
# locate the data deposit).
local({
  leak_vars <- c("COSTTAG", "COSTS", "FITS", "A", "FTYPE", "MARKER", "MK",
                 "GRIDTAG", "RULEVER", "HOPDIR", "BUDGET_AWARE", "TRIGGER",
                 "GRID", "LEFTOVER", "VG", "PC_EVENT", "V6RULE", "SKIP_REFCACHE",
                 "CHECK_DIR", "TEST_B", "SCREEN", "CI_LABEL", "NCORE",
                 "TURNOVER_DEFS_ONLY", "GRIDDIR")
  Sys.unsetenv(leak_vars)
})

# ---------------------------------------------------------------------------
# Paths and the step runner
# ---------------------------------------------------------------------------
ANALYSIS_DIR <- normalizePath(getwd())          # .../marine_r12c (targets runs here)
R_DIR        <- file.path(ANALYSIS_DIR, "R")
OUT_LABEL    <- "r12c"                           # common.R RUN_LABEL when SCREEN=1
OUT_DIR      <- file.path(ANALYSIS_DIR, paste0("out_", OUT_LABEL))

# Locate the r12c deposit explicitly and pass it to every step as CI_DEPOSIT, so
# the run never depends on the launching shell exporting CI_DEPOSIT, nor on
# common.R's getwd()-relative fallback (which resolves to R/ under callr and
# fails -- "cannot open .../R/dryad_2026-10-24/MANIFEST.csv"). Honour a
# CI_DEPOSIT/MULTIDETECTR_DATA the user did set; otherwise walk up from the
# analysis dir to the repo root that holds dryad_2026-10-24/.
DEPOSIT_DIR <- local({
  env_dep <- Sys.getenv("CI_DEPOSIT", "")
  if (nzchar(env_dep) && file.exists(file.path(env_dep, "MANIFEST.csv"))) return(env_dep)
  md <- Sys.getenv("MULTIDETECTR_DATA", "")
  cand <- character(0)
  if (nzchar(md)) cand <- c(file.path(md, "input", "marine", "dryad_2026-10-24"),
                            file.path(md, "dryad_2026-10-24"))
  d <- ANALYSIS_DIR
  for (i in 1:6) { cand <- c(cand, file.path(d, "dryad_2026-10-24")); d <- dirname(d) }
  hit <- cand[file.exists(file.path(cand, "MANIFEST.csv"))]
  if (!length(hit)) stop("could not locate dryad_2026-10-24/MANIFEST.csv; set CI_DEPOSIT")
  normalizePath(hit[1])
})

# ---------------------------------------------------------------------------
# Parallelism (opt-in). Plain tar_make() is SEQUENTIAL: it honours the
# dependency graph but runs one target at a time, so the two independent lanes
# (A = template -> designs -> visit bands; B = disc/capsule + unseen + misc)
# run back to back instead of together. Set the env var R12C_WORKERS to a small
# integer to run independent targets concurrently via a `crew` controller.
#
#   R12C_WORKERS unset or "1"  -> sequential (safe default).
#   R12C_WORKERS=3             -> up to 3 targets at once.
#
# SAFE TO PARALLELISE: every design step writes into its own
# out_r12c/costs_<COSTTAG>/ subdir (the pc_points_*.rds / coast_points_*.rds
# intermediates included), so independent steps touch disjoint files; and
# common.R's wcsv() hard-stops on an existing file, so any accidental overlap
# fails loudly instead of corrupting a table.
#
# DON'T OVERSUBSCRIBE: t_bands (script 29) already uses NCORE=7 cores inside one
# target, and 02/23/24 are memory-hungry. Keep R12C_WORKERS small (2-3) unless
# the machine has many cores and RAM to spare. crew spreads the OTHER targets
# across lane B while 29 runs; it does not speed up 29 itself.
.workers <- suppressWarnings(as.integer(Sys.getenv("R12C_WORKERS", "1")))
if (is.na(.workers) || .workers < 1L) .workers <- 1L

# Cores for 02_fit's internal grid search (NCORE). Default 7 = the widest grid
# (the CG run-multiplier sweep); more cores do not help a single fit. Set
# R12C_FIT_CORES=1 to reproduce the original serial fit exactly.
FIT_CORES <- suppressWarnings(as.integer(Sys.getenv("R12C_FIT_CORES", "7")))
if (is.na(FIT_CORES) || FIT_CORES < 1L) FIT_CORES <- 1L

.controller <- NULL
if (.workers > 1L) {
  if (!requireNamespace("crew", quietly = TRUE))
    stop("R12C_WORKERS=", .workers, " needs the `crew` package. ",
         "install.packages(\"crew\") or unset R12C_WORKERS to run sequentially.")
  .controller <- crew::crew_controller_local(workers = .workers, seconds_idle = 60)
}

tar_option_set(
  packages = c("callr"),
  controller = .controller,   # NULL => sequential tar_make()
  format = "rds"
)

# run_step(): execute one of Kevin's scripts in R_DIR with the given env vars.
# `script`  -- file name in R/ (e.g. "02_fit.R").
# `env`     -- named character vector of extra environment variables for this step.
# `produces`-- paths (relative to the analysis dir) the step is expected to write,
#              used to build the returned manifest and to fail fast if missing.
# `needs`   -- (ignored at runtime) upstream manifests passed in only to force
#              targets to order this step after them.
run_step <- function(script, env = character(0), produces = character(0),
                     needs = NULL, cores = 1L) {
  # SKIP-IF-COMPLETE. Kevin's scripts run in fresh callr processes and their
  # common.R::wcsv() refuses to overwrite an existing CSV (a deliberate guard
  # against silent clobbering). So if a step's outputs are already on disk --
  # e.g. build/fit were run by hand before `tar_make()`, or a prior run got part
  # way -- re-running the script would just error on the first wcsv(). Treat the
  # on-disk outputs as the source of truth: when `produces` is non-empty and
  # every listed file exists, skip the callr run. This makes the pipeline
  # idempotent and lets a by-hand build/fit feed straight into tar_make().
  # To FORCE a rebuild, set R12C_FORCE=1 (or delete the step's outputs). Targets
  # with produces = character(0) (the handoff writers 12/30/31) are never skipped
  # here -- their files are not tracked, so they always run.
  .force <- nzchar(Sys.getenv("R12C_FORCE"))
  if (!.force && length(produces) &&
      all(file.exists(file.path(ANALYSIS_DIR, produces)))) {
    return(data.frame(
      script = script, costtag = NA_character_,
      finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      produced = paste(produces, collapse = "; "),
      status = "skipped (outputs present)",
      stringsAsFactors = FALSE))
  }
  # COSTTAG handling. common.R uses Sys.getenv("COSTTAG", "hours_2026-10-11"),
  # i.e. the default applies only when COSTTAG is UNSET -- an empty-string
  # COSTTAG="" means the TOP of out_r12c/, NOT the default cost subdir. The child
  # process inherits the ambient environment (callr merges env with the parent,
  # and the shell that launched tar_make may have COSTTAG="" exported from an
  # earlier hand-run). So we must set COSTTAG EXPLICITLY for every step, never
  # leave it to chance: default it to "hours_2026-10-11" here (reproducing the
  # driver's reliance on common.R's default), and let a caller override via
  # env -- including COSTTAG="" for 01/02, which legitimately write to the top.
  # (This empty-vs-default distinction was the bug behind "03_designs.R did not
  # produce ...costs_hours_2026-10-11/...".)
  # callr merges `env` with the PARENT environment (verified). The scripts read
  # their knobs with Sys.getenv(x, <default>), where <default> applies ONLY when
  # x is UNSET -- an empty string "" suppresses the default. So we must NOT set
  # these to ""; the leak-prone vars are instead UNSET in the parent once at load
  # time (see the Sys.unsetenv block at the top of this file), so callr inherits
  # a clean environment and each script gets its intended default. Here we set
  # only the few vars a step always needs; COSTTAG is set explicitly (default
  # "hours_2026-10-11") so a stale shell value cannot send outputs to the wrong
  # dir, and a caller overrides it (incl. "" for 01/02).
  base_env <- c(
    MULTIDETECTR_R12C_ROOT = ANALYSIS_DIR,
    CI_DEPOSIT = DEPOSIT_DIR,                     # explicit, resolved at load time
    COSTS = file.path(ANALYSIS_DIR, "costs.csv"),
    COSTTAG = "hours_2026-10-11",
    NCORE = as.character(cores)
  )
  full_env <- modifyList(as.list(base_env), as.list(env))
  # callr runs a clean R session; set wd to R/ so the scripts' relative
  # source("common.R")/source("model.R") resolve.
  callr::rscript(
    script = file.path(R_DIR, script),
    wd = R_DIR,
    env = unlist(full_env),
    show = TRUE,
    fail_on_status = TRUE
  )
  missing <- produces[!file.exists(file.path(ANALYSIS_DIR, produces))]
  if (length(missing)) {
    stop("step ", script, " did not produce: ",
         paste(missing, collapse = ", "))
  }
  data.frame(
    script = script,
    costtag = full_env$COSTTAG %||% "",
    finished = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    produced = paste(produces, collapse = "; "),
    status = "ran",
    stringsAsFactors = FALSE
  )
}
`%||%` <- function(a, b) if (is.null(a) || !length(a) || !nzchar(a)) b else a

# convenience: output paths
o  <- function(...) file.path(paste0("out_", OUT_LABEL), ...)               # top of out_r12c/
oh <- function(...) o("costs_hours_2026-10-11", ...)                        # disc/default cost tag
ot <- function(...) o("costs_hours_2026-10-11_template", ...)              # template cost tag

# ---------------------------------------------------------------------------
# The pipeline (mirrors run_r12c.sh; see README_port.md "Run order")
# ---------------------------------------------------------------------------
list(
  # ===== Stage 1: build + base fits =====
  # COSTTAG="" (empty) for 01/02: their outputs (build.rds, fits.rds, the base
  # t_*.csv) live at the TOP of out_r12c/, matching the driver's `COSTTAG=`.
  tar_target(build,  run_step("01_build.R", env = c(COSTTAG = ""),
                              produces = c(o("build.rds"),
                                           o("t_effort_by_visit.csv"),
                                           o("t_species_detected.csv")))),
  # 02_fit's grid search parallelises over NCORE (see R/02_fit.R). The widest
  # grid is 7 (the CG run-multiplier sweep), so >7 cores does not help one fit;
  # default 7, override with R12C_FIT_CORES. NCORE=1 reproduces the old serial
  # run exactly. Independent of the lane-level R12C_WORKERS (crew) parallelism.
  tar_target(fit,    run_step("02_fit.R", env = c(COSTTAG = ""),
                              cores = FIT_CORES,
                              produces = c(o("fits.rds"), o("t_fit_summary.csv"),
                                           o("t_run_multipliers.csv"),
                                           o("t_level_tests.csv"),
                                           o("t_calibration.csv")),
                              needs = build)),
  tar_target(fit_a0833, run_step("02_fit.R", env = c(COSTTAG = "", A = "0.0833"),
                              cores = FIT_CORES,
                              produces = o("fits_a0.0833.rds"), needs = build)),
  tar_target(pooling,   run_step("16_pooling_penalty.R",
                              produces = o("t_pooling_penalty_summary.csv"), needs = fit)),
  tar_target(a_singles, run_step("18_a_from_singles.R",
                              produces = o("t_a_from_singles_summary.csv"), needs = fit)),

  # ===== Stage 2 lane A: template fit -> template designs -> visit bands =====
  # 23 sets OUT <- FITDIR, so BOTH fits_template.rds and its t_*.csv go to the
  # top of out_r12c/ (not the cost subdir). So both are o().
  tar_target(template, run_step("23_template_fit.R",
                              produces = c(o("fits_template.rds"),
                                           o("t_template_summary.csv")),
                              needs = fit)),
  tar_target(t_designs, run_step("03_designs.R",
                              env = c(FITS = "fits_template.rds", A = "0.03",
                                      COSTTAG = "hours_2026-10-11_template"),
                              produces = ot("t_pc_best_by_budget_disc_a0.03.csv"),
                              needs = template)),
  tar_target(t_rule, run_step("25_rule_designs.R",
                              env = c(FITS = "fits_template.rds",
                                      COSTTAG = "hours_2026-10-11_template"),
                              produces = ot("t_rule_designs_template_v2.csv"),
                              needs = t_designs)),
  tar_target(t_rule_check, run_step("27_rule_check_fine.R",
                              env = c(FITS = "fits_template.rds",
                                      COSTTAG = "hours_2026-10-11_template"),
                              produces = ot("t_rule_check_fine_template_v2.csv"),
                              needs = t_rule)),
  tar_target(t_em, run_step("24_em_convergence.R",
                              env = c(FITS = "fits_template.rds"),
                              produces = o("t_em_convergence_template.csv"),
                              needs = template)),
  # 26 sources 20_turnover.R, which sets OUT <- FITDIR, so 26's t_*.csv land at
  # the TOP of out_r12c/, not the cost subdir -> o(). (Confirmed by the run:
  # c_loc wrote t_location_curve_crossing.csv to the top level.)
  tar_target(t_loc, run_step("26_location_curve.R",
                              env = c(FITS = "fits_template.rds"),
                              produces = o("t_location_curve_crossing_template.csv"),
                              needs = template)),
  tar_target(t_bands, run_step("29_rule_visit_bands.R",
                              env = c(COSTTAG = "hours_2026-10-11_template", NCORE = "7"),
                              cores = 7L,
                              produces = ot("t_rule_visit_bands_all_template.csv"),
                              needs = t_rule)),

  # ===== Stage 2 lane B: disc/capsule designs, unseen profile, misc checks =====
  tar_target(d_designs, run_step("03_designs.R", env = c(FTYPE = "disc"),
                              produces = oh("t_pc_best_by_budget_disc.csv"),
                              needs = fit)),
  tar_target(d_cap33,  run_step("13_extract_cap33.R",
                              produces = oh("t_extract_cap33_designs.csv"), needs = d_designs)),
  tar_target(d_siteknow, run_step("14_site_knowledge.R",
                              produces = oh("t_site_knowledge.csv"), needs = d_designs)),
  tar_target(d_unseen_sens, run_step("07_unseen_sensitivity.R", env = c(A = "0.03"),
                              produces = oh("t_unseen_sensitivity_a0.03.csv"), needs = d_designs)),
  # 04 writes t_common_species_candidates (not t_common_species).
  tar_target(d_common, run_step("04_common_species.R",
                              produces = oh("t_common_species_candidates.csv"), needs = d_designs)),
  tar_target(d_rule,   run_step("25_rule_designs.R",
                              produces = oh("t_rule_designs_v2.csv"), needs = d_designs)),
  tar_target(d_rule_check, run_step("27_rule_check_fine.R",
                              produces = oh("t_rule_check_fine_v2.csv"), needs = d_rule)),
  # unseen profile per marker, then merge the lo/hi fits (as run_r12c.sh does)
  tar_target(d_unseen_12S, run_step("15_unseen_profile.R", env = c(MARKER = "12S"),
                              produces = o("fits_unseen_lo_12S.rds"), needs = d_designs)),
  tar_target(d_unseen_COI, run_step("15_unseen_profile.R", env = c(MARKER = "COI"),
                              produces = o("fits_unseen_lo_COI.rds"), needs = d_designs)),
  tar_target(unseen_merge, {
      needs <- list(d_unseen_12S, d_unseen_COI)                 # force order
      for (e in c("lo", "hi")) {
        x <- c(readRDS(file.path(OUT_DIR, paste0("fits_unseen_", e, "_12S.rds"))),
               readRDS(file.path(OUT_DIR, paste0("fits_unseen_", e, "_COI.rds"))))
        p <- file.path(OUT_DIR, paste0("fits_unseen_", e, ".rds"))
        if (!file.exists(p)) saveRDS(x, p)
      }
      data.frame(script = "unseen_merge", finished = format(Sys.time()), stringsAsFactors = FALSE)
  }),
  tar_target(u03_lo, run_step("03_designs.R",
                              env = c(FITS = "fits_unseen_lo.rds",
                                      COSTTAG = "hours_2026-10-11_unseen_lo",
                                      A = "0.03", FTYPE = "disc"),
                              produces = o("costs_hours_2026-10-11_unseen_lo",
                                           "t_pc_best_by_budget_disc_a0.03.csv"),
                              needs = unseen_merge)),
  tar_target(u03_hi, run_step("03_designs.R",
                              env = c(FITS = "fits_unseen_hi.rds",
                                      COSTTAG = "hours_2026-10-11_unseen_hi",
                                      A = "0.03", FTYPE = "disc"),
                              produces = o("costs_hours_2026-10-11_unseen_hi",
                                           "t_pc_best_by_budget_disc_a0.03.csv"),
                              needs = unseen_merge)),
  tar_target(d_unseen_designs, run_step("17_unseen_profile_designs.R",
                              produces = oh("t_unseen_profile_designs.csv"),
                              needs = list(u03_lo, u03_hi))),
  tar_target(c_designs, run_step("03_designs.R", env = c(A = "0.03", FTYPE = "capsule"),
                              produces = oh("t_pc_best_by_budget_capsule_a0.03.csv"),
                              needs = fit)),
  tar_target(c_em,  run_step("24_em_convergence.R",
                              produces = o("t_em_convergence.csv"), needs = fit)),
  tar_target(c_loc, run_step("26_location_curve.R",
                              produces = o("t_location_curve_crossing.csv"), needs = fit)),
  # 11_enumerate_small_levels: the wide-grid trigger sweep + the base variants
  tar_target(enum_wide, {
      needs <- d_designs
      for (t in c("3", "6", "10")) {
        run_step("11_enumerate_small_levels.R",
                 env = c(BUDGET_AWARE = "1", GRID = "wide", TRIGGER = t))
        run_step("11_enumerate_small_levels.R",
                 env = c(BUDGET_AWARE = "1", GRID = "wide", LEFTOVER = "unspent", TRIGGER = t))
      }
      data.frame(script = "11 wide sweep", finished = format(Sys.time()), stringsAsFactors = FALSE)
  }),
  tar_target(enum_base, {
      needs <- d_designs
      run_step("11_enumerate_small_levels.R")
      run_step("11_enumerate_small_levels.R", env = c(BUDGET_AWARE = "1"))
      run_step("11_enumerate_small_levels.R", env = c(BUDGET_AWARE = "1", TRIGGER = "6"))
      run_step("11_enumerate_small_levels.R", env = c(BUDGET_AWARE = "1", TRIGGER = "10"))
      data.frame(script = "11 base variants", finished = format(Sys.time()), stringsAsFactors = FALSE)
  }),
  # 05 writes three files to the cost subdir (not a single t_export_* table):
  # rule_input_species_types, rule_input_observed_species_posterior_means,
  # rule_input_costs. (The earlier produces name "t_export_rule_inputs.csv" was
  # invented and never produced -- that is why the check failed though 05 ran.)
  tar_target(d_export, run_step("05_export_rule_inputs.R",
                              produces = c(oh("rule_input_species_types.csv"),
                                           oh("rule_input_observed_species_posterior_means.csv"),
                                           oh("rule_input_costs.csv")),
                              needs = d_designs)),
  tar_target(d_branch, run_step("10_branch_test.R",
                              produces = oh("t_branch_test_sites.csv"), needs = d_designs)),
  tar_target(d_rvg,    run_step("09_rule_vs_grid.R",
                              produces = oh("t_rule_vs_grid_single_target_v3.csv"), needs = d_designs)),
  tar_target(a_within, run_step("19_a_within_run_profile.R",
                              produces = o("t_a_within_run_summary.csv"), needs = fit)),
  tar_target(turnover, run_step("20_turnover.R",
                              produces = o("t_turnover_summary.csv"), needs = fit)),
  tar_target(turnover_chk, run_step("21_turnover_checks.R",
                              produces = o("t_turnover_step_checks.csv"), needs = turnover)),
  tar_target(template_sens, run_step("28_template_sensitivity.R",
                              produces = o("t_template_sensitivity.csv"), needs = template)),
  # 22_single_reaction_detections is SCOPED OUT of this pipeline (not a TODO --
  # a decision). It needs CalIntertidal_reads_long_hopflagged_<marker>.rds, which
  # hop_wrapper_r12c.R builds only with (1) code_2026-10-24/ [present],
  # (2) raw_inputs.zip unzipped [present, zipped], and (3) the TaxaFlag package's
  # flag_hopped_detections() [NOT installed; the wrapper hard-stops without it].
  # TaxaFlag is not in this workspace and the wrapper's own note says this was
  # gated on a package reinstall that may never have landed. Per the wrapper, the
  # hop flag is DIAGNOSTIC ONLY -- it "has never filtered a detection or a call in
  # any marker" -- and script 22's three t_single_reaction_* tables are a leaf
  # check, not on the design critical path and not cited in the manuscript. So
  # pulling in a fragile external-package build for it is not worth the reproduction
  # risk. To add it later, see README_port.md "Script 22 (scoped out)".

  # ===== Handoff tables (write Kevin's handoff CSVs; verify the port) =====
  # 12_handoff reads a wide set of lane-B outputs (the cmp_* loop + the
  # rule-vs-grid / enumerate / pc_best / coast_best / site_knowledge / branch /
  # unseen-sensitivity / cap33 tables), so it must run AFTER all of them -- in
  # run_r12c.sh it is dead last. Depend on every lane-B target it consumes, not
  # just fit (the earlier `needs = fit` let targets run it too early, which
  # failed reading t_rule_vs_grid_single_target_v3.csv before d_rvg wrote it).
  tar_target(handoff_main, run_step("12_handoff.R", produces = character(0),
                              needs = list(fit, d_designs, c_designs, d_rvg,
                                           d_branch, d_unseen_sens, d_cap33,
                                           d_siteknow, enum_base))),
  # 30 reads fits_template.rds and the template_* handoff tables that 23 writes,
  # so it needs `template` (not just `fit`).
  tar_target(handoff_species, run_step("30_handoff_species.R",
                              env = c(SKIP_REFCACHE = "1"),
                              produces = character(0), needs = template)),
  # 31 reads handoff_out/template_species_type_weights.csv, which 30 writes, so
  # it must run AFTER handoff_species.
  tar_target(handoff_inverts, run_step("31_detection_histories_inverts.R",
                              produces = character(0), needs = handoff_species)),
  # 32 DATE-STAMPS its output filenames with format(Sys.Date()), so the produced
  # name depends on the run date and cannot be predicted here (Kevin's handoff
  # file is dated 2026-10-25). So produces=character(0) (always runs, like the
  # other handoff writers). Script 32 has its own file-exists guards. NOTE for
  # the verify: the regenerated rule_visits_fixed_* / rule_visit_bands_* carry
  # TODAY's date, so they will NOT name-match Kevin's 2026-10-25 files -- compare
  # them by content, or pass the date, rather than expecting a name match.
  tar_target(handoff_visits, run_step("32_rule_visits_handoff.R",
                              produces = character(0),
                              needs = t_bands))
)
