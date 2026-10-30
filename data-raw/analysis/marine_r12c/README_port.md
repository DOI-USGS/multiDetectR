# marine_r12c -- port of Kevin's r12c sampling-optimization analysis

This directory ports Kevin Lafferty's Pt. Conception sampling-optimization analysis
(`handoff/scripts_2026-10-25/`, release r12c) into the multiDetectR analysis tree. It is the
marine OPTIMIZATION example for the manuscript, replacing the dropped occStan marine analysis
(scope decision: `notes/scoping-marine-example-r12c.md`).

## Status (2026-10-26)

Port is mostly built; awaiting the full run + acceptance diff. Completed:
- `R/` holds Kevin's 32 numbered scripts + `model.R` + `common.R` + `hop_wrapper_r12c.R`,
  copied verbatim EXCEPT for de-hardcoded paths (below). `costs.csv` at the dir root.
- Paths de-hardcoded (no `/Users/lafferty/...` remains); all 34 scripts parse.
- **Verified:** ported `01_build.R` runs on the real r12c deposit and reproduces the brief's
  frame exactly -- 12S intertidal species = 61, COI = 265 -- matching Kevin's
  `community_summary.csv` (`species_observed`). Outputs land in `out_r12c/`.
- **`_targets.R` written and validated:** `tar_manifest()` lists 41 targets and
  `tar_network()` resolves 42 dependency edges from the `needs =` arguments, so ordering is
  correct; `callr` is installed. It mirrors `run_r12c_original.sh` order. Script 22 (hop-flag)
  is intentionally excluded until its `hop_wrapper_r12c.R` input is wired (one leaf result,
  off the design critical path).
- **`verify_against_handoff.R` written** -- the acceptance test (see "Verification plan").
- `02_fit.R` being run by hand (~90 min) to confirm the base fits reproduce
  (12S ll = -9686.57, COI ll = -15458.73) before the full run.

Not yet done:
- Full `tar_make()` run + the acceptance diff against the handoff CSVs.
- Wire script 22's hop-flag input (or leave it documented as out of scope).
- Retire the occStan marine targets (decision 1) after this pipeline is the marine source.

## Package dependencies

Needed on the analysis machine (beyond base R):
- `targets`, `callr` -- the pipeline driver.
- `dplyr` -- used throughout Kevin's scripts.
- `parallel` -- base R; `02_fit.R` grid search and `29_rule_visit_bands.R` use `mclapply`.
- `sf` + `rnaturalearth` + **`rnaturalearthhires`** -- `20_turnover.R` (and, via it,
  `26_location_curve.R`) use the Natural Earth 1:10m coastline (`ne_coastline(scale = 10)`) to
  compute along-coast distances. `rnaturalearthhires` is a data package and was NOT installed by
  default on this box; install it from `ropensci/rnaturalearthhires`
  (`remotes::install_github("ropensci/rnaturalearthhires")`) or the CRAN tarball. On an offline /
  locked-down machine, obtain the tarball out-of-band and `R CMD INSTALL` it. See
  `notes/r12c-rnaturalearthhires-blocker.md`.
- `crew` -- ONLY for the optional `R12C_WORKERS` parallelism; not required for a sequential run.

## What was changed from Kevin's originals (ONLY paths; numerics untouched)

- `common.R`:
  - `ROOT` now comes from `MULTIDETECTR_R12C_ROOT` (the driver sets it to this dir), falling
    back to `getwd()`. Keeps the trailing slash and all original roles.
  - `DEPOSIT` default now resolves the r12c release via the `MULTIDETECTR_DATA` convention,
    then a repo-root `dryad_2026-10-24/` fallback. `CI_DEPOSIT` still overrides (as Kevin's
    code intended). `SCREEN=1` (screened calls) kept as the default; the MANIFEST md5 check
    kept intact.
  - `CI_ARCHIVE` fallback default made local (only used if `CI_DEPOSIT` is empty; not in this
    workflow).
- All 31 scripts: the hardcoded `source(".../frame2/common.R")` (and `; source(paste0(ROOT,
  "model.R"))`) rewritten to relative `source("common.R")` / `source("model.R")`; the two
  turnover scripts' `source(".../20_turnover.R")` to relative `source("20_turnover.R")`.
- Leaf handoff scripts' Mac output/cache paths made env-overridable with local defaults:
  `12_handoff.R`, `30_handoff_species.R`, `31_detection_histories_inverts.R` use
  `R12C_HANDOFF_OUT` (default `ROOT/handoff_out/`) and `R12C_REFCACHE`; `22_single_reaction`'s
  `HOPDIR` default made local.

## How scripts are run (the contract the _targets.R must honor)

- Working directory = `R/` (so the relative `source("common.R")`/`source("model.R")` resolve).
- Env vars the driver sets: `MULTIDETECTR_R12C_ROOT` (this dir), `CI_DEPOSIT` (the r12c deposit)
  or rely on `MULTIDETECTR_DATA`, `COSTS` (= `ROOT/costs.csv`), plus per-step tags exactly as in
  `run_r12c_original.sh`: `COSTTAG`, `FITS`, `A`, `FTYPE`, `GRIDTAG`, `MARKER`, `BUDGET_AWARE`,
  `GRID`, `TRIGGER`, `LEFTOVER`, `VG`, `NCORE`, `RULEVER`, `HOPDIR`.
- Scripts hand off via `.rds` in `out_r12c/` (`build.rds`, `fits.rds`, `fits_template.rds`,
  `fits_unseen_{lo,hi}.rds`, `pc_points_*.rds`, `coast_points_*.rds`) and write `t_*.csv` tables.
- NOTE / open structural item: `common.R` sets `COSTS_FILE <- paste0(ROOT, "../costs.csv")` and
  `model.R` is expected relative to the scripts. In this port `costs.csv` is at `ROOT/costs.csv`
  (not `ROOT/../costs.csv`), so the driver sets `COSTS` explicitly (which `common.R` honors via
  `Sys.getenv("COSTS", ...)`). The `01_build.R` test used `COSTS=$ROOT/costs.csv`. Confirm every
  design step is likewise given `COSTS`, or move `costs.csv` to `ROOT/../` to match the default.
  A copy of `model.R` was also placed at `ROOT/model.R` defensively; the scripts source it
  relatively, so that copy is belt-and-suspenders -- the `_targets.R` should standardize on one.

## Run order (from run_r12c_original.sh)

Stage 1 (fits): `01_build` (COSTTAG=) -> `02_fit` (+ `02_fit` A=0.0833) ; `16`, `18` in parallel.
Stage 2 lane A (template): `23_template_fit` -> `03_designs`(FITS=fits_template) -> `25_rule_designs`
  -> `27_rule_check_fine` -> `24_em_convergence` -> `26_location_curve` -> `29_rule_visit_bands`.
Stage 2 lane B (disc/capsule + unseen + misc): `03`(disc) -> `13`,`14`,`07`,`04`,`25`,`27`,
  `15`(12S),`15`(COI) -> unseen lo/hi merge -> `u03` lo/hi -> `17` -> `03`(capsule) -> `24`,`26`
  -> `11` wide trigger {3,6,10} x {fill,unspent} -> `05`,`10`,`09`,`11` variants -> `19`,`20`,`21`,
  `28` -> `22`(HOPDIR). Then `12_handoff`, `30_handoff_species`(SKIP_REFCACHE=1), `31`, `32`.
The machine-specific MEMORY GATE in the shell driver is dropped; targets manages dependencies.
Deterministic (no seeds) -> a correct port reproduces the handoff CSVs exactly.

## Script 22 (scoped out)

`22_single_reaction_detections.R` is intentionally NOT in `_targets.R`. It is a leaf diagnostic
(j = 1 single-reaction detections: are they low-read, blank- or index-hop-explained?) that writes
`t_single_reaction_detections_{rows,summary,by_species}_v2.csv`. It is not on the design critical
path and none of its outputs are cited in the manuscript.

It reads `CalIntertidal_reads_long_hopflagged_<marker>.rds` from `HOPDIR` (default `ROOT/hop_out/`).
`hop_wrapper_r12c.R` builds those files, but only with all three of:

1. `code_2026-10-24/` -- the Dryad data-production code (`ptcon_ingest.R`, `paths.R`,
   `library_drops.R`, `pc_build_inputs()`). PRESENT in the repo root.
2. The deposit's `raw_inputs.zip` unzipped (needs `derived/` crosswalks and `decisions/` CSVs).
   PRESENT but still zipped at `dryad_2026-10-24/raw_inputs.zip`.
3. The `TaxaFlag` package exporting `flag_hopped_detections()`. NOT INSTALLED -- and the wrapper
   hard-`stop()`s without it. Its own note says this was gated on a package-review reinstall that
   may never have landed.

Per the wrapper's header, the hop flag is DIAGNOSTIC ONLY -- it "has never filtered a detection or
a call in any marker" (`DROP_FLAGS` is empty; validity is carried as provenance). So wiring a
fragile external-package build into the pipeline for an uncited leaf table is not worth the
reproduction risk. Decision recorded in `notes/scoping-marine-example-r12c.md`.

To enable it later: install `TaxaFlag` (commit 69feacd per the wrapper), unzip `raw_inputs.zip`,
confirm `hop_wrapper_r12c.R`'s hardcoded scratch paths are repointed at the unzipped `raw_inputs/`
and `code_2026-10-24/.../Multi_site_multi_marker` (it still carries Kevin's `/private/tmp/...`
paths -- de-hardcode those first), run it once per marker (`MK=12S`, `COI`) to write
`hop_out/CalIntertidal_reads_long_hopflagged_<marker>.rds`, then add a `tar_target` running
`22_single_reaction_detections.R` with `HOPDIR` pointing at `hop_out/`, downstream of `fit`.

## Run and verify

From this directory (`marine_r12c/`):

```sh
Rscript -e 'targets::tar_make()'      # SEQUENTIAL; ~6 h; longest job is 29 (~2.5 h on 7 cores)
Rscript verify_against_handoff.R      # acceptance diff; exits 0 if all matched files agree
```

Hand `tar_make()` to Richie -- the WSL terminal is flaky for long runs, and it returns exit -1
even when a command succeeded (read the output files, not the exit code).

### Running it faster (optional parallelism)

Plain `tar_make()` is sequential: it respects the dependency graph but runs one target at a
time, so the two independent lanes (A = template -> designs -> visit bands; B = disc/capsule +
unseen + misc) run back to back. To overlap them, set `R12C_WORKERS` and `targets` uses a `crew`
controller:

```sh
R12C_WORKERS=3 Rscript -e 'targets::tar_make()'   # up to 3 targets at once (needs `crew`)
```

This spreads lane B's many independent leaf steps across the window while lane A's long
`29_rule_visit_bands` runs, cutting wall-clock toward the critical path
(`01_build -> 02_fit -> 23_template -> 03_designs -> 25_rule -> 29_bands`), roughly ~3-4 h.

Notes:
- It is SAFE: each design step writes into its own `out_r12c/costs_<COSTTAG>/` subdir (the
  `pc_points_*.rds` / `coast_points_*.rds` intermediates included), so independent steps touch
  disjoint files; and `common.R`'s `wcsv()` hard-stops on an existing file, so any accidental
  overlap fails loudly rather than corrupting a table. The numbers are identical either way
  (deterministic, no seeds).
- Keep `R12C_WORKERS` small (2-3). `29_rule_visit_bands` already uses 7 cores inside one target,
  and 02/23/24 are memory-hungry -- too many workers oversubscribes CPU/RAM. `crew` does not
  speed up script 29 itself, only the other targets running alongside it.
- `crew` must be installed (`install.packages("crew")`); otherwise `_targets.R` stops with that
  instruction. Leaving `R12C_WORKERS` unset runs sequentially with no extra dependency.

### `02_fit.R` is now multi-core (the biggest serial bottleneck)

`02_fit.R` was the single longest serial step (~3 h observed), because its capsule (KC, 5 values)
and per-batch run-multiplier (CG, 7 values) grid searches call `fit_model()` ~75 times per sweep
-- 10 non-ref batches x 7 + 5 -- x2 sweeps x2 markers, i.e. ~300 EM fits, all on one core.

Those grid evaluations are independent (each is an `llf()` call, then `which.max`), so the port
now runs each grid over `NCORE` fork workers (`parallel::mclapply`) instead of `sapply`. This is
numerically identical: same grid, same values, same `which.max` selection, same coordinate-ascent
order (the batch loop stays sequential so each `cv[b]` still feeds the next batch). Verified that
`mclapply` returns results in input order even with `mc.preschedule = FALSE`, so the selection
matches `sapply` exactly.

`_targets.R` passes `R12C_FIT_CORES` (default 7 = the widest grid; more does not help one fit) to
the two `02_fit` targets. Set `R12C_FIT_CORES=1` to reproduce the original serial fit exactly.
This should cut `02_fit` from ~3 h toward well under an hour. It is independent of the lane-level
`R12C_WORKERS`: `R12C_FIT_CORES` speeds up each fit internally; `R12C_WORKERS` overlaps whole
lanes. On a many-core box you can use both (mind total cores: the two `02_fit` targets could run
on separate crew workers, each using `R12C_FIT_CORES` cores).

`verify_against_handoff.R` is the acceptance test. Kevin's analysis is deterministic (no seeds),
so a correct port reproduces his result CSVs value-for-value. The script pairs each produced CSV
with its handoff/ counterpart -- by exact name for the handoff CSVs written into `handoff_out/`
(scripts 12/30/31/32), and by dropping the `t_` prefix / mapping `t_` -> `cmp_` for the summary
tables in `out_r12c/` -- then compares numeric columns within `--tol` (default 1e-6) and
character columns exactly. It prints per-file OK/DIFF, lists anything matched on only one side,
writes `verify_report.csv`, and exits non-zero if any matched pair differs. The key acceptance
files are `rule_designs_template_v2.csv`, `template_summary.csv`,
`location_curve_crossing_template.csv`,
`rule_visits_fixed_12S_COI_template_2026-10-25_r12c_V6low.csv`, and
`rule_check_fine_template_v2.csv`. `01_build` already matches (61 / 265).

Expected "not produced" entries (NOT failures -- do not chase these):
- `abalone_limpet_COI_reference.csv` -- script 30 is run with `SKIP_REFCACHE=1`; this one table
  reads a reference-sequence cache (`cache_reference/`) that is not in the deposit and does not
  depend on the r12c release (Kevin carries it forward unchanged). Everything else in script 30
  is produced.
- `single_reaction_detections_{rows,summary,by_species}_v2.csv` -- produced only by script 22,
  which is scoped out (see "Script 22 (scoped out)").
- the handoff writers (12/30/31) write into `handoff_out/` with Kevin's exact filenames; the
  verify matches those by exact name. They use an md5 guard (identical file = left alone, only a
  DIFFERING file stops them), so they are safe to re-run.
