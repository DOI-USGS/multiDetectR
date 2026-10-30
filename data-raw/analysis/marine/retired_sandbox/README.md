# Retired marine sandbox scripts (2026-10-29)

These eight loose scripts are the earliest exploratory marine occupancy work, from before the
analysis was formalized into a `targets` pipeline.
They were `git mv`'d here (history preserved) to keep the `marine/` top level clean; the live
`marine/_targets.R` runs only the Pareto cost/richness analysis, and the marine optimization
example now lives in `../../marine_r12c/`.

They are NOT part of any workflow: nothing in `_targets.R`, the `multiDetectR` package, or the
manuscript sources them. They read superseded input CSVs by bare filename (e.g.
`ptconception_occupancy_*level_may.csv`, `model_estimates_*.csv`), several of which are themselves
flagged for removal in `notes/plan-repo-split.md`, and two of them use the older `occstanhm`
package rather than the `occStan` the live analysis uses.

Kept for history, not deletion, in case any detail needs to be consulted.

Files:
- `two_level.R`, `three_level.R` -- early two-/three-level marine occupancy fits via `occstanhm`
  (superseded by the `occStan` `R/` functions, now in `../retired_occstan/`, and ultimately by
  `../../marine_r12c/`).
- `occstan.R`, `single_species_models.R` -- early `occStan` occupancy sandboxes.
- `bootstrap_curves.R` + `helper_files.R` -- an early species-accumulation prototype;
  `helper_files.R` holds a documented, `@export`-ready `species_from_event()` that was once
  considered for the package but is not currently used anywhere. If that accumulation helper is
  ever wanted in `multiDetectR`, lift it from here.
- `plot_single_spp_out.R`, `plot_occStan.R` -- plot the superseded single-event model-estimate
  CSVs.
