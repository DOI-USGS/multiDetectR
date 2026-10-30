# Retired occStan marine functions (2026-10-27)

These eight functions were the occStan occupancy analysis of the marine (Point Conception)
eDNA data: per-species two- and three-level occupancy fits, the two/three-level precision
comparison, the per-species occupancy distributions, the coastwide accumulation and turnover
curves, and the Dryad reader that fed them.

They were retired from the marine workflow on 2026-10-27 under decision 1 in
`notes/scoping-marine-example-r12c.md`: occStan handles at most three levels, but the marine
survey design is deeper, so the marine optimization example was moved to Kevin Lafferty's NPMLE
rule-v2 method, ported in `../../marine_r12c/`.
The occupancy material and its figures were also removed from the manuscript.

They are kept here (not deleted) for history and in case any detail needs to be consulted.
They are NOT sourced by `../_targets.R` (which now runs only the Pareto cost/richness analysis),
and they depend on packages the trimmed workflow no longer loads (`occStan`, `rstan`).

occStan itself is still used, unchanged, for the FRESHWATER example in `../../carp/`.

Files:
- `get_dryad_marine.R` -- reader: Dryad `<marker>_reads.csv` + `samples.csv` -> occupancy frame.
- `fit_cw_2.R`, `fit_cw_3.R` -- two- and three-level per-species occupancy fits.
- `compare_occ_levels_cw.R` -- precision gain, three-level vs two-level.
- `plot_occ_dist_cw.R` -- per-species occupancy (psi/theta/p) distributions.
- `check_convergence_cw.R` -- convergence across the per-species fits.
- `accum_curve_cw.R`, `turnover_curve_cw.R` -- coastwide accumulation and turnover curves.
