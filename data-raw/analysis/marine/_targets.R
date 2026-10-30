library(targets)
library(tarchetypes)
library(here)

# Marine workflow -- Pareto cost/richness analysis ONLY.
#
# SCOPE CHANGE 2026-10-27 (notes/scoping-marine-example-r12c.md, decision 1):
# the marine OCCUPANCY example no longer uses occStan. occStan caps at three
# levels, but the marine design is deeper, so the marine optimization example was
# moved to Kevin Lafferty's NPMLE rule-v2 method, ported in
# ../marine_r12c/ (its own targets pipeline on the r12c Dryad release). The
# occStan occupancy fits, the two/three-level precision comparison, the
# per-species occupancy distributions, and the coastwide accumulation/turnover
# curves were REMOVED from this workflow and from the manuscript. occStan is kept
# for the FRESHWATER example (../carp/), so the paper still shows the method with
# two modelling approaches.
#
# What remains here is the Pareto front (decision 2): the cost / richness /
# Pareto-frontier visualization from the precomputed Lafferty design table. It
# does not use occStan or the Dryad reads, so it stands alone.
#
# The retired functions were moved to retired_occstan/ (not deleted) for history;
# they are no longer sourced. If the marine occupancy analysis is ever needed
# again, see ../marine_r12c/ (the current marine method) rather than reviving these.

release_dir <- Sys.getenv(
  "MULTIDETECTR_DATA",
  unset = file.path(dirname(here::here()), "data-release")
)
data_dir <- file.path(release_dir, "input", "marine")
r_dir <- here::here("data-raw", "analysis", "marine", "R")

richness_path <- file.path(data_dir, "Lafferty_Pareto_Front_data.csv")

tar_source(files = r_dir)
# Packages the Pareto functions use (via ::): MCMCglmm + coda for the power-law
# richness model, colorspace for the frontier palette, ggplot2/ggrepel for the
# plot, readr/dplyr/tibble for data. The occStan-only deps (rstan, occStan) were
# dropped with the occupancy targets.
tar_option_set(packages = c("tibble", "readr", "dplyr", "ggplot2", "here",
  "scales", "multiDetectR", "tidyr", "ggrepel", "MCMCglmm", "coda", "colorspace"
))

list(
  # Cost / richness / Pareto (precomputed Lafferty design table).
  tar_target(name = richness_file, command = richness_path, format = "file"),
  tar_target(name = cost_grid, command = make_cost_grid()),
  tar_target(name = richness_model,
             command = fit_richness_powerlaw(richness_file, cost_grid)),
  tar_target(name = pareto_plot,
             command = pareto_front(richness_file, cost_grid))
)
