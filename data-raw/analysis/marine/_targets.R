library(targets)
library(tarchetypes)
library(here)

# Consolidated marine workflow. Two data sources:
#  1. The coastwide Dryad release (occupancy, accumulation, turnover) -- read in
#     place from the REPO ROOT (parent of the multiDetectR package root that
#     here::here() anchors to), to avoid duplicating ~70 MB of reads files.
#  2. The precomputed Lafferty Pareto richness-by-design table (cost / richness /
#     Pareto), staged in data/. The cost-optimization analysis uses this because
#     the Dryad release's Government Point subset is a single site and cannot
#     reconstruct the X x N x K design grid the cost model needs (see
#     notes/plan-marine-coastwide-pipeline.md, decision A2).
# The input data (Dryad reads + the Lafferty Pareto design table) lives outside
# the package repo. By default we read it from the `data-release` repo checked
# out alongside the package; override with MULTIDETECTR_DATA. The coastwide
# Dryad reads are read in place from the release to avoid duplicating ~70 MB.
# See data-raw/analysis/marine/data/README.md.
release_dir <- Sys.getenv(
  "MULTIDETECTR_DATA",
  unset = file.path(dirname(here::here()), "data-release")
)
dryad_dir <- file.path(dirname(here::here()), "marine_dryad", "dryad_2026-10-13")
data_dir <- file.path(release_dir, "input", "marine")
r_dir <- here::here("data-raw", "analysis", "marine", "R")

samples_path <- file.path(dryad_dir, "samples.csv")
reads_12S_path <- file.path(dryad_dir, "12S_reads.csv")
reads_18S_path <- file.path(dryad_dir, "18S_reads.csv")
reads_COI_path <- file.path(dryad_dir, "COI_reads.csv")
richness_path <- file.path(data_dir, "Lafferty_Pareto_Front_data.csv")

tar_source(files = r_dir)
tar_option_set(packages = c("tibble", "readr", "dplyr", "ggplot2", "here",
  "tidyverse", "rstan", "occStan", "ggthemes", "scales", "multiDetectR",
  "tidyr", "forcats", "ggrepel", "colorspace", "MCMCglmm", "coda"
))

list(
  # input files (tracked for change detection)
  tar_target(name = samples_file, command = samples_path, format = "file"),
  tar_target(name = reads_12S_file, command = reads_12S_path, format = "file"),
  tar_target(name = reads_18S_file, command = reads_18S_path, format = "file"),
  tar_target(name = reads_COI_file, command = reads_COI_path, format = "file"),

  # Government Point, Illumina era: the depth subset (K varies) for occupancy,
  # per marker.
  tar_target(name = gp_12S,
             command = get_dryad_marine(reads_12S_file, samples_file,
                                        platform = "illumina",
                                        site = "Government Point")),
  tar_target(name = gp_18S,
             command = get_dryad_marine(reads_18S_file, samples_file,
                                        platform = "illumina",
                                        site = "Government Point")),
  tar_target(name = gp_COI,
             command = get_dryad_marine(reads_COI_file, samples_file,
                                        platform = "illumina",
                                        site = "Government Point")),

  # Coastwide, Nanopore era: the across-site subset for accumulation and
  # turnover, per marker.
  tar_target(name = cw_12S,
             command = get_dryad_marine(reads_12S_file, samples_file,
                                        platform = "nanopore")),
  tar_target(name = cw_18S,
             command = get_dryad_marine(reads_18S_file, samples_file,
                                        platform = "nanopore")),
  tar_target(name = cw_COI,
             command = get_dryad_marine(reads_COI_file, samples_file,
                                        platform = "nanopore")),

  # Occupancy ---------------------------------------------------------------
  # Government Point (Illumina): two- and three-level per marker (K varies).
  tar_target(name = gp_occ_2_12S, command = fit_cw_2(gp_12S, "12S")),
  tar_target(name = gp_occ_2_18S, command = fit_cw_2(gp_18S, "18S")),
  tar_target(name = gp_occ_2_COI, command = fit_cw_2(gp_COI, "COI")),
  tar_target(name = gp_occ_3_12S, command = fit_cw_3(gp_12S, "12S")),
  tar_target(name = gp_occ_3_18S, command = fit_cw_3(gp_18S, "18S")),
  tar_target(name = gp_occ_3_COI, command = fit_cw_3(gp_COI, "COI")),

  # Coastwide (Nanopore): two-level per marker across the 24 sites.
  tar_target(name = cw_occ_2_12S, command = fit_cw_2(cw_12S, "12S")),
  tar_target(name = cw_occ_2_18S, command = fit_cw_2(cw_18S, "18S")),
  tar_target(name = cw_occ_2_COI, command = fit_cw_2(cw_COI, "COI")),

  # Convergence check across every per-species fit.
  tar_target(name = convergence_cw,
             command = check_convergence_cw(list(
               gp_2_12S = gp_occ_2_12S, gp_2_18S = gp_occ_2_18S,
               gp_2_COI = gp_occ_2_COI,
               gp_3_12S = gp_occ_3_12S, gp_3_18S = gp_occ_3_18S,
               gp_3_COI = gp_occ_3_COI,
               cw_2_12S = cw_occ_2_12S, cw_2_18S = cw_occ_2_18S,
               cw_2_COI = cw_occ_2_COI))),

  # Occupancy figures (Government Point, per marker) ------------------------
  # Precision gain of the three-level over the two-level model (Suppl Fig 0).
  tar_target(name = compare_occ_12S,
             command = compare_occ_levels_cw(gp_occ_2_12S, gp_occ_3_12S, "12S")),
  tar_target(name = compare_occ_18S,
             command = compare_occ_levels_cw(gp_occ_2_18S, gp_occ_3_18S, "18S")),
  tar_target(name = compare_occ_COI,
             command = compare_occ_levels_cw(gp_occ_2_COI, gp_occ_3_COI, "COI")),
  # Distribution of per-species occupancy estimates (psi, theta, p).
  tar_target(name = occ_dist_12S, command = plot_occ_dist_cw(gp_occ_3_12S, "12S")),
  tar_target(name = occ_dist_18S, command = plot_occ_dist_cw(gp_occ_3_18S, "18S")),
  tar_target(name = occ_dist_COI, command = plot_occ_dist_cw(gp_occ_3_COI, "COI")),

  # Accumulation + turnover (coastwide, Nanopore, per marker) ---------------
  tar_target(name = accum_cw_12S, command = accum_curve_cw(cw_12S, "12S")),
  tar_target(name = accum_cw_18S, command = accum_curve_cw(cw_18S, "18S")),
  tar_target(name = accum_cw_COI, command = accum_curve_cw(cw_COI, "COI")),
  tar_target(name = turnover_cw_12S,
             command = turnover_curve_cw(cw_12S, "12S")),
  tar_target(name = turnover_cw_18S,
             command = turnover_curve_cw(cw_18S, "18S")),
  tar_target(name = turnover_cw_COI,
             command = turnover_curve_cw(cw_COI, "COI")),

  # Cost / richness / Pareto (precomputed Lafferty design table; see A2) -----
  tar_target(name = richness_file, command = richness_path, format = "file"),
  tar_target(name = cost_grid, command = make_cost_grid()),
  tar_target(name = richness_model,
             command = fit_richness_powerlaw(richness_file, cost_grid)),
  tar_target(name = pareto_plot,
             command = pareto_front(richness_file, cost_grid))
)
