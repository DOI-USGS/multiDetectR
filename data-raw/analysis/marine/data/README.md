# Marine input data

The input data for this analysis is **not tracked in this package repository**.
It is maintained in a separate data-release repository and read from a
`data-release` checkout placed alongside this package by default.

## Data sources

The marine analysis uses two inputs:

1. The coastwide Point Conception Dryad release (per-marker reads and samples),
   read in place from `../marine_dryad/dryad_2026-10-13/` at the repository root
   to avoid duplicating the large reads files. This is an external, citable
   release (Altstatt et al.); obtain it from its DOI.
2. The precomputed Lafferty Pareto richness-by-design table,
   `Lafferty_Pareto_Front_data.csv`, used by the cost / Pareto analysis.

The Pareto table lives in the data-release repository:

```
<parent>/
  multiDetectR/        <- this package
  data-release/
    input/marine/
      Lafferty_Pareto_Front_data.csv
```

## Pointing the pipeline at the data

By default the marine `_targets.R` reads the Pareto table from
`../data-release/input/marine/`. Override with the `MULTIDETECTR_DATA`
environment variable to point at the directory that contains `input/marine/`:

```r
Sys.setenv(MULTIDETECTR_DATA = "/path/to/data-release")
targets::tar_make()
```
