# multiDetectR



#### Authors: Kevin D. Lafferty, Danelle M. Larson, Stephen F. Spear, Nick Schulte, Jason A. Coombs, Patrick DeHaan, Cari-Ann Hayer, Aaron P. Maloy, Richard A. Erickson
#### Point of contact: Richard A. Erickson (rerickson@usgs.gov)
#### Repository Type: _R package_
#### Year of Origin:   2026
#### Year of Version:  2026
#### Version: 0.1.0
#### Digital Object Identifier (DOI): https://doi.org/10.5066/XXXXXXXXXX
#### USGS Information Product Data System (IPDS) no.: IP-XXXXXXX (internal agency tracking)

***

_Suggested Citation:_

Lafferty, KD; Larson, DM; Spear, SF; Schulte, N; Coombs, JA; DeHaan, P; Hayer, C-A; Maloy, AP; Erickson, RA.
2026.
multiDetectR.
Version 0.1.0.
U.S. Geological Survey software release. Reston, Va.
https://doi.org/10.5066/XXXXXXXX


_Authors' [ORCID](https://orcid.org) nos.:_

- Kevin D. Lafferty [0000-0001-7583-4593](https://orcid.org/0000-0001-7583-4593)
- Danelle M. Larson [0000-0001-6349-6267](https://orcid.org/0000-0001-6349-6267)
- Stephen F. Spear [0000-0001-8351-9382](https://orcid.org/0000-0001-8351-9382)
- Richard A. Erickson [0000-0003-4649-482X](https://orcid.org/0000-0003-4649-482X)


***
***

GitLab Pipeline status:
[![pipeline status](https://code.usgs.gov/umesc/quant-ecology/multidetectr/badges/main/pipeline.svg)](https://code.usgs.gov/umesc/quant-ecology/multidetectr/-/pipelines)

***
***

Exact allocation of detection effort across the levels of a nested survey, at
any depth.

A survey that looks for a target (a species, an infectious agent, a
contaminant) often nests its sampling: sites hold samples, samples hold
replicates, and so on. Each level has its own cost and its own detectability,
and a fixed budget must be split among them. `multiDetectR` implements the
closed-form rule that solves this allocation exactly, for any number of levels,
and reproduces the paper's own numerical checks as unit tests.

The package models a nested survey as a cascade of conditional Bernoulli trials
and measures detection on a log (nat) scale, so the design problem separates
into a chain of one-variable problems solved from the innermost level outward.
Function and argument names mirror the equations of the accompanying
manuscript.

## Installation

```r
if (!require("remotes")) install.packages("remotes")
remotes::install_gitlab("umesc/quant-ecology/multidetectr", host = "code.usgs.gov")
```

## Quick start

```r
library(multiDetectR)

## Allocate a budget across a three-level eDNA survey
## (sites -> samples -> PCR replicates).
des <- allocate(
  theta  = c(0.3, 0.5, 0.6),   # per-level conditional detectabilities
  cost   = c(200, 20, 2),      # per-level unit costs
  goal   = "detection",
  budget = 5000
)
des
summary(des)
plot(des)                      # the budget-detection frontier
```

Supply `P_target` instead of `budget` to find the cheapest design that reaches a
required detection probability. Use `goal = "count"` for the count goal,
`enumerate_designs()` for the exact integer optimum at small budgets,
`optimality_identity()` to audit a costed design, and `allocate_branches()` to
choose among methods or places.

See `vignette("multiDetectR")` for a walkthrough.

## What the functions compute

| Function | What it does | Manuscript |
|---|---|---|
| `detection_recursion()` | detection probability of a nested design | eqs 1-2 |
| `total_cost()`, `unit_cost_full()` | cost of a design | eqs 3, 8 |
| `nats()`, `prob()` | detection <-> hazard (nats) | log scale |
| `solve_level()` | optimal effort at one level | eqs 14, 15 |
| `allocate()` | the exact innermost-outward pass | eqs 4-6, 10-13 |
| `enumerate_designs()` | exact integer optimum | Supplement S2 |
| `optimality_identity()` | audit a design (A = B) | Supplement S1.8 |
| `allocate_branches()` | concentrate vs spread | Supplement S1.5 |

## Repository organization

This repository contains the code for an R package and the standard R repository files
(see the official R documentation [_Writing R Extensions_][r_ext] or [_R Packages_][r_pkgs]).
In addition to the R package source files, the main level of the repository contains:

- `README.md`, which is this file.
- `.gitignore`, a Git ignore file for the repository.
- `.Rbuildignore`, the list of files ignored when building the R package.
- `LICENSE.md`, the license for the software (public domain / Unlicense).
- `Contributing.md`, which describes how to contribute to this project.
- `DISCLAIMER.md`, the standard USGS disclaimer.
- `code.json`, the required machine-readable metadata about the code.

The reproducible analyses that feed the manuscript, the vignettes, and the Shiny app live under `data-raw/analysis/`.
That directory is excluded from the built package by `.Rbuildignore`, so the installable package stands alone.

## Disclaimer

See `DISCLAIMER.md`. This software is preliminary or provisional and is subject
to revision.

## License

This project is in the public domain in the United States and is additionally
released worldwide under the Unlicense. See `LICENSE.md`.

[r_ext]: https://cran.r-project.org/doc/manuals/r-release/R-exts.html
[r_pkgs]: https://r-pkgs.org/
