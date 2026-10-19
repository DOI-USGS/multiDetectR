# multiDetectR

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

From source (development):

```r
# install.packages("remotes")
remotes::install_git("https://code.usgs.gov/umesc/quant-ecology/multidetectr")
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

## Disclaimer

See `DISCLAIMER.md`. This software is preliminary or provisional and is subject
to revision.

## License

CC0 1.0 Universal (public domain dedication). See `LICENSE`.
