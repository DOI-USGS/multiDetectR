# Freshwater (carp) input data

The input data for this analysis is **not tracked in this package repository**.
The raw U.S. Fish and Wildlife Service eDNA records contain sensitive content
(precise sampling coordinates, site identifiers, and internal identifiers) and
were removed from the public package per the DGEC code review.

## Where the data lives

The data is maintained in a separate data-release repository and will be
published as a citable archive (DOI) once authorized by the data owner (USFWS).
The pipeline reads it from a `data-release` checkout placed alongside this
package by default:

```
<parent>/
  multiDetectR/        <- this package
  data-release/
    input/carp/
      all_species.csv
      Surrogate species matrix_Master.csv
      fws_inputs/
        metabar_2023.csv
        metabar_2024.csv
        00_format_raw_data.R
```

## Pointing the pipeline at the data

By default the carp `_targets.R` reads from `../data-release/input/carp/`
(a sibling of the package root). To read the data from somewhere else (for
example, a downloaded DOI archive), set the `MULTIDETECTR_DATA` environment
variable to the directory that contains `input/carp/`:

```r
Sys.setenv(MULTIDETECTR_DATA = "/path/to/data-release")
targets::tar_make()
```

The expected input is a presence/absence grid; see `get_data()` in
`../R/get_data.R` for the exact columns the reader expects and the reshaping it
performs on the raw release.
