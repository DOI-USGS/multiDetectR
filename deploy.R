#!/usr/bin/env Rscript
## Prepare the multiDetectR Shiny app for deployment to Posit Connect
## (rconnect.chs.usgs.gov).
##
## WHAT THIS DOES
##   1. Installs multiDetectR from GitHub (DOI-USGS/multiDetectR) so the
##      installed package records its GitHub remote + commit SHA. rsconnect
##      reads those fields to pin the deployed app to a reproducible commit.
##   2. Writes inst/shiny/manifest.json via rsconnect::writeManifest(), which
##      resolves the app's library() dependencies (shiny, DT, multiDetectR, and
##      their transitive deps) against the configured repositories.
##
## RUN THIS ON A MACHINE WITH INTERNET ACCESS to github.com and the R package
## repository (CRAN / Posit Package Manager). It cannot run in a network-
## restricted environment.
##
## Usage:
##   Rscript deploy.R            # install from the default branch tip
##   Rscript deploy.R <sha|ref>  # pin to a specific commit/tag/branch
##
## After this writes the manifest, deploy manually (see notes/deploy-connect.md):
##   rsconnect::deployApp("inst/shiny", server = "rconnect.chs.usgs.gov")
## or publish via Connect's Git-backed deployment pointing at inst/shiny.

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
ref <- if (length(args) >= 1) args[1] else NULL   # NULL = default branch (HEAD)

repo_spec <- if (is.null(ref)) {
  "DOI-USGS/multiDetectR"
} else {
  paste0("DOI-USGS/multiDetectR@", ref)
}

## --- 1. Install multiDetectR from GitHub so the remote/SHA are recorded. -----
message("Installing ", repo_spec, " from GitHub ...")
if (requireNamespace("pak", quietly = TRUE)) {
  pak::pak(repo_spec, ask = FALSE)
} else if (requireNamespace("remotes", quietly = TRUE)) {
  remotes::install_github(repo_spec, upgrade = "never")
} else {
  stop("Need either 'pak' or 'remotes' installed to fetch from GitHub. ",
       "install.packages('pak')")
}

## Confirm the GitHub remote was recorded (this is what pins the deploy).
d <- packageDescription("multiDetectR")
if (is.null(d$RemoteType) || !identical(tolower(d$RemoteType), "github")) {
  warning("multiDetectR does not show RemoteType: github after install; the ",
          "manifest may not record a GitHub source. Check the install.")
} else {
  message(sprintf("Recorded source: github %s/%s @ %s",
                  d$RemoteUsername, d$RemoteRepo,
                  substr(d$RemoteSha %||% "", 1, 12)))
}

## --- 2. Write the manifest from the app directory. ---------------------------
if (!requireNamespace("rsconnect", quietly = TRUE)) {
  stop("Need 'rsconnect' installed. install.packages('rsconnect')")
}
app_dir <- "inst/shiny"
manifest_path <- file.path(app_dir, "manifest.json")
message("Writing ", manifest_path, " ...")
rsconnect::writeManifest(appDir = app_dir)

## --- 3. Sanity-check the manifest BEFORE committing/pushing. -----------------
## Connect's Git-backed deploy installs multiDetectR from the source recorded
## here; if it was recorded as a local path, the Connect build will fail.
mf <- jsonlite::fromJSON(manifest_path, simplifyVector = FALSE)
if (is.null(mf$metadata$appmode) || mf$metadata$appmode != "shiny") {
  warning("manifest appmode is not 'shiny' (got: ",
          mf$metadata$appmode %||% "<none>", ").")
}
md <- mf$packages[["multiDetectR"]]
if (is.null(md)) {
  warning("multiDetectR is NOT in the manifest packages. Connect will not ",
          "install it; the app's library(multiDetectR) will fail.")
} else {
  src <- md$Source %||% md$description$Source %||% "<unknown>"
  message("multiDetectR recorded with Source = ", src)
  if (!grepl("github", tolower(src))) {
    warning("multiDetectR Source is '", src, "', not GitHub. Connect cannot ",
            "install it from a local path. Re-run after installing the package ",
            "from GitHub (that is what records the GitHub source).")
  }
}
message("Done. If the checks above are clean, commit and push the manifest:")
message("  git add ", manifest_path)
message("  git commit -m \"Add Connect manifest for the Shiny app\"")
message("  git push origin main")
message("Then point Connect's Git-backed publishing at the 'inst/shiny' dir.")
