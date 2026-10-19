#' Launch the multiDetectR Shiny application
#'
#' Opens the graphical front end to the allocation rule, a wrapper around
#' [allocate()] for users who prefer not to write R. The app lets you set the
#' number of levels, each level's detectability and cost, the goal, and either a
#' budget or a required detection probability, and it can export a standalone,
#' reproducible R script for the current design.
#'
#' @param ... Passed to [shiny::runApp()] (for example `launch.browser` or
#'   `port`).
#'
#' @return Called for its side effect of starting the app; does not return a
#'   useful value.
#'
#' @details The app lives in `inst/shiny/app.R` and ships with the package.
#'   Requires the \pkg{shiny} and \pkg{DT} packages, which are listed under
#'   Suggests.
#'
#' @examplesIf interactive() && requireNamespace("shiny", quietly = TRUE)
#' run_app()
#' @export
run_app <- function(...) {
  missing <- c("shiny", "DT")[!vapply(c("shiny", "DT"), requireNamespace,
                                      logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("The app needs these packages: ", paste(missing, collapse = ", "),
         ". Install them with install.packages(c(",
         paste(sprintf('\"%s\"', missing), collapse = ", "), ")).",
         call. = FALSE)
  }
  app_dir <- system.file("shiny", package = "multiDetectR")
  if (!nzchar(app_dir) || !file.exists(file.path(app_dir, "app.R"))) {
    stop("Could not find the Shiny app. Try reinstalling multiDetectR.",
         call. = FALSE)
  }
  shiny::runApp(app_dir, ...)
}
