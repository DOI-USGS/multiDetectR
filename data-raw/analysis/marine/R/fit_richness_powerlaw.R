fit_richness_powerlaw <- function(richness_file, cost_grid, seed = 20260101) {
  # Write the regenerated parameter table to the data-release repo (output/),
  # not the package repo. Mirrors the MULTIDETECTR_DATA convention in _targets.R.
  release_dir <- Sys.getenv(
    "MULTIDETECTR_DATA",
    unset = file.path(dirname(here::here()), "data-release")
  )
  data_dir_out <- file.path(release_dir, "output", "marine")
  dir.create(data_dir_out, recursive = TRUE, showWarnings = FALSE)
  txt_dir <- here::here("data-raw", "analysis", "marine", "outputs")

  X_sites <- N_samples <- K_subsamples <- richness <- #nolint
    total_cost <- effort <- NULL #nolint

  # Bayesian power-law richness model (ported from cost_benefit_analysis_DL.R):
  #   richness = a * effort^b * total_cost^c
  # Fit on the log scale as log(richness) = log(a) + b log(effort) + c log(cost)
  # with MCMCglmm, so a = exp(Intercept), b = logeffort, c = logcost.

  richness_data <- readr::read_csv(richness_file, col_types = readr::cols())

  # One row per design: mean richness across bootstrap replicates, joined to the
  # design's total cost, with effort = X * N * K.
  rich_by_design <-
    richness_data |>
    dplyr::group_by(X_sites, N_samples, K_subsamples) |>
    dplyr::summarize(richness = mean(richness, na.rm = TRUE),
                     .groups = "drop")

  dat <-
    cost_grid |>
    dplyr::select(X_sites, N_samples, K_subsamples, total_cost) |>
    dplyr::inner_join(rich_by_design,
                      by = c("X_sites", "N_samples", "K_subsamples")) |>
    dplyr::mutate(effort = X_sites * N_samples * K_subsamples) |>
    dplyr::filter(richness > 0, effort > 0, total_cost > 0)

  dat$logrich <- log(dat$richness)
  dat$logeffort <- log(dat$effort)
  dat$logcost <- log(dat$total_cost)

  # Weakly informative priors: intercept sd 5, slopes sd 2; residual inverse-Wishart.
  prior <- list(
    B = list(mu = c(0, 0, 0), V = diag(c(25, 4, 4))),
    R = list(V = 1, nu = 0.002)
  )

  set.seed(seed)
  fit <- MCMCglmm::MCMCglmm(
    logrich ~ logeffort + logcost,
    data = dat,
    family = "gaussian",
    prior = prior,
    nitt = 60000, burnin = 10000, thin = 25,
    verbose = FALSE
  )

  sol <- as.data.frame(fit$Sol)
  draws <- list(
    `a (scaling)` = exp(sol$`(Intercept)`),
    `b (effort exponent)` = sol$logeffort,
    `c (cost exponent)` = sol$logcost
  )

  param_table <-
    tibble::tibble(
      parameter = names(draws),
      median = sapply(draws, stats::median),
      CrI_2.5 = sapply(draws, stats::quantile, 0.025),
      CrI_97.5 = sapply(draws, stats::quantile, 0.975),
      P_gt_0 = sapply(draws, function(x) mean(x > 0))
    ) |>
    dplyr::mutate(dplyr::across(dplyr::where(is.numeric),
                                function(x) round(x, 3)))

  # Bayesian R^2 on the log(richness) scale, per posterior draw.
  x_mat <- stats::model.matrix(~ logeffort + logcost, data = dat)
  fitted_draws <- x_mat %*% t(as.matrix(fit$Sol))
  var_fit <- apply(fitted_draws, 2, stats::var)
  var_res <- as.numeric(fit$VCV[, "units"])
  r2_draws <- var_fit / (var_fit + var_res)

  # Factor by which effort and cost must scale together to double richness:
  # f^(b+c) = 2, so f = 2^(1/(b+c)), per draw.
  f_double <- 2 ^ (1 / (sol$logeffort + sol$logcost))

  param_table <-
    dplyr::bind_rows(
      param_table,
      tibble::tibble(parameter = "Bayesian R^2",
                     median = round(stats::median(r2_draws), 3),
                     CrI_2.5 = round(stats::quantile(r2_draws, 0.025), 3),
                     CrI_97.5 = round(stats::quantile(r2_draws, 0.975), 3),
                     P_gt_0 = NA_real_),
      tibble::tibble(parameter = "factor to double richness",
                     median = round(stats::median(f_double), 3),
                     CrI_2.5 = round(stats::quantile(f_double, 0.025), 3),
                     CrI_97.5 = round(stats::quantile(f_double, 0.975), 3),
                     P_gt_0 = NA_real_)
    )

  readr::write_csv(x = param_table,
                   file = file.path(data_dir_out,
                                    "richness_model_parameters.csv"))

  # Convergence diagnostics for the single MCMCglmm chain: effective sample
  # size per fixed effect, and a Heidelberger-Welch stationarity test.
  eff <- coda::effectiveSize(fit$Sol)
  heidel <- coda::heidel.diag(fit$Sol)
  sink(file.path(txt_dir, "richness_model_convergence.txt"))
  cat("MCMCglmm power-law richness model convergence\n")
  cat("Stored posterior draws:", nrow(sol), "\n\n")
  cat("Effective sample size (per fixed effect):\n")
  print(round(eff, 1))
  cat("\nHeidelberger-Welch stationarity test:\n")
  print(heidel)
  sink()

  param_table
}
