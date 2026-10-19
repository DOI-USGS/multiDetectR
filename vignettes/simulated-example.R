## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 4)
set.seed(2026)

## ----setup--------------------------------------------------------------------
library(multiDetectR)

## ----truth--------------------------------------------------------------------
truth <- list(
  theta = c(site = 0.35, sample = 0.55, replicate = 0.70),
  cost  = c(site = 300,  sample = 25,   replicate = 3)
)
truth

## ----simulator----------------------------------------------------------------
simulate_survey <- function(theta, n_sites, n_samples, n_reps) {
  rows <- list()
  for (s in seq_len(n_sites)) {
    site_has <- stats::runif(1) < theta[1]
    for (m in seq_len(n_samples)) {
      sample_has <- site_has && stats::runif(1) < theta[2]
      detections <- if (sample_has) {
        stats::rbinom(1, size = n_reps, prob = theta[3])
      } else {
        0L
      }
      rows[[length(rows) + 1]] <- data.frame(
        site = s, sample = m,
        n_reps = n_reps, detections = detections
      )
    }
  }
  do.call(rbind, rows)
}

## A generous pilot survey to estimate the per-level detectabilities.
pilot <- simulate_survey(truth$theta, n_sites = 60, n_samples = 8, n_reps = 10)
head(pilot)

## ----estimate-----------------------------------------------------------------
## Replicate-level p: positive replicates among samples with any detection.
pos_samples <- pilot[pilot$detections > 0, ]
p_hat <- sum(pos_samples$detections) / sum(pos_samples$n_reps)

## Sample-level: positive samples among sites with any detection.
site_tab <- aggregate(detections ~ site, pilot, FUN = function(x) sum(x > 0))
sites_with <- site_tab$site[site_tab$detections > 0]
samp_in_pos_sites <- pilot[pilot$site %in% sites_with, ]
sample_pos_rate <-
  mean(tapply(samp_in_pos_sites$detections, samp_in_pos_sites$site,
              function(x) mean(x > 0)))

## Site-level: observed occupancy, corrected for missed-but-occupied sites.
n_samples_per_site <- max(site_tab$detections[site_tab$detections > 0], 1)
obs_site_rate <- mean(site_tab$detections > 0)
## Probability an occupied site is detected at least once given the pilot design,
## used to inflate the observed rate back toward true occupancy.
detect_given_occ <- 1 - (1 - sample_pos_rate)^8
psi_hat <- min(obs_site_rate / detect_given_occ, 0.999)

theta_hat <- c(site = psi_hat, sample = sample_pos_rate, replicate = p_hat)
rbind(truth = truth$theta, estimated = round(theta_hat, 3))

## ----allocate-----------------------------------------------------------------
design <- allocate(
  theta  = theta_hat,
  cost   = truth$cost,
  goal   = "detection",
  budget = 2000
)
design
summary(design)

## ----frontier-----------------------------------------------------------------
plot(design)

## ----check--------------------------------------------------------------------
brute <- enumerate_designs(
  theta  = theta_hat,
  cost   = truth$cost,
  budget = 2000,
  goal   = "detection"
)
c(one_pass = round(design$P1, 4),
  enumerated = round(brute$P1, 4))

## ----curve-check--------------------------------------------------------------
c(reported   = round(design$P1, 4),
  from_curve = round(detection_from_budget(design$kappa1, design$C), 4))

## ----naive--------------------------------------------------------------------
naive_counts <- function(theta, cost, budget, per_level = 3) {
  n <- c(1, per_level, per_level)
  C1 <- unit_cost_full(cost, n)[1]
  n[1] <- floor(budget / C1)
  detection_recursion(theta, n)$P1
}

c(optimized = round(design$P1, 4),
  naive_3x3 = round(naive_counts(theta_hat, truth$cost, 2000), 4))

## ----identity-----------------------------------------------------------------
optimality_identity(
  allocate(theta_hat, truth$cost, goal = "detection",
           budget = 2000, integer = FALSE)
)

