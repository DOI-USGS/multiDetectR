#' Allocate a budget across branches (methods or places)
#'
#' Branches are options that are not nested within one another: separate places
#' to sample, or separate methods applied to the same material. Given the
#' detection hazard each branch delivers per unit of budget against each target,
#' this returns the budget split that maximizes the weighted chance of detecting
#' the targets, together with the concentrate-versus-spread test of Theorem 2
#' (Supplement S1.5).
#'
#' @param rate A targets-by-branches numeric matrix `r_tk >= 0`: the detection
#'   hazard (nats) branch `k` delivers per unit of budget against target `t`.
#' @param budget Total budget `B`, `> 0`.
#' @param weight Optional per-target weights `w_t` (default equal).
#' @param prevalence Optional per-target prevalences `psi_t` (default equal).
#'
#' @return A list with `share` (the optimal fraction of the budget to each
#'   branch), `concentrate` (logical, whether the optimum puts the whole budget
#'   on one branch), `detection` (the weighted expected detection at the
#'   optimum), and `regime` (`"dominance"`, `"small-budget"`, or `"large-budget"`
#'   describing which part of Theorem 2 decides the corner test here).
#'
#' @details The objective is
#'   `D(x) = sum_t w_t psi_t (1 - exp(-sum_k r_tk x_k))` (equation (S1.7)). Its
#'   gradient (equation (S1.8)) shares the undetected-weight `lambda_t` across
#'   branches, so a branch takes the whole budget exactly when its gradient is
#'   largest at that corner. Theorem 2: dominance (a branch best for every
#'   target) holds the corner at every budget; at small budgets the test is a
#'   weighted average over targets; at large budgets it is the branch's own
#'   worst target. Separate places never take the whole budget, since a branch
#'   at one place has rate zero against a target elsewhere.
#'
#' @examples
#' ## One target, two methods: the better method takes everything.
#' allocate_branches(rate = matrix(c(0.9, 0.1), nrow = 1), budget = 20)
#' @export
allocate_branches <- function(rate, budget, weight = NULL, prevalence = NULL) {
  if (!is.matrix(rate) || !is.numeric(rate) || any(rate < 0) ||
        any(!is.finite(rate))) {
    stop("`rate` must be a non-negative numeric targets-by-branches matrix.",
         call. = FALSE)
  }
  if (!is.numeric(budget) || length(budget) != 1L || budget <= 0) {
    stop("`budget` must be a single positive number.", call. = FALSE)
  }
  n_t <- nrow(rate)
  n_k <- ncol(rate)
  if (is.null(weight)) weight <- rep(1, n_t)
  if (is.null(prevalence)) prevalence <- rep(1, n_t)
  wpsi <- weight * prevalence

  ## Objective and gradient on the simplex sum(x) = budget, x >= 0.
  negD <- function(x) {
    hz <- as.vector(rate %*% x)
    -sum(wpsi * -expm1(-hz))
  }

  if (n_k == 1L) {
    share <- 1
    best_x <- budget
  } else {
    ## Optimize over the simplex via a softmax reparameterization, from several
    ## starts (including each corner) to avoid a local optimum.
    obj <- function(z) {
      x <- budget * exp(z - max(z)) / sum(exp(z - max(z)))
      negD(x)
    }
    starts <- c(list(rep(0, n_k)),
                lapply(seq_len(n_k), function(k) {
                  z <- rep(-5, n_k); z[k] <- 5; z
                }))
    best <- Inf
    best_x <- NULL
    for (z0 in starts) {
      fit <- stats::optim(z0, obj, method = "Nelder-Mead",
                          control = list(reltol = 1e-12, maxit = 3000))
      if (fit$value < best) {
        best <- fit$value
        z <- fit$par
        best_x <- budget * exp(z - max(z)) / sum(exp(z - max(z)))
      }
    }
    share <- best_x / budget
  }

  concentrate <- max(share) > 1 - 1e-3
  detection <- -negD(best_x)
  list(
    share = share,
    concentrate = concentrate,
    detection = detection,
    regime = branch_regime(rate, budget, wpsi)
  )
}

#' Which Theorem 2 regime decides the corner test
#'
#' Internal helper reporting whether, for the current rates and budget, the
#' concentrate test is governed by dominance, the small-budget weighted average,
#' or the large-budget worst-target rule (Supplement S1.5).
#'
#' @inheritParams allocate_branches
#' @param wpsi Combined weight-times-prevalence vector.
#' @return A single string, one of `"dominance"`, `"small-budget"`,
#'   `"large-budget"`.
#' @keywords internal
#' @noRd
branch_regime <- function(rate, budget, wpsi) {
  n_k <- ncol(rate)
  ## Dominance: some branch is >= every rival for every target.
  dominates <- vapply(seq_len(n_k), function(k) {
    all(rate[, k] >= rate[, -k, drop = FALSE])
  }, logical(1))
  if (any(dominates)) {
    return("dominance")
  }
  ## Otherwise classify by how far the budget has pushed lambda toward the
  ## worst-target regime: compare the per-target undetected weights at the
  ## best single-branch spend.
  best_k <- which.max(colSums(wpsi * rate))
  hz <- rate[, best_k] * budget
  lambda <- wpsi * exp(-hz)
  if (max(lambda) / sum(lambda) > 0.9) "large-budget" else "small-budget"
}
