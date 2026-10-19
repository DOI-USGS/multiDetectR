#' Price of a detection from budget, and cost to reach a target
#'
#' Once a design is reduced to the price of a detection at the outermost level,
#' `kappa_1` (equation (9)), the whole budget-detection relationship is a single
#' exponential. `detection_from_budget()` is equation (10),
#' `P(B) = 1 - exp(-B / kappa_1)`; `cost_for_target()` is equation (11),
#' `kappa_1 * (-log(1 - P*))`.
#'
#' @param kappa1 Price of a detection at the outermost level (nats per unit
#'   cost), `> 0`.
#' @param budget Total budget `B`, `> 0`.
#' @param P_target Required detection probability `P*` in `(0, 1)`.
#'
#' @return `detection_from_budget()` returns a detection probability;
#'   `cost_for_target()` returns a cost.
#'
#' @examples
#' detection_from_budget(kappa1 = 50, budget = 150)
#' cost_for_target(kappa1 = 50, P_target = 0.95)
#' @export
detection_from_budget <- function(kappa1, budget) {
  if (kappa1 <= 0 || budget < 0) {
    stop("`kappa1` must be > 0 and `budget` >= 0.", call. = FALSE)
  }
  -expm1(-budget / kappa1)
}

#' @rdname detection_from_budget
#' @export
cost_for_target <- function(kappa1, P_target) {
  if (kappa1 <= 0 || P_target <= 0 || P_target >= 1) {
    stop("`kappa1` must be > 0 and `P_target` in (0, 1).", call. = FALSE)
  }
  kappa1 * nats(P_target)
}

#' Allocate detection effort across nested survey levels
#'
#' The main entry point. Given per-level detectabilities and costs, returns the
#' design that best spends a budget (or most cheaply reaches a required
#' detection probability) by running the exact innermost-outward pass of the
#' manuscript: price the innermost level with equation (9), then for each level
#' outward form the cost ratio (equation (13)), solve one scalar problem
#' (equation (14) or (15)), and update the price of a detection (equation (12)).
#'
#' @param theta Numeric vector of per-level detectabilities in `(0, 1)`, ordered
#'   outermost (level 1) to innermost (level `d`).
#' @param cost Numeric vector of per-level marginal unit costs, same length as
#'   `theta`.
#' @param goal `"detection"` (maximize the chance of at least one detection,
#'   Problems 1 and 2) or `"count"` (maximize the expected number of positive
#'   level-1 samples, Problem 3).
#' @param budget Total budget `B` (Problem 1 or 3). Supply exactly one of
#'   `budget` or `P_target`.
#' @param P_target Required detection probability `P*` (Problem 2).
#' @param integer Logical; if `TRUE` (default) inner replicate counts are
#'   integers, obtained by rounding in place (Supplement S2). If `FALSE`,
#'   effort is treated as continuous and needs no rounding.
#'
#' @return An object of class `multidetect_design`: a list with the per-level
#'   counts `n`, the realized detection probability `P1`, the total cost `C`,
#'   the outermost price of a detection `kappa1`, the per-level vectors `kappa`,
#'   `h`, and `u`, plus the `goal`, `theta`, `cost`, and whether the design is
#'   `integer`.
#'
#' @details The budget enters only at the last step, setting the number of
#'   outermost units `n_1`; replication at every inner level is budget-free
#'   (manuscript, "One number summarises a hierarchy"). Under the count goal,
#'   `n_1` is left as the budget-feasible count and level 1 uses equation (15);
#'   under the detection goal every level uses equation (14).
#'
#' @examples
#' allocate(theta = c(0.3, 0.5, 0.6), cost = c(200, 20, 2),
#'          goal = "detection", budget = 5000)
#' @export
allocate <- function(theta, cost, goal = c("detection", "count"),
                     budget = NULL, P_target = NULL, integer = TRUE) {
  goal <- match.arg(goal)
  check_levels(theta, cost, prob = TRUE)
  if (is.null(budget) == is.null(P_target)) {
    stop("Supply exactly one of `budget` or `P_target`.", call. = FALSE)
  }
  if (!is.null(P_target) && goal == "count") {
    stop("`P_target` applies to the detection goal; use `budget` for the ",
         "count goal.", call. = FALSE)
  }
  d <- length(theta)

  ## Innermost level: no choice. kappa_d = c_d / h_d, h_d = -log(1 - theta_d).
  h <- numeric(d)
  kappa <- numeric(d)
  u <- numeric(d)
  n <- rep(1, d)
  h[d] <- nats(theta[d])
  kappa[d] <- cost[d] / h[d]

  ## Outward pass over the inner levels d-1, ..., 1.
  if (d >= 2) {
    for (i in (d - 1):1) {
      r_i <- cost[i] / kappa[i + 1]                       # equation (13)
      level_goal <- if (goal == "count" && i == 1) "count" else "detection"
      u_star <- solve_level(r_i, theta = theta[i], goal = level_goal)
      u[i + 1] <- u_star
      ## Convert effort to a count of the units nested one level in.
      n[i + 1] <- u_star / h[i + 1]
      if (integer) {
        n[i + 1] <- max(1, round(n[i + 1]))               # round in place (S2)
        u[i + 1] <- n[i + 1] * h[i + 1]
      }
      ## Recompute this level's hazard and price from the (rounded) count.
      P_in <- -expm1(-u[i + 1])                           # 1 - exp(-u_{i+1})
      h[i] <- nats(theta[i] * P_in)                       # equation (12) hazard
      C_i <- cost[i] + n[i + 1] * (kappa[i + 1] * h[i + 1])
      kappa[i] <- C_i / h[i]                              # equation (12)
    }
  }

  ## Outermost count n_1 from the budget (equation (10)/(11) territory).
  C_full <- unit_cost_full(cost, n)
  C1 <- C_full[1]
  if (!is.null(budget)) {
    n[1] <- if (integer) floor(budget / C1) else budget / C1
    if (integer && n[1] < 1) {
      warning("Budget buys fewer than one outermost unit; returning n_1 = 0. ",
              "Consider enumerate_designs() at small budgets (Supplement S2).",
              call. = FALSE)
      n[1] <- 0
    }
  } else {
    ## Problem 2: n_1 implied by the required detection probability.
    need_nats <- nats(P_target)
    n[1] <- if (integer) ceiling(need_nats / h[1]) else need_nats / h[1]
  }

  C <- total_cost(cost, n)
  P1 <- if (n[1] >= 1) detection_recursion(theta, n)$P1 else 0
  kappa1 <- kappa[1]

  structure(
    list(n = n, P1 = P1, C = C, kappa1 = kappa1,
         kappa = kappa, h = h, u = u,
         goal = goal, theta = theta, cost = cost, integer = integer),
    class = "multidetect_design"
  )
}

#' @export
print.multidetect_design <- function(x, ...) {
  d <- length(x$theta)
  cat(sprintf("<multidetect_design> depth %d, goal = %s, %s counts\n",
              d, x$goal, if (x$integer) "integer" else "continuous"))
  cat(sprintf("  per-level n:     %s\n",
              paste(format(x$n, digits = 4), collapse = ", ")))
  cat(sprintf("  detection P1:    %.4f\n", x$P1))
  cat(sprintf("  total cost C:    %.4g\n", x$C))
  cat(sprintf("  price kappa_1:   %.4g (cost per nat at the outermost level)\n",
              x$kappa1))
  invisible(x)
}

#' @export
summary.multidetect_design <- function(object, ...) {
  d <- length(object$theta)
  data.frame(
    level = seq_len(d),
    theta = object$theta,
    cost = object$cost,
    n = object$n,
    hazard = object$h,
    kappa = object$kappa
  )
}

#' Plot the budget-detection frontier of a design
#'
#' Draws the detection probability a design achieves as a function of the
#' budget, `P(B) = 1 - exp(-B / kappa_1)` (equation (10) in the manuscript),
#' using base graphics. The design's own budget and realized detection are
#' marked. The curve is smooth and strictly concave with no natural stopping
#' point, so any "knee" is a property of the axes rather than of the design
#' (manuscript, "A continuous frontier has no natural stopping point").
#'
#' @param x A `multidetect_design` object.
#' @param max_budget Right-hand limit of the budget axis. Defaults to the budget
#'   that reaches 0.999 detection, so the whole informative range is shown.
#' @param n_points Number of points along the curve.
#' @param ... Passed to [graphics::plot()].
#'
#' @return The design, invisibly; called for the plot it draws.
#'
#' @examples
#' des <- allocate(theta = c(0.3, 0.5, 0.6), cost = c(200, 20, 2),
#'                 goal = "detection", budget = 5000)
#' plot(des)
#' @export
plot.multidetect_design <- function(x, max_budget = NULL, n_points = 200, ...) {
  kappa1 <- x$kappa1
  if (is.null(max_budget)) {
    max_budget <- kappa1 * nats(0.999)
  }
  B <- seq(0, max_budget, length.out = n_points)
  P <- -expm1(-B / kappa1)
  graphics::plot(B, P, type = "l",
                 xlab = "Budget", ylab = "Detection probability",
                 ylim = c(0, 1), ...)
  if (x$C <= max_budget && x$P1 > 0) {
    graphics::points(x$C, x$P1, pch = 19)
    graphics::segments(x$C, 0, x$C, x$P1, lty = 3)
    graphics::segments(0, x$P1, x$C, x$P1, lty = 3)
  }
  invisible(x)
}
