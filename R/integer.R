#' Exhaustive integer optimum of a nested design
#'
#' Finds the best integer design by evaluating every feasible combination of
#' inner replicate counts, rather than rounding the continuous solution. This is
#' the reference computation of Supplement S2: at small budgets the integer
#' design surface has many local optima (11 to 20 percent of feasible four-level
#' designs in the manuscript), so a hill climb or a rounded continuous solution
#' can miss the best design, and only enumeration is guaranteed correct.
#'
#' Use this when the budget buys only a handful of outermost units, which is the
#' regime where rounding in place (the default in [multiDetectR::allocate()]) is unreliable.
#' The manuscript finds the two agree once the budget buys about ten or more
#' outermost units.
#'
#' @param theta Numeric vector of per-level detectabilities in `(0, 1)`, ordered
#'   outermost to innermost.
#' @param cost Numeric vector of per-level marginal unit costs, same length as
#'   `theta`.
#' @param budget Total budget `B`, `> 0`.
#' @param goal `"detection"` or `"count"`, as in [multiDetectR::allocate()].
#' @param max_n Integer upper bound on each inner count searched. The default
#'   grows with the budget; raise it if the optimum appears to sit at the bound.
#'
#' @return A `multidetect_design` object (see [multiDetectR::allocate()]) for the enumerated
#'   optimum, with an added `enumerated = TRUE` flag.
#'
#' @examples
#' enumerate_designs(theta = c(0.3, 0.5, 0.6), cost = c(200, 20, 2),
#'                   budget = 1500, goal = "detection")
#' @export
enumerate_designs <- function(theta, cost, budget,
                              goal = c("detection", "count"), max_n = NULL) {
  goal <- match.arg(goal)
  check_levels(theta, cost, prob = TRUE)
  if (!is.numeric(budget) || length(budget) != 1L || budget <= 0) {
    stop("`budget` must be a single positive number.", call. = FALSE)
  }
  d <- length(theta)

  objective <- function(n) {
    if (goal == "detection") {
      detection_recursion(theta, n)$P1
    } else {
      ## Count goal: expected positive level-1 samples up to target properties,
      ## proportional to n_1 * P_2 (manuscript Problem 3).
      n[1] * detection_recursion(theta, n)$P[2]
    }
  }

  if (d == 1L) {
    n1 <- floor(budget / cost[1])
    if (n1 < 1) stop("Budget buys fewer than one unit.", call. = FALSE)
    best_n <- n1
  } else {
    ## Reasonable cap on inner counts: no inner level can exceed what the
    ## budget could buy if everything else were a single unit.
    if (is.null(max_n)) {
      max_n <- max(2, ceiling(budget / min(cost)))
      max_n <- min(max_n, 200L)   # keep the grid tractable
    }
    inner_grid <- expand.grid(
      rep(list(seq_len(max_n)), d - 1),
      KEEP.OUT.ATTRS = FALSE
    )
    best_val <- -Inf
    best_n <- NULL
    for (row in seq_len(nrow(inner_grid))) {
      n_inner <- as.integer(inner_grid[row, ])
      n <- c(1L, n_inner)
      C1 <- unit_cost_full(cost, n)[1]
      n1 <- floor(budget / C1)
      if (n1 < 1) next
      n[1] <- n1
      val <- objective(n)
      if (val > best_val) {
        best_val <- val
        best_n <- n
      }
    }
    if (is.null(best_n)) {
      stop("No feasible integer design within the budget; the smallest ",
           "design already exceeds it.", call. = FALSE)
    }
  }

  design_from_counts(theta, cost, best_n, goal)
}

#' Build a design object from fixed integer counts
#'
#' Internal helper that assembles a `multidetect_design` (prices, hazards, cost,
#' detection) from a given count vector, so enumeration and the outward pass
#' return the same shape of object.
#'
#' @inheritParams enumerate_designs
#' @param n Integer count vector.
#' @param goal Design goal, carried onto the object.
#' @return A `multidetect_design` with `enumerated = TRUE`.
#' @keywords internal
#' @noRd
design_from_counts <- function(theta, cost, n, goal) {
  d <- length(theta)
  h <- numeric(d)
  kappa <- numeric(d)
  u <- numeric(d)
  h[d] <- nats(theta[d])
  kappa[d] <- cost[d] / h[d]
  C_full <- unit_cost_full(cost, n)
  if (d >= 2) {
    for (i in (d - 1):1) {
      u[i + 1] <- n[i + 1] * h[i + 1]
      P_in <- -expm1(-u[i + 1])
      h[i] <- nats(theta[i] * P_in)
      kappa[i] <- C_full[i] / h[i]
    }
  }
  structure(
    list(n = n, P1 = detection_recursion(theta, n)$P1,
         C = total_cost(cost, n), kappa1 = kappa[1],
         kappa = kappa, h = h, u = u,
         goal = goal, theta = theta, cost = cost,
         integer = TRUE, enumerated = TRUE),
    class = "multidetect_design"
  )
}
