#' Detection probability of a nested design
#'
#' Computes the probability of at least one detection somewhere in a nested
#' survey, given that the target is present, from the per-level detectabilities
#' and replicate counts. This is the recursion of equations (1) and (2) in the
#' manuscript, run from the innermost level outward.
#'
#' Levels are ordered from the outermost (reporting) level, index 1, to the
#' innermost level, index `d`. `theta[i]` is the probability that a level-`i`
#' unit carries the target given that the unit containing it does, and `n[i]`
#' is the number of level-`i` units taken within one level-`(i - 1)` unit.
#'
#' @param theta Numeric vector of per-level conditional detectabilities, each in
#'   `(0, 1)`, ordered outermost to innermost.
#' @param n Integer or numeric vector of per-level replicate counts, the same
#'   length as `theta`, each `>= 1`.
#'
#' @return A list with elements `P` (the per-level vector `P_1, ..., P_d`, where
#'   `P[i]` is the probability of at least one detection among the `n[i]` units
#'   within one carrying level-`(i - 1)` unit) and `P1` (the overall detection
#'   probability `P[1]`).
#'
#' @details Equation (1): `P_d = 1 - (1 - theta_d)^{n_d}`. Equation (2):
#'   `P_i = 1 - (1 - theta_i P_{i + 1})^{n_i}` for `i = d - 1, ..., 1`.
#'
#' @examples
#' detection_recursion(theta = c(0.3, 0.5, 0.6), n = c(5, 3, 2))
#' @export
detection_recursion <- function(theta, n) {
  check_levels(theta, n)
  d <- length(theta)
  P <- numeric(d)
  ## Equation (1): innermost level.
  P[d] <- 1 - (1 - theta[d])^n[d]
  ## Equation (2): work outward.
  if (d >= 2) {
    for (i in (d - 1):1) {
      P[i] <- 1 - (1 - theta[i] * P[i + 1])^n[i]
    }
  }
  list(P = P, P1 = P[1])
}

#' Total cost of a nested design
#'
#' Equation (3) in the manuscript: the total number of units purchased times
#' their per-level costs, where reaching any level-`i` unit requires paying for
#' every unit that contains it, so costs compound outward.
#'
#' @param cost Numeric vector of per-level marginal unit costs (excluding the
#'   cost of contents), ordered outermost to innermost.
#' @param n Integer or numeric vector of per-level replicate counts, the same
#'   length as `cost`.
#'
#' @return The scalar total cost `C = sum_i cost_i * prod_{j <= i} n_j`.
#'
#' @examples
#' total_cost(cost = c(200, 20, 2), n = c(10, 3, 2))
#' @export
total_cost <- function(cost, n) {
  check_levels(cost, n, prob = FALSE)
  sum(cost * cumprod(n))
}

#' Full cost of a unit including its contents
#'
#' Equation (8) in the manuscript: the cost of one level-`i` unit together with
#' everything nested inside it, built from the innermost level outward.
#'
#' @inheritParams total_cost
#'
#' @return A numeric vector `C_1, ..., C_d`, where `C[i]` is the full cost of one
#'   level-`i` unit including its contents. `C[1]` is the cost of one outermost
#'   unit and the whole survey below it.
#'
#' @examples
#' unit_cost_full(cost = c(200, 20, 2), n = c(10, 3, 2))
#' @export
unit_cost_full <- function(cost, n) {
  check_levels(cost, n, prob = FALSE)
  d <- length(cost)
  C_full <- numeric(d)
  C_full[d] <- cost[d]
  if (d >= 2) {
    for (i in (d - 1):1) {
      C_full[i] <- cost[i] + n[i + 1] * C_full[i + 1]
    }
  }
  C_full
}

#' Detection hazard in nats
#'
#' Converts a detection probability to the detection hazard `-log(1 - P)`,
#' measured in nats, and back. On the nat scale independent replicates add
#' rather than multiply, which is the property the allocation rule relies on
#' (manuscript, "Detection on a log scale").
#'
#' @param P Numeric detection probability in `[0, 1)`.
#' @param nats Numeric detection hazard in nats (`>= 0`).
#'
#' @return `nats()` returns the hazard `-log(1 - P)`; `prob()` returns the
#'   probability `1 - exp(-nats)`.
#'
#' @examples
#' nats(0.95)   # about 3.00
#' prob(3)      # about 0.95
#' @export
nats <- function(P) {
  if (any(P < 0 | P >= 1)) {
    stop("`P` must be in [0, 1).", call. = FALSE)
  }
  -log1p(-P)
}

#' @rdname nats
#' @export
prob <- function(nats) {
  if (any(nats < 0)) {
    stop("`nats` must be >= 0.", call. = FALSE)
  }
  -expm1(-nats)
}

#' Validate per-level inputs
#'
#' Internal helper. Checks that two per-level vectors are numeric, the same
#' length, non-empty, and (for probabilities) in the open interval `(0, 1)`;
#' counts and costs must be positive.
#'
#' @param x,y The two per-level vectors to compare (e.g. `theta` and `n`).
#' @param prob Logical; if `TRUE`, `x` is treated as probabilities and checked
#'   to lie in `(0, 1)`.
#'
#' @return Invisibly `TRUE`; called for its side effect of stopping on invalid
#'   input.
#' @keywords internal
#' @noRd
check_levels <- function(x, y, prob = TRUE) {
  nm_x <- deparse(substitute(x))
  nm_y <- deparse(substitute(y))
  if (!is.numeric(x) || !is.numeric(y)) {
    stop(sprintf("`%s` and `%s` must be numeric.", nm_x, nm_y), call. = FALSE)
  }
  if (length(x) == 0L || length(x) != length(y)) {
    stop(sprintf("`%s` and `%s` must be non-empty and the same length.",
                 nm_x, nm_y), call. = FALSE)
  }
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("Inputs must be finite.", call. = FALSE)
  }
  if (prob && any(x <= 0 | x >= 1)) {
    stop(sprintf("`%s` must be in the open interval (0, 1).", nm_x),
         call. = FALSE)
  }
  if (!prob && any(x <= 0)) {
    stop(sprintf("`%s` must be positive.", nm_x), call. = FALSE)
  }
  if (any(y <= 0)) {
    stop(sprintf("`%s` (counts) must be positive.", nm_y), call. = FALSE)
  }
  invisible(TRUE)
}
