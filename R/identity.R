#' Audit a design against the optimality identity A = B
#'
#' At every level that contains another, an optimal design satisfies one
#' identity (Proposition S1.8 in the manuscript): the elasticity of the level's
#' detection hazard with respect to the effort inside it, `A_i`, equals the
#' share of the level's cost that sits inside it, `B_i`. In words, an optimal
#' design spends on a unit's contents the same share of that unit's budget as
#' those contents contribute to its detection power.
#'
#' This lets a planner check a costed design, level by level, without re-solving
#' the allocation: where `A_i` and `B_i` diverge, that level is over- or
#' under-replicated relative to the optimum.
#'
#' @param design A `multidetect_design` object from [multiDetectR::allocate()] or
#'   [multiDetectR::enumerate_designs()].
#'
#' @return A data frame with one row per level that has contents (levels 1 to
#'   `d - 1`): the elasticity `A`, the inside-cost share `B`, and their
#'   difference `A - B` (zero at the continuous optimum).
#'
#' @details `A_i = H_u * u_{i+1} / h_i`, where `H_u` is the derivative of the
#'   hazard with respect to inner effort, and `B_i = n_{i+1} C_{i+1} / C_i`.
#'   For a continuous optimal design `A - B` is zero to numerical tolerance; for
#'   an integer design it is small but nonzero because the counts are rounded.
#'
#' @examples
#' des <- allocate(theta = c(0.3, 0.5, 0.6), cost = c(200, 20, 2),
#'                 goal = "detection", budget = 5000, integer = FALSE)
#' optimality_identity(des)
#' @export
optimality_identity <- function(design) {
  if (!inherits(design, "multidetect_design")) {
    stop("`design` must be a multidetect_design object.", call. = FALSE)
  }
  theta <- design$theta
  cost <- design$cost
  n <- design$n
  d <- length(theta)
  if (d < 2) {
    return(data.frame(level = integer(0), A = numeric(0),
                      B = numeric(0), diff = numeric(0)))
  }

  h <- design$h
  C_full <- unit_cost_full(cost, n)
  A <- numeric(d - 1)
  B <- numeric(d - 1)
  for (i in seq_len(d - 1)) {
    u <- n[i + 1] * h[i + 1]                 # effort inside a level-i unit
    w <- 1 - theta[i] + theta[i] * exp(-u)   # as in equation (14)
    H_u <- theta[i] * exp(-u) / w            # derivative of the hazard
    A[i] <- H_u * u / h[i]
    B[i] <- n[i + 1] * C_full[i + 1] / C_full[i]
  }
  data.frame(level = seq_len(d - 1), A = A, B = B, diff = A - B)
}
