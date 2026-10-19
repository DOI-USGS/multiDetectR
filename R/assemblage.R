#' Allocate detection effort for an assemblage of targets
#'
#' Extends [allocate()] from a single target to several targets (for example,
#' species) that share one nested survey but detect differently. All targets are
#' sampled by the same design, so a single set of per-level counts must serve
#' them all. The rule returns the one design that maximizes the weighted chance
#' of detecting the targets, which the manuscript shows is governed by the
#' hardest-to-detect targets rather than by an average target.
#'
#' @param theta A targets-by-levels numeric matrix of conditional
#'   detectabilities in `(0, 1)`: row `t` is target `t`'s detectability at each
#'   level, ordered outermost (column 1) to innermost (column `d`). A single row
#'   reproduces [allocate()].
#' @param cost Numeric vector of per-level marginal unit costs, length `d`
#'   (shared across targets).
#' @param goal Currently `"detection"`: maximize the weighted probability of at
#'   least one detection per target (equation (S1.7)).
#' @param budget Total budget `B`. Supply exactly one of `budget` or `D_target`.
#' @param D_target Required weighted detection in `(0, 1)`: the cheapest design
#'   whose weighted detection reaches this value.
#' @param weight Optional per-target weights `w_t` (default equal). Rescaled to
#'   sum to one.
#' @param integer Logical; if `TRUE` (default) inner counts are integers.
#'
#' @return An object of class `multidetect_assemblage`: the shared per-level
#'   counts `n`, the weighted detection `D`, the total cost `C`, a per-target
#'   data frame `targets` (weight, detection `P1`, hazard, and the miss-weight
#'   `lambda = w * (1 - P1)` that measures how much each target pulls the
#'   design), the index `binding` of the target with the largest miss-weight,
#'   and the inputs.
#'
#' @details Unlike the single-target problem, the targets couple through the
#'   shared counts, so the inner efforts are chosen by optimizing the weighted
#'   detection `D(n) = sum_t w_t (1 - exp(-H_t(n)))` over the shared design, with
#'   `H_t` target `t`'s detection hazard. The budget then sets the number of
#'   outermost units, exactly as in the single-target rule.
#'
#' @examples
#' theta <- rbind(easy = c(0.6, 0.7, 0.8), hard = c(0.2, 0.3, 0.4))
#' allocate_assemblage(theta, cost = c(200, 20, 2), budget = 5000)
#' @export
allocate_assemblage <- function(theta, cost, goal = "detection",
                                budget = NULL, D_target = NULL,
                                weight = NULL, integer = TRUE) {
  goal <- match.arg(goal)
  theta <- as.matrix(theta)
  if (!is.numeric(theta) || any(theta <= 0 | theta >= 1) ||
        any(!is.finite(theta))) {
    stop("`theta` must be a numeric matrix with entries in (0, 1).",
         call. = FALSE)
  }
  d <- ncol(theta)
  n_t <- nrow(theta)
  if (!is.numeric(cost) || length(cost) != d || any(cost <= 0)) {
    stop("`cost` must be a positive numeric vector of length ncol(theta).",
         call. = FALSE)
  }
  if (is.null(budget) == is.null(D_target)) {
    stop("Supply exactly one of `budget` or `D_target`.", call. = FALSE)
  }
  if (is.null(weight)) weight <- rep(1, n_t)
  if (length(weight) != n_t || any(weight < 0) || all(weight == 0)) {
    stop("`weight` must be a non-negative vector, one per target, not all 0.",
         call. = FALSE)
  }
  w <- weight / sum(weight)

  ## Weighted detection of a given integer design n across all targets.
  weighted_D <- function(n) {
    p1 <- vapply(seq_len(n_t),
                 function(t) detection_recursion(theta[t, ], n)$P1, numeric(1))
    sum(w * p1)
  }

  if (d == 1L) {
    ## No inner levels to choose; only n_1 from the budget / target.
    best_inner <- integer(0)
  } else {
    ## The targets couple through the shared inner counts, so the inner design
    ## cannot be solved target-by-target. Choose it to maximize weighted
    ## detection AT the spending level, with the outermost count set by the
    ## budget: with n_1 = B / C_1 outermost units, each target t accrues hazard
    ## n_1 * h_1t, so D(B) = sum_t w_t (1 - exp(-(B / C_1) * h_1t)). Optimizing
    ## this over the inner counts balances detection against cost (more inner
    ## replication raises each h_1t but also C_1), which is the assemblage
    ## analogue of minimizing kappa_1 in the single-target rule. The inner
    ## design that results does not depend on the absolute budget; we use a
    ## reference budget here and set the real n_1 below.
    ref_budget <- if (!is.null(budget)) {
      budget
    } else {
      ## No budget given (D_target path): use a reference budget that a proxy
      ## weighted-average target would need to reach D_target, so the inner
      ## design is chosen near the spending level the target implies.
      avg_theta <- as.vector(w %*% theta)
      avg_theta <- pmin(pmax(avg_theta, 1e-3), 1 - 1e-3)
      proxy <- allocate(avg_theta, cost, goal = "detection",
                        P_target = min(D_target, 0.999), integer = FALSE)
      proxy$C
    }
    neg_DB <- function(z) {
      n <- c(1, exp(z))
      C1 <- unit_cost_full(cost, n)[1]
      n1 <- ref_budget / C1
      ## per-target outermost hazard h_1t under this inner design
      h1 <- vapply(seq_len(n_t), function(t) {
        th <- theta[t, ]
        P2 <- 1
        if (d >= 2) {
          Pd <- 1 - (1 - th[d])^exp(z)[d - 1]
          P <- Pd
          if (d >= 3) for (i in (d - 1):2) P <- 1 - (1 - th[i] * P)^exp(z)[i - 1]
          P2 <- P
        }
        -log1p(-th[1] * P2)
      }, numeric(1))
      -sum(w * -expm1(-n1 * h1))
    }
    if (d == 2L) {
      ## One inner count: a 1-D search; use optimize() (Brent) over log-effort.
      opt <- stats::optimize(function(z) neg_DB(z), interval = c(log(1), log(1e4)))
      best_u <- exp(opt$minimum)
    } else {
      best <- Inf
      best_u <- rep(2, d - 1)
      for (s in 1:8) {
        set.seed(s)
        st <- log(stats::runif(d - 1, 1, 5))
        fit <- stats::optim(st, neg_DB, method = "Nelder-Mead",
                            control = list(reltol = 1e-10, maxit = 4000))
        if (fit$value < best) {
          best <- fit$value
          best_u <- exp(fit$par)
        }
      }
    }
    best_inner <- best_u
    if (integer) best_inner <- pmax(1, round(best_inner))
  }

  ## Assemble counts with n_1 chosen from the budget or the detection target.
  n <- c(1, best_inner)
  C1 <- unit_cost_full(cost, n)[1]
  if (!is.null(budget)) {
    n[1] <- if (integer) floor(budget / C1) else budget / C1
    if (n[1] < 1) {
      warning("Budget buys fewer than one outermost unit; returning n_1 = 0.",
              call. = FALSE)
      n[1] <- 0
    }
  } else {
    ## Grow n_1 until weighted detection reaches D_target.
    n1 <- 1
    while (weighted_D(c(n1, best_inner)) < D_target && n1 < 1e6) n1 <- n1 + 1
    n[1] <- n1
  }

  ## Per-target summary at the chosen design.
  p1 <- vapply(seq_len(n_t),
               function(t) if (n[1] >= 1) detection_recursion(theta[t, ], n)$P1
                           else 0, numeric(1))
  haz <- ifelse(p1 < 1, -log1p(-p1), Inf)
  lambda <- w * (1 - p1)                 # miss-weight: how much each pulls design
  targets <- data.frame(
    target = if (!is.null(rownames(theta))) rownames(theta) else seq_len(n_t),
    weight = w,
    P1 = p1,
    hazard = haz,
    lambda = lambda,
    row.names = NULL
  )

  structure(
    list(n = n, D = sum(w * p1), C = total_cost(cost, n),
         targets = targets, binding = which.max(lambda),
         goal = goal, theta = theta, cost = cost,
         weight = w, integer = integer),
    class = "multidetect_assemblage"
  )
}

#' Weighted detection for continuous shared counts
#'
#' Internal helper: weighted detection `sum_t w_t P1_t(n)` where the inner
#' counts may be continuous (used while optimizing the shared inner design).
#'
#' @param theta Targets-by-levels detectability matrix.
#' @param n Shared count vector (may be continuous).
#' @param w Normalized target weights.
#' @return Scalar weighted detection.
#' @keywords internal
#' @noRd
weighted_D_cont <- function(theta, n, w) {
  n_t <- nrow(theta)
  d <- ncol(theta)
  p1 <- numeric(n_t)
  for (t in seq_len(n_t)) {
    th <- theta[t, ]
    P <- 1 - (1 - th[d])^n[d]
    if (d >= 2) {
      for (i in (d - 1):1) P <- 1 - (1 - th[i] * P)^n[i]
    }
    p1[t] <- P
  }
  sum(w * p1)
}

#' @export
print.multidetect_assemblage <- function(x, ...) {
  d <- ncol(x$theta)
  cat(sprintf("<multidetect_assemblage> %d targets, depth %d, %s counts\n",
              nrow(x$theta), d, if (x$integer) "integer" else "continuous"))
  cat(sprintf("  shared per-level n: %s\n",
              paste(format(x$n, digits = 4), collapse = ", ")))
  cat(sprintf("  weighted detection D: %.4f\n", x$D))
  cat(sprintf("  total cost C:         %.4g\n", x$C))
  bt <- x$targets$target[x$binding]
  cat(sprintf("  binding target:       %s (largest miss-weight)\n", bt))
  invisible(x)
}

#' @export
summary.multidetect_assemblage <- function(object, ...) {
  object$targets
}
