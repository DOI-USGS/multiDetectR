#' Cost ratio as a function of effort (detection goal)
#'
#' Evaluates the right-hand side of the one-level first-order condition,
#' equation (14) in the manuscript, giving the cost ratio `r` as a function of
#' the dimensionless effort `u` for a given per-level detectability `theta`.
#' Because the function is strictly increasing in `u`, inverting it recovers the
#' optimal effort `u*(r, theta)`.
#'
#' @details This uses the numerically stable form of equation (14) derived in
#'   Supplement S1.4, equation (S1.6). Equation (14) must **not** be evaluated as
#'   written: its denominator `w - 1 + theta` is the difference of order-one
#'   quantities whose true value is `theta * exp(-u)`, so it loses all
#'   significance once `theta * exp(-u)` approaches machine epsilon (roughly
#'   `u > 23 + log(theta)`), which makes the naive form appear non-monotone and
#'   breaks the inversion. Substituting the identity `w - 1 + theta =
#'   theta * exp(-u)` and taking the logarithm with [log1p()] cures it.
#'
#' @param u Numeric dimensionless effort, `>= 0`.
#' @param theta Per-level detectability in `(0, 1)`.
#'
#' @return The cost ratio `r` corresponding to effort `u` at detectability
#'   `theta`.
#' @keywords internal
#' @noRd
r_of_u <- function(u, theta) {
  ## a = theta * (1 - exp(-u)) = 1 - w, kept via expm1 for small u.
  a <- theta * -expm1(-u)
  w <- 1 - a
  ## -w * log(w) / (theta * exp(-u)) - u, with log(w) = log1p(-a).
  -(w * log1p(-a)) / (theta * exp(-u)) - u
}

#' Solve one level of the allocation problem
#'
#' Returns the dimensionless effort `u*` placed inside one unit at a level,
#' given the level's cost ratio `r` and (for the detection goal) its
#' detectability `theta`. This is the scalar step of the outward pass.
#'
#' @details Under the detection goal, and at every level except level 1 under
#'   the count goal, `u*` inverts equation (14) (manuscript), solved here with
#'   a numerically stable form (Supplement S1.4). Under the count goal at
#'   level 1, `u*` solves equation (15), `exp(u) - u = 1 + r`, in which `theta`
#'   does not appear. Both are solved by a monotone root find; equation (15) is
#'   the limit of equation (14) as `theta` approaches zero.
#'
#' @param r Non-negative dimensionless cost ratio (cost divided by the price of
#'   a detection at the level within, equation (13)).
#' @param theta Per-level detectability in `(0, 1)`; ignored for
#'   `goal = "count"`.
#' @param goal Either `"detection"` (equation (14)) or `"count"` (equation
#'   (15)).
#' @param upper Upper bound on `u` for the root search. The default spans the
#'   range a real design can reach.
#'
#' @return The optimal dimensionless effort `u*`, a non-negative scalar.
#'
#' @examples
#' solve_level(r = 10, theta = 0.3, goal = "detection")
#' solve_level(r = 10, goal = "count")
#' @export
solve_level <- function(r, theta = NULL, goal = c("detection", "count"),
                        upper = 50) {
  goal <- match.arg(goal)
  if (!is.numeric(r) || length(r) != 1L || !is.finite(r) || r < 0) {
    stop("`r` must be a single non-negative finite number.", call. = FALSE)
  }
  if (r == 0) {
    return(0)
  }

  if (goal == "count") {
    ## Equation (15): exp(u) - u = 1 + r, theta-free.
    f <- function(u) expm1(u) - u - r   # exp(u) - 1 - u - r
  } else {
    if (is.null(theta) || theta <= 0 || theta >= 1) {
      stop("`theta` in (0, 1) is required for the detection goal.",
           call. = FALSE)
    }
    ## Equation (14) via the stable S1.6 form: find u with r_of_u(u) = r.
    f <- function(u) r_of_u(u, theta) - r
  }

  ## Both targets are strictly increasing in u and f(0) < 0, so grow the
  ## bracket until it contains the root, then solve.
  hi <- upper
  while (f(hi) < 0 && hi < 1e6) {
    hi <- hi * 2
  }
  stats::uniroot(f, lower = 0, upper = hi, tol = .Machine$double.eps^0.5)$root
}

#' Optimal effort via the Lambert W form (count goal)
#'
#' Closed-form evaluation of equation (15) using the lower branch of the
#' Lambert W function, equation (S1.5) in the manuscript:
#' `u* = -(1 + r) - W_{-1}(-exp(-(1 + r)))`. Provided for interpretation and
#' speed; [multiDetectR::solve_level()] with `goal = "count"` solves equation (15) directly
#' and is preferred, because `exp(-(1 + r))` underflows for `r` beyond about
#' 700 while the direct solve continues to work.
#'
#' @param r Non-negative dimensionless cost ratio.
#'
#' @return The optimal dimensionless effort `u*`, or `NA_real_` when the
#'   argument underflows and the Lambert form cannot be evaluated.
#'
#' @examples
#' solve_count_lambert(10)
#' @export
solve_count_lambert <- function(r) {
  if (!is.numeric(r) || length(r) != 1L || !is.finite(r) || r < 0) {
    stop("`r` must be a single non-negative finite number.", call. = FALSE)
  }
  z <- -exp(-(1 + r))
  if (z == 0) {
    ## Underflow: fall back to the direct solve rather than returning Inf.
    return(NA_real_)
  }
  -(1 + r) - lambert_wm1(z)
}

#' Lower branch of the Lambert W function
#'
#' Self-contained evaluation of `W_{-1}(x)` for `x` in `[-1/e, 0)`, the branch
#' with `W_{-1}(x) <= -1`. Implemented here so the package carries no external
#' dependency for a single special function. Uses an asymptotic initial guess
#' refined by Halley's method.
#'
#' @param x Numeric scalar in `[-1/e, 0)`.
#'
#' @return `W_{-1}(x)`.
#' @keywords internal
#' @noRd
lambert_wm1 <- function(x) {
  if (x < -exp(-1) || x >= 0) {
    stop("`x` must be in [-1/e, 0) for the W_{-1} branch.", call. = FALSE)
  }
  ## Asymptotic series near 0^- gives a good starting value (Corless et al.).
  L1 <- log(-x)
  L2 <- log(-L1)
  w <- L1 - L2 + L2 / L1
  ## Halley refinement of w * exp(w) = x.
  for (i in 1:50) {
    ew <- exp(w)
    f <- w * ew - x
    denom <- ew * (w + 1) - (w + 2) * f / (2 * w + 2)
    step <- f / denom
    w <- w - step
    if (abs(step) <= 1e-14 * (abs(w) + 1e-14)) {
      break
    }
  }
  w
}
