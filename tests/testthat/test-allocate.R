## Full numerical optimum of the continuous problem, sharing none of the
## package's decomposition: directly minimize kappa_1 over the inner efforts
## u_2, ..., u_d, then compare to the one-pass design.
kappa1_from_u <- function(u_inner, theta, cost) {
  d <- length(theta)
  h <- numeric(d)
  C_full <- numeric(d)
  h[d] <- -log1p(-theta[d])
  C_full[d] <- cost[d]
  if (d >= 2) {
    for (i in (d - 1):1) {
      u <- u_inner[i]                       # effort inside a level-i unit
      n_in <- u / h[i + 1]
      P_in <- -expm1(-u)
      h[i] <- -log1p(-theta[i] * P_in)
      C_full[i] <- cost[i] + n_in * C_full[i + 1]
    }
  }
  C_full[1] / h[1]
}

numeric_optimum <- function(theta, cost) {
  d <- length(theta)
  obj <- function(log_u) kappa1_from_u(exp(log_u), theta, cost)
  best <- Inf
  best_u <- NULL
  ## Multiple starts to avoid any local trap in the direct optimizer.
  for (s in 1:8) {
    set.seed(s)
    start <- log(stats::runif(d - 1, 0.1, 5))
    fit <- stats::optim(start, obj, method = "Nelder-Mead",
                        control = list(reltol = 1e-12, maxit = 5000))
    if (fit$value < best) {
      best <- fit$value
      best_u <- exp(fit$par)
    }
  }
  list(kappa1 = best, u_inner = best_u)
}

test_that("one-pass kappa_1 matches full numerical optimization (250 designs)", {
  set.seed(7)
  n_designs <- 250
  max_rel <- 0
  for (k in seq_len(n_designs)) {
    d <- sample(2:5, 1)
    theta <- stats::runif(d, 0.01, 0.97)
    cost <- 10^stats::runif(d, 0, 4)                 # log-uniform, 4 orders
    design <- allocate(theta, cost, goal = "detection",
                       budget = 1e6, integer = FALSE)
    num <- numeric_optimum(theta, cost)
    rel <- abs(design$kappa1 - num$kappa1) / num$kappa1
    max_rel <- max(max_rel, rel)
    expect_lt(rel, 1e-4)
  }
  expect_lt(max_rel, 1e-4)
})

test_that("budget invariance: inner counts do not depend on the budget", {
  theta <- c(0.3, 0.5, 0.6)
  cost <- c(200, 20, 2)
  d1 <- allocate(theta, cost, goal = "detection", budget = 2000)
  d2 <- allocate(theta, cost, goal = "detection", budget = 20000)
  ## Inner replication is budget-free; only n_1 should change.
  expect_equal(d1$n[-1], d2$n[-1])
  expect_gt(d2$n[1], d1$n[1])
})

test_that("the two design problems give the same inner design", {
  theta <- c(0.3, 0.5, 0.6)
  cost <- c(200, 20, 2)
  by_budget <- allocate(theta, cost, goal = "detection", budget = 5000,
                        integer = FALSE)
  by_target <- allocate(theta, cost, goal = "detection", P_target = 0.9,
                        integer = FALSE)
  ## Problems 1 and 2 both minimize kappa_1, so inner efforts coincide.
  expect_equal(by_budget$u[-1], by_target$u[-1], tolerance = 1e-6)
  expect_equal(by_budget$kappa1, by_target$kappa1, tolerance = 1e-8)
})

test_that("cost_for_target and detection_from_budget invert each other", {
  theta <- c(0.3, 0.5, 0.6)
  cost <- c(200, 20, 2)
  des <- allocate(theta, cost, goal = "detection", budget = 5000,
                  integer = FALSE)
  cost_90 <- cost_for_target(des$kappa1, 0.90)
  expect_equal(detection_from_budget(des$kappa1, cost_90), 0.90,
               tolerance = 1e-8)
})

test_that("realized P1 agrees with the budget-detection curve, equation (10)", {
  theta <- c(0.3, 0.5, 0.6)
  cost <- c(200, 20, 2)
  des <- allocate(theta, cost, goal = "detection", budget = 5000,
                  integer = FALSE)
  ## With continuous n_1 the budget is fully spent, so P1 = 1 - exp(-B/kappa1).
  expect_equal(des$P1, detection_from_budget(des$kappa1, des$C),
               tolerance = 1e-6)
})

test_that("input guards fire", {
  expect_error(allocate(c(0.3, 0.5), c(200, 2), goal = "detection"),
               "exactly one")
  expect_error(allocate(c(0.3, 0.5), c(200, 2), goal = "detection",
                        budget = 1000, P_target = 0.9),
               "exactly one")
  expect_error(allocate(c(0.3, 0.5), c(200, 2), goal = "count",
                        P_target = 0.9),
               "count goal")
})
