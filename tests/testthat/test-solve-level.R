test_that("r_of_u is strictly increasing (Supplement S1.4 monotonicity)", {
  ## The naive form of equation (14) produces spurious decreases once
  ## theta * exp(-u) underflows; the stable S1.6 form must stay monotone.
  r_of_u <- multiDetectR:::r_of_u
  for (theta in c(0.001, 0.01, 0.1, 0.5, 0.9, 0.99)) {
    u <- seq(0.01, 40, length.out = 4000)
    r <- vapply(u, r_of_u, numeric(1), theta = theta)
    expect_true(all(diff(r) > 0),
                info = sprintf("non-monotone r(u) at theta = %g", theta))
  }
})

test_that("solve_level (detection) inverts a bounded numerical minimization", {
  ## u* from equation (14) should minimize kappa_i = (c + u * kappa_next) / H(u).
  for (theta in c(0.1, 0.4, 0.8)) {
    for (r in c(0.1, 1, 10, 100)) {
      u_star <- solve_level(r, theta = theta, goal = "detection")
      H <- function(u) -log1p(-theta * -expm1(-u))
      f <- function(u) (r + u) / H(u)             # objective in units of kappa_next
      u_num <- stats::optimize(f, c(1e-6, 60))$minimum
      expect_equal(u_star, u_num, tolerance = 1e-3,
                   info = sprintf("theta = %g, r = %g", theta, r))
    }
  }
})

test_that("solve_level (count) solves exp(u) - u = 1 + r", {
  for (r in c(0.1, 1, 10, 100)) {
    u_star <- solve_level(r, goal = "count")
    expect_equal(exp(u_star) - u_star, 1 + r, tolerance = 1e-6)
  }
})

test_that("detection goal approaches the count goal as theta -> 0", {
  ## Equation (15) is the theta -> 0 limit of equation (14).
  r <- 5
  u_count <- solve_level(r, goal = "count")
  u_small <- solve_level(r, theta = 1e-4, goal = "detection")
  expect_equal(u_small, u_count, tolerance = 1e-2)
})

test_that("Lambert W form agrees with the direct count solve", {
  for (r in c(0.01, 1, 10, 100, 500)) {
    expect_equal(solve_count_lambert(r), solve_level(r, goal = "count"),
                 tolerance = 1e-8, info = sprintf("r = %g", r))
  }
})

test_that("solve_level handles the r = 0 boundary", {
  expect_equal(solve_level(0, theta = 0.5, goal = "detection"), 0)
  expect_equal(solve_level(0, goal = "count"), 0)
})
