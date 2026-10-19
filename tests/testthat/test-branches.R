test_that("a dominant branch takes the whole budget (Theorem 2a)", {
  ## One target, method A nine times better: concentrate on A.
  res <- allocate_branches(rate = matrix(c(0.9, 0.1), nrow = 1), budget = 20)
  expect_true(res$concentrate)
  expect_equal(res$share[1], 1, tolerance = 1e-2)
  expect_equal(res$regime, "dominance")
})

test_that("two places (zero cross-rates) are never fully concentrated", {
  ## Target 1 only reachable by branch 1, target 2 only by branch 2.
  rate <- matrix(c(0.5, 0,
                   0,   0.5), nrow = 2, byrow = TRUE)
  res <- allocate_branches(rate, budget = 20)
  expect_false(res$concentrate)
  ## By symmetry the split is even.
  expect_equal(res$share[1], 0.5, tolerance = 1e-2)
})

test_that("a method strong overall but blind to one target eventually shares", {
  ## Branch 1 is far better on target 1 but zero on target 2; branch 2 covers
  ## target 2. At a large budget the money must move toward branch 2.
  rate <- matrix(c(0.9, 0.2,
                   0.0, 0.2), nrow = 2, byrow = TRUE)
  small <- allocate_branches(rate, budget = 1)
  large <- allocate_branches(rate, budget = 200)
  expect_gt(large$share[2], small$share[2])
})

test_that("detection is bounded in (0, 1) and improves with budget", {
  rate <- matrix(c(0.3, 0.1), nrow = 1)
  d_small <- allocate_branches(rate, budget = 5)$detection
  d_large <- allocate_branches(rate, budget = 50)$detection
  expect_gt(d_large, d_small)
  expect_lt(d_large, 1)
})

test_that("input guards fire", {
  expect_error(allocate_branches(rate = c(0.3, 0.1), budget = 10), "matrix")
  expect_error(allocate_branches(matrix(c(-0.1, 0.2), nrow = 1), budget = 10),
               "non-negative")
  expect_error(allocate_branches(matrix(0.3, nrow = 1), budget = -1),
               "positive")
})
