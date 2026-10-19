test_that("A = B at the continuous optimum (Proposition S1.8)", {
  set.seed(11)
  for (k in 1:40) {
    d <- sample(2:5, 1)
    theta <- stats::runif(d, 0.05, 0.95)
    cost <- 10^stats::runif(d, 0, 3)
    des <- allocate(theta, cost, goal = "detection", budget = 1e6,
                    integer = FALSE)
    id <- optimality_identity(des)
    expect_lt(max(abs(id$diff)), 1e-5,
              label = sprintf("max |A - B| for d = %d", d))
  }
})

test_that("optimality_identity returns one row per level with contents", {
  des <- allocate(c(0.3, 0.5, 0.6), c(200, 20, 2), goal = "detection",
                  budget = 5000, integer = FALSE)
  id <- optimality_identity(des)
  expect_equal(nrow(id), 2L)
  expect_named(id, c("level", "A", "B", "diff"))
})

test_that("single-level design has no identity rows", {
  des <- allocate(0.5, 100, goal = "detection", budget = 1000)
  expect_equal(nrow(optimality_identity(des)), 0L)
})

test_that("integer design shows small but nonzero A - B", {
  des <- allocate(c(0.3, 0.5, 0.6), c(200, 20, 2), goal = "detection",
                  budget = 5000, integer = TRUE)
  id <- optimality_identity(des)
  ## Rounding perturbs the identity, but only modestly.
  expect_lt(max(abs(id$diff)), 0.5)
})
