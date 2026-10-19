test_that("round-in-place matches enumeration once the budget is large", {
  ## Manuscript S2: the two agree once the budget buys ~10+ outermost units.
  theta <- c(0.3, 0.5, 0.6)
  cost <- c(200, 20, 2)
  budget <- 1e5                       # buys many outermost units
  rounded <- allocate(theta, cost, goal = "detection", budget = budget)
  enum <- enumerate_designs(theta, cost, budget = budget, goal = "detection")
  ## Detection should match closely; counts should be identical or adjacent.
  expect_equal(rounded$P1, enum$P1, tolerance = 1e-3)
})

test_that("enumeration is at least as good as rounding at a small budget", {
  ## Where rounding is unreliable, the enumerated optimum cannot be worse.
  theta <- c(0.3, 0.5, 0.6)
  cost <- c(200, 20, 2)
  budget <- 1200                      # only a few outermost units
  rounded <- allocate(theta, cost, goal = "detection", budget = budget)
  enum <- enumerate_designs(theta, cost, budget = budget, goal = "detection")
  expect_gte(enum$P1 + 1e-9, rounded$P1)
})

test_that("enumerated design respects the budget", {
  theta <- c(0.4, 0.6)
  cost <- c(100, 5)
  enum <- enumerate_designs(theta, cost, budget = 1000, goal = "detection")
  expect_lte(enum$C, 1000)
  expect_true(enum$enumerated)
})

test_that("enumerate_designs errors when nothing fits", {
  expect_error(
    enumerate_designs(c(0.3, 0.5), c(1e6, 1e5), budget = 10,
                      goal = "detection"),
    "feasible|exceeds"
  )
})
