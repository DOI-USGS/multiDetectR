test_that("single-row assemblage reproduces single-target allocate()", {
  theta <- c(0.3, 0.5, 0.6)
  cost <- c(200, 20, 2)
  single <- allocate(theta, cost, goal = "detection", budget = 5000,
                     integer = FALSE)
  asm <- allocate_assemblage(matrix(theta, nrow = 1), cost, budget = 5000,
                             integer = FALSE)
  ## Same shared inner design and detection (to optimization tolerance).
  expect_equal(asm$targets$P1[1], single$P1, tolerance = 1e-3)
  expect_equal(asm$n[-1], single$n[-1], tolerance = 0.05)
})

test_that("identical targets give the single-target design regardless of weights", {
  theta1 <- c(0.3, 0.5, 0.6)
  cost <- c(200, 20, 2)
  single <- allocate(theta1, cost, goal = "detection", budget = 5000,
                     integer = FALSE)
  theta <- rbind(theta1, theta1, theta1)
  asm <- allocate_assemblage(theta, cost, budget = 5000,
                             weight = c(1, 5, 2), integer = FALSE)
  expect_equal(asm$targets$P1[1], single$P1, tolerance = 1e-3)
})

test_that("the elusive target binds the design", {
  theta <- rbind(easy = c(0.6, 0.7, 0.8), hard = c(0.15, 0.25, 0.35))
  cost <- c(200, 20, 2)
  asm <- allocate_assemblage(theta, cost, budget = 5000)
  ## The harder target has the larger miss-weight and should be the binding one.
  expect_equal(asm$binding, 2L)
  expect_lt(asm$targets$P1[2], asm$targets$P1[1])
})

test_that("assemblage design beats the average-target design on weighted detection", {
  theta <- rbind(c(0.6, 0.7, 0.8), c(0.15, 0.25, 0.35))
  cost <- c(200, 20, 2)
  w <- c(0.5, 0.5)
  asm <- allocate_assemblage(theta, cost, budget = 5000, weight = w,
                             integer = FALSE)

  ## Design built for the average target, evaluated on the true targets.
  avg_theta <- colMeans(theta)
  avg_design <- allocate(avg_theta, cost, goal = "detection", budget = 5000,
                         integer = FALSE)
  n_avg <- avg_design$n
  p1_avg <- vapply(1:2, function(t) {
    detection_recursion(theta[t, ], n_avg)$P1
  }, numeric(1))
  D_avg <- sum(w * p1_avg)

  expect_gte(asm$D + 1e-6, D_avg)
})

test_that("upweighting a hard target shifts the design toward it", {
  theta <- rbind(easy = c(0.7, 0.8, 0.9), hard = c(0.1, 0.2, 0.3))
  cost <- c(200, 20, 2)
  equal <- allocate_assemblage(theta, cost, budget = 5000,
                               weight = c(1, 1), integer = FALSE)
  hardwt <- allocate_assemblage(theta, cost, budget = 5000,
                                weight = c(1, 10), integer = FALSE)
  ## Caring more about the hard target should not lower its detection.
  expect_gte(hardwt$targets$P1[2] + 1e-6, equal$targets$P1[2])
})

test_that("D_target path reaches the required weighted detection", {
  theta <- rbind(c(0.4, 0.6), c(0.3, 0.5))
  cost <- c(100, 5)
  asm <- allocate_assemblage(theta, cost, D_target = 0.8)
  expect_gte(asm$D, 0.8)
})

test_that("input guards fire", {
  expect_error(allocate_assemblage(matrix(c(0.3, 1.1), nrow = 1), c(200, 2),
                                   budget = 1000), "in \\(0, 1\\)")
  expect_error(allocate_assemblage(matrix(c(0.3, 0.5), nrow = 1), c(200),
                                   budget = 1000), "length")
  expect_error(allocate_assemblage(matrix(c(0.3, 0.5), nrow = 1), c(200, 2)),
               "exactly one")
  expect_error(allocate_assemblage(matrix(c(0.3, 0.5), nrow = 1), c(200, 2),
                                   budget = 1000, weight = c(1, 1)),
               "one per target")
})
