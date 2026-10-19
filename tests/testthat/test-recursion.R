test_that("detection_recursion matches a direct Bernoulli-cascade simulation", {
  ## The analytic recursion (eqs 1-2) should match the fraction of simulated
  ## cascades in which the target is found, which shares none of its algebra.
  simulate_cascade <- function(theta, n, n_sim = 2e5, seed = 1) {
    set.seed(seed)
    d <- length(theta)
    detected <- logical(n_sim)
    for (s in seq_len(n_sim)) {
      detected[s] <- cascade_present(theta, n, 1)
    }
    mean(detected)
  }
  ## Recursive presence check: at level i, any of the n[i] units is present
  ## (prob theta[i]) AND, if an inner level exists, detection succeeds within it.
  cascade_present <- function(theta, n, i) {
    d <- length(theta)
    for (k in seq_len(n[i])) {
      if (stats::runif(1) < theta[i]) {
        if (i == d) return(TRUE)
        if (cascade_present(theta, n, i + 1)) return(TRUE)
      }
    }
    FALSE
  }

  theta <- c(0.3, 0.5, 0.6)
  n <- c(4, 3, 2)
  analytic <- detection_recursion(theta, n)$P1
  simulated <- simulate_cascade(theta, n, n_sim = 2e5, seed = 42)
  expect_equal(analytic, simulated, tolerance = 0.01)
})

test_that("detection_recursion edges are correct", {
  ## Single level reduces to 1 - (1 - theta)^n.
  expect_equal(detection_recursion(0.4, 3)$P1, 1 - (1 - 0.4)^3)
  ## More replicates never lower detection.
  p2 <- detection_recursion(c(0.3, 0.5), c(5, 2))$P1
  p3 <- detection_recursion(c(0.3, 0.5), c(5, 4))$P1
  expect_gte(p3, p2)
})

test_that("total_cost and unit_cost_full satisfy C = n_1 * C_1", {
  cost <- c(200, 20, 2)
  n <- c(10, 3, 2)
  expect_equal(total_cost(cost, n), n[1] * unit_cost_full(cost, n)[1])
  ## Three-level cost expands as in equation (3).
  expect_equal(total_cost(cost, n),
               cost[1] * n[1] + cost[2] * n[1] * n[2] +
                 cost[3] * n[1] * n[2] * n[3])
})

test_that("nats and prob invert each other and hit known thresholds", {
  expect_equal(prob(nats(0.73)), 0.73)
  expect_equal(nats(0.95), 2.9957, tolerance = 1e-3)
  expect_equal(nats(0.99), 4.6052, tolerance = 1e-3)
  expect_equal(nats(0.999), 6.9078, tolerance = 1e-3)
})

test_that("invalid inputs are rejected", {
  expect_error(detection_recursion(c(0.3, 1.2), c(2, 2)), "open interval")
  expect_error(detection_recursion(c(0.3, 0.5), c(2)), "same length")
  expect_error(total_cost(c(-1, 2), c(2, 2)), "positive")
})
