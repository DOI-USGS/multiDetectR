## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6,
                      fig.height = 4)

## ----setup--------------------------------------------------------------------
library(multiDetectR)

## -----------------------------------------------------------------------------
des <- allocate(
  theta  = c(0.3, 0.5, 0.6),   # site, sample, replicate detectabilities
  cost   = c(200, 20, 2),      # cost of a site, a sample, a replicate
  goal   = "detection",
  budget = 5000
)
des

## -----------------------------------------------------------------------------
summary(des)

## -----------------------------------------------------------------------------
allocate(theta = c(0.3, 0.5, 0.6), cost = c(200, 20, 2),
         goal = "detection", P_target = 0.95)

## -----------------------------------------------------------------------------
plot(des)

## -----------------------------------------------------------------------------
allocate(theta = c(0.3, 0.5), cost = c(200, 2),
         goal = "count", budget = 5000)

## -----------------------------------------------------------------------------
enumerate_designs(theta = c(0.3, 0.5, 0.6), cost = c(200, 20, 2),
                  budget = 1500, goal = "detection")

## -----------------------------------------------------------------------------
optimality_identity(allocate(theta = c(0.3, 0.5, 0.6), cost = c(200, 20, 2),
                             goal = "detection", budget = 5000,
                             integer = FALSE))

## -----------------------------------------------------------------------------
## One target, two methods; the first is nine times better.
allocate_branches(rate = matrix(c(0.9, 0.1), nrow = 1), budget = 20)

