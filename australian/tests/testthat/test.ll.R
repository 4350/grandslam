library(australian)
library(foreach)
library(mvtnorm)

context("Test Correct Log Likelihood")

copula_spec <- function(dynamics, distribution) {
  CopulaSpecification(
    dynamics = do.call(CopulaDynamics, dynamics),
    distribution = do.call(CopulaDistribution, distribution)
  )
}

setup_uniform <- function() {
  sigma <- matrix(0.50, 2, 2)
  diag(sigma) <- 1

  set.seed(403)
  x <- rmvnorm(2000, mean = rep(0, 2), sigma = sigma)
  pnorm(x)
}

expect_ll <- function(dynamics, distribution, ll, X = NULL) {
  spec <- copula_spec(dynamics, distribution)
  u <- setup_uniform()

  registerDoSEQ()
  filtered <- australian::copula_filter(spec, u, X)
  expect_equal(filtered$ll, ll, tolerance = 1e-4)
}

test_that('LL(constant, normal)', {
  expect_ll(
    list(alpha = 0, beta = 0),
    list(nu = Inf, gamma = rep(0, 2)),
    269.775
  )
})

test_that('LL(constant, t(8))', {
  expect_ll(
    list(alpha = 0, beta = 0),
    list(nu = 8, gamma= rep(0, 2)),
    257.9517
  )
})

test_that('LL(constant, ghst(8, c(0.25, 0.25)))', {
  expect_ll(
    list(alpha = 0, beta = 0),
    list(nu = 8, gamma= rep(0.25, 2)),
    256.1900
  )
})

# Dynamic Copulas --------------------------------------------------------

test_that('LL(0.06 / 0.91, normal)', {
  expect_ll(
    list(alpha = 0.06, beta = 0.91),
    list(nu = Inf, gamma = rep(0, 2)),
    247.8804
  )
})

test_that('LL(0.06 / 0.91, t(8))', {
  expect_ll(
    list(alpha = 0.06, beta = 0.91),
    list(nu = 8, gamma = rep(0, 2)),
    240.8095
  )
})

test_that('LL(0.06 / 0.91, ghst(8, c(0.25, 0.25)))', {
  expect_ll(
    list(alpha = 0.06, beta = 0.91),
    list(nu = 8, gamma = rep(0.25, 2)),
    237.8483
  )
})

# Dynamic Copulas with Upsilon -------------------------------------------

time_trend <- function() {
  t <- 1:2000
  array(rbind(t, t), c(2, 1, length(t)))
}

test_that('LL(0.06 / 0.91, normal, with time trend)', {
  X <- time_trend()

  # It's a bit worrisome that this has a higher likelihood than the copula
  # without time trend. But maybe it's the dynamicity too.
  expect_ll(
    list(alpha = 0.06, beta = 0.91, phi = 0.50, theta = rbind(0.50, 0.50)),
    list(nu = Inf, gamma = rep(0, 2)),
    248.5022,
    X = X
  )
})

test_that('LL(0.06 / 0.91, t(8), with time trend)', {
  X <- time_trend()

  # It's a bit worrisome that this has a higher likelihood than the copula
  # without time trend. But maybe it's the dynamicity too.
  expect_ll(
    list(alpha = 0.06, beta = 0.91, phi = 0.50, theta = rbind(0.50, 0.50)),
    list(nu = 8, gamma = rep(0, 2)),
    241.3514,
    X = X
  )
})

test_that('LL(0.06 / 0.91, ghst(8, c(0.25, 0.25)), with time trend)', {
  X <- time_trend()

  # It's a bit worrisome that this has a higher likelihood than the copula
  # without time trend. But maybe it's the dynamicity too.
  expect_ll(
    list(alpha = 0.06, beta = 0.91, phi = 0.50, theta = rbind(0.50, 0.50)),
    list(nu = 8, gamma = rep(0.25, 2)),
    238.5208,
    X = X
  )
})
