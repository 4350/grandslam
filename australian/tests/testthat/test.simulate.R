library(australian)
library(foreach)
library(mvtnorm)

context("Test Simulation")

copula_spec <- function(dynamics, distribution) {
  CopulaSpecification(
    dynamics = do.call(CopulaDynamics, dynamics),
    distribution = do.call(CopulaDistribution, distribution)
  )
}

test_that('Simulate(constant, normal)', {
  Omega <- matrix(0.50, ncol = 2, nrow = 2)
  diag(Omega) <- 1

  spec <- copula_spec(list(alpha = 0, beta = 0, phi = 0, Omega = Omega),
                      list(nu = Inf, gamma = rep(0, 2)))

  registerDoSEQ()
  simulations <- copula_simulate(spec, 100, 50)
  expect_length(simulations, 50)
  expect_equal(dim(simulations[[1]]), c(100, 2))
})

test_that('Simulate(dynamic, normal)', {
  Omega <- matrix(0.50, ncol = 2, nrow = 2)
  diag(Omega) <- 1

  spec <- copula_spec(list(alpha = 0.06, beta = 0.91, phi = 0, Omega = Omega),
                      list(nu = Inf, gamma = rep(0, 2)))

  registerDoSEQ()
  simulations <- copula_simulate(spec, 100, 5, shocks_T = cbind(0, 0), Q_T = Omega)
  expect_length(simulations, 5)
  expect_equal(dim(simulations[[1]]), c(100, 2))
})

test_that('Dynamic simulation without proper data gives warnings', {
  Omega <- matrix(0.50, ncol = 2, nrow = 2)
  diag(Omega) <- 1

  spec <- copula_spec(list(alpha = 0.06, beta = 0.91, phi = 0, Omega = Omega),
                      list(nu = Inf, gamma = rep(0, 2)))

  registerDoSEQ()
  expect_warning(copula_simulate(spec, 1, 1, Q_T = Omega))
  expect_warning(copula_simulate(spec, 1, 1, shocks_T = cbind(0, 0)))
})

test_that('Dynamic simulation with X works', {
  Omega <- matrix(0.50, ncol = 2, nrow = 2)
  diag(Omega) <- 1

  spec <- copula_spec(list(alpha = 0.06, beta = 0.91, phi = 0,
                           theta = rbind(0.50, 0.50), Omega = Omega),
                      list(nu = Inf, gamma = rep(0, 2)))

  X <- array(rbind(1:5, 1:5), c(2, 1, 5))
  registerDoSEQ()
  simulations <- copula_simulate(spec, 5, 1,
                                 Q_T = Omega, shocks_T = cbind(0, 0), X = X)

  expect_length(simulations, 1)
  expect_equal(dim(simulations[[1]]), c(5, 2))
})
