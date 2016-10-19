library(australian)

context("Test Copula Fit")

test_that('copula_fit_pars', {
  spec <- CopulaSpecification(
    dynamics = CopulaDynamics(alpha = 0, beta = 0, phi = 0, theta = rbind(0, 0)),
    distribution = CopulaDistribution(nu = Inf, gamma = c(0, 0))
  )

  r <- copula_fit_pars(spec, 'norm', T, F)
  expect_length(r, 3)
  expect_null(r$pars)

  spec@distribution@nu <- 8
  r <- copula_fit_pars(spec, 't', T, F)
  expect_length(r, 3)
  expect_equal(names(r$pars), c('distribution.nu'))
  expect_equivalent(r$pars, c(8))
  expect_length(r$ci, 2)
  expect_equal(dim(r$ui), c(2, 1))

  spec@distribution@gamma <- c(0.25, 0.10)
  r <- copula_fit_pars(spec, 'ghst', T, F)
  expect_length(r, 3)
  expect_equal(names(r$pars), c('distribution.nu',
                                'distribution.gamma1', 'distribution.gamma2'))
  expect_equivalent(r$pars, c(8, 0.25, 0.10))
  expect_length(r$ci, 2 + 2 * 2)
  expect_equal(dim(r$ui), c(2 + 2 * 2, 3))

  spec@dynamics@alpha <- 0.25
  spec@dynamics@beta <- 0.10
  spec@distribution@nu <- Inf
  spec@distribution@gamma <- c(0, 0)
  r <- copula_fit_pars(spec, 'norm', F, F)
  expect_length(r, 3)
  expect_equal(names(r$pars), c('alphabeta.alpha', 'alphabeta.beta'))
  expect_equivalent(r$pars, c(0.25, 0.10))
  expect_length(r$ci, 3)
  expect_equal(dim(r$ui), c(3, 2))

  spec@dynamics@phi <- 0.25
  spec@dynamics@theta <- rbind(0.10, 0.10)
  r <- copula_fit_pars(spec, 'norm', F, T)
  expect_length(r, 3)
  expect_equal(names(r$pars), c('alphabeta.alpha', 'alphabeta.beta',
                                'upsilon.phi', 'upsilon.theta'))
  expect_equivalent(r$pars, c(0.25, 0.10, 0.25, 0.10))
  expect_length(r$ci, 3 + 2)
  expect_equal(dim(r$ui), c(3 + 2, 4))

  # Don't support different theta
  spec@dynamics@theta <- rbind(0.10, 0.20)
  expect_error(copula_fit_pars(spec, 'norm', F, T))
})
