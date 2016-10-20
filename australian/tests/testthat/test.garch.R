library(australian)
library(mvtnorm)
library(doParallel)

context("Test GARCH functionality")

test_that('garch_fit works', {
  registerDoSEQ()

  garch <- list(garch_specgen(0, 0),
                garch_specgen(1, 0))

  x <- rmvnorm(100, c(0, 0))
  fits <- garch_fit(garch, x)

  expect_length(fits, 2)
})

test_that('garch_fit works with data.frame', {
  registerDoSEQ()

  garch <- list(garch_specgen(0, 0))
  x <- data.frame(x = rnorm(100))
  fits <- garch_fit(garch, x)
  expect_length(fits, 1)
})
