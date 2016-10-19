library(australian)
library(foreach)

context('Test Upsilon Construction')

test_that('Upsilon has correct properties', {
  theta <- rbind(0.25, 0.25)
  X <- array(rbind(1:10, 1:10), c(2, 1, 10))

  Upsilon <- copula_Upsilon(theta, X)

  # Upsilon should have the right format
  expect_equal(dim(Upsilon), c(2, 2, 10))

  # Upsilon should be a correlation matrix
  foreach(i = 1:dim(Upsilon)[3], .combine = 'c') %do% {
    Upsilon_t <- Upsilon[,,i]

    # Unit diagonals
    expect_equal(diag(Upsilon_t), rep(1, 2))

    # Symmetric
    expect_equal(Upsilon_t[lower.tri(Upsilon_t)],
                 Upsilon_t[upper.tri(Upsilon[,,i])])
  }
})
