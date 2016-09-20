#' Estimate a GARCH-COPULA model on Weekly Returns Data a la
#' Christoffersen and Langlois (2013)

# Libraries ----
library(rugarch)
library(sn)
library(ghyp)
library(dplyr)

# Reset Workspace to just Returns Data ----
rm(list = ls())
load("data/ff-weekly.RData")
df <- select(df, -RF)

# Build spec and fit GARCH model ----
garch.spec <- ugarchspec(
  mean.model = list(
    armaOrder = c(3, 0)
  ),
  variance.model = list(
    model = "fGARCH",
    submodel = "NGARCH",
    variance.targeting = TRUE
  ),
  distribution.model = 'sstd'
)

# Fit GARCH to each factor, storing result in a list for each factor,
# where each factor in turn has a list with element "fit" holding the results
# (so we can hold more stuff next to it later)
garch.fit <- apply(
  select(df, -Date), 2,
  function(ts) {
    list(fit = ugarchfit(garch.spec, ts))
  }
)

# Find the eta elements for the copula from the standardized residuals ----
# d = pdf, p = cdf, q = quantile

get.u <- function(ts) {
  residuals.st <- ts$fit@fit$residuals / ts$fit@fit$sigma
  
  # Generate empirical CDF for the standardized residuals and use it to
  # compute union margins (eta)
  #cdf <- ecdf(residuals.st)
  
  # Ugly hack. PS get your shit together Victor
  cdf <- rugarch:::psstd(
    residuals.st,
    skew = ts$fit@fit$coef[['skew']],
    shape = ts$fit@fit$coef[['shape']]
  )
  
  list(
    fit = ts$fit,
    residuals.st = residuals.st,
    u = cdf#(residuals.st)
  )
}

garch.fit <- lapply(garch.fit, get.u)

# ML estimation of dynamic copula ----
# We need to compute the likelihood function for each iteration of parameters
#
# Step 1: Generate Q series
#
# Given an estimate of kappa and lambda (df and skewness), we can construct
# the standardized fractile series z_t. The fractile series needs to be
# standardized; Appendix has details about how to do that (subtract expectation
# divide by standard deviation)
#
# We should then be able to recursively generate Q_t series. Q is estimated
# as the sample correlation of z_t over the entire sample period.

get.epsilon <- function(u, df, skew) {
  epsilon.c.mean <- df / (df - 2) * skew
  
  # The covariance is dervied from the correlation matrix Psi, which has unit
  # diagonals. If we had had normals for our copula "innovations" then this
  # reduces to just rep(1, ncol(u)).
  epsilon.c.variance <-
    sd(
      df / (df - 2) * rep(1, ncol(u)) +
      (2 * df ^ 2 * skew ^ 2) / ((df - 2)^2 * (df - 4))
    )
  
  sapply(seq(ncol(u)), function(i) {
    # Invert the univariate t distribution to get skew t-distributed
    # innovations
    
    # USE THE RIGHT QUANTILE FUNCTION
    epsilon.c <- qst(u[,i], alpha = skew[i], nu = df, method = 2, tol = 1/100000)
    
    # Standardize them so that they have expectation 0 and sd 1 according to
    # the conditional distribtuion they're assumed to be generated from
    (epsilon.c - epsilon.c.mean) / (epsilon.c.sd)
  })
}

get.zstar_bar <- function(zstar, alpha, beta) {
  # Recursively build the diagonal Gamma by going over the series
  
  # Add an additional observation
  zstar <- rbind(rep(0, ncol(zstar)), zstar)
  
  # The normalized zstars
  zstar_bar <- zstar * NA
  
  # The diagonals of gamma
  gamma <- zstar * NA
  gamma[1,] <- 1
  
  for (t in 2:nrow(zstar)) {
    # Compute gamma diagonal for this period
    gamma[t,] <-
      (1 - alpha - beta) +
      beta * gamma[t - 1,] +
      alpha * (zstar[t - 1,] / gamma[t - 1,]) ^ 2
    
    # Compute zstar_bar for this period by normalizing by gamma of this period
    zstar_bar[t,] <- zstar[t,] / sqrt(gamma[t,])
  }
  
  zstar_bar[-1,]
}

get.Omega <- function(zstar_bar) {
  # Sum the cross products of zstar_bar over time
  Omega <- matrix(0, ncol(zstar_bar), ncol(zstar_bar))
  for (t in 1:nrow(zstar_bar)) {
    Omega = Omega + (zstar_bar[t,] %*% t(zstar_bar[t,]))
  }
  
  # Divide by T to get average
  Omega / nrow(zstar_bar)
}

get.Gamma <- function(zstar_bar, Omega, alpha, beta) {
  N <- ncol(zstar_bar)
  T <- nrow(zstar_bar)
  
  # Create observations for t0
  Gamma <- array(dim = c(N, N, T + 1))
  zstar_bar <- rbind(rep(0, N), zstar_bar)
  
  # Initialize t0 Gamma to Omega (alternatively: diag(1, N, N))
  Gamma[,,1] <- Omega
  
  for (t in 2:(T + 1)) {
    Gamma[,,t] <-
      (1 - alpha - beta) * Omega +
      beta * Gamma[,,t-1] +
      alpha * zstar_bar[t-1,] %*% t(zstar_bar[t-1,])
  }
  
  Gamma[,,-1]
}

normalize.Gamma <- function(Gamma) {
  N <- dim(Gamma)[1]
  T <- dim(Gamma)[3]
  
  Psi <- Gamma * NA
  for (t in seq(T)) {
    for (i in seq(N)) {
      for (j in seq(N)) {
        Psi[i,j,t] <- Gamma[i,j,t] / sqrt(Gamma[i,i,t] * Gamma[j,j,t])
      }
    }
  }
  Psi
}

logLikelihood <- function(eta, df, skew, alpha, beta) {
  # Construct correlation series according to Christoffersen
  zstar <- get.zstar(eta, df, skew)
  zstar_bar <- get.zstar_bar(zstar, alpha, beta)
  Omega <- get.Omega(zstar_bar)
  Gamma <- get.Gamma(zstar_bar, Omega, alpha, beta)
  Correlation <- normalize.Gamma(Gamma)
  
  # zstar is our x
  marginalLL <- sapply(
    seq(ncol(zstar)),
    function(i) dst(zstar[,i], alpha = skew[i], nu = df, log = TRUE)
  )
  marginalLL <- rowSums(marginalLL)
  
  # Join likelihood depends on time-varying correlation matrix
  jointLL <- sapply(seq(T), function(t) {
    dmst(zstar[t,], Omega = Correlation[,,t], alpha = skew, nu = df, log = TRUE)
  })
  
  # Apply Sklar's theorem
  sum(jointLL - marginalLL)
}

# Optimized LL function with a single skewness parameter
optimize1Fn <- function(params, eta) {
  df <- params[1]
  skew <- rep(params[2], ncol(eta)) # Repeat!
  alpha <- params[3]
  beta <- params[4]
  
  print(params)
  
  if (alpha < 0 | beta < 0 | alpha + beta >= 1) {
    return(100000)
  }
  
  -logLikelihood(eta, df = df, skew = skew, alpha = alpha, beta = beta)
}

# Do Optimization ----
# Input data
params <- c(5, 1, 0.03, 0.96)
eta <- sapply(garch.fit, function(f) f$eta)

param.nlm <- nlm(
  optimize1Fn,
  p = params,
  eta = eta,
  print.level = 2,
  iterlim = 500
)

# Get Correlation ----
get.series.Correlation <- function(eta, df, skew, alpha, beta) {
  zstar <- get.zstar(eta, df, skew)
  zstar_bar <- get.zstar_bar(zstar, alpha, beta)
  Omega <- get.Omega(zstar_bar)
  Gamma <- get.Gamma(zstar_bar, Omega, alpha, beta)
  Correlation <- normalize.Gamma(Gamma)
  Correlation
}

corr <- get.series.Correlation(
  eta,
  df = df,
  skew = rep(param.nlm$estimate[1], ncol(eta)),
  alpha = param.nlm$estimate[2],
  beta = param.nlm$estimate[3]
)
plot(corr[1,2,], type = 'l')

# Other ----

# Initial parameters: df, skew, alpha, beta
# params <- c(5, 1, 0.03, 0.96)
# 
# params.optim <- optim(
#   params,
#   optimize1Fn,
#   gr = NULL,
#   eta,
#   lower = c(3, -Inf, 0, 0),
#   upper = c(Inf, Inf, 0.9999, 0.9999),
#   control = list(
#     trace = 10,
#     maxit = 2
#   ),
#   method = 'SANN'
# )
# 
# # Constrained optimization ----
# ui <- rbind(
#   c( 1,  0,  0,  0),
#   c( 0,  0,  1,  0),
#   c( 0,  0,  0,  1),
#   c( 0,  0, -1, -1)
# )
# ci <- rbind(
#   3,
#   0,
#   0,
#   -0.9999
# )
# params.constrOptim <- constrOptim(
#   params,
#   optimize1Fn,
#   grad = NULL,
#   ui = ui,
#   ci = ci,
#   control = list(
#     trace = 5,
#     maxit = 1,
#     fnscale = -1
#   ),
#   eta = eta
# )
