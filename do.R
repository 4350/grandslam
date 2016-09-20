#' Estimate a GARCH-COPULA model on Weekly Returns Data a la
#' Christoffersen and Langlois (2013)

# Libraries ----
library(rugarch)
library(sn)
library(ghyp)
library(dplyr)
library(parallel)

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

# Dynamic Copula Estimation ----

source('func/dynamicCopula.R')

#' Optimize log-likelihood with a single skewness parameter
#'
#' @param params (df, skew, alpha, beta)
#' @param u 
#'
#' @return negative log likelihood
optimize1Fn <- function(params, u) {
  df <- params[1]
  skew <- rep(params[2], ncol(u))
  alpha <- params[3]
  beta <- params[4]
  
  print(params)
  
  -logLikelihood(u, df, skew, alpha, beta)
}

#' Optimize conditional log-likelihood with a single skewness parameter
#' 
#' @param 
optimize1FnCond <- function(params, u) {
  df <- params[1]
  skew <- rep(params[2], ncol(u))
  alpha <- params[3]
  beta <- params[4]
  
  print(params)
  -condLogLikelihood(u, df, skew, alpha, beta)
}

# Do Optimization ----
# Input data
params <- c(5, 1, 0.03, 0.96)
u <- sapply(garch.fit, function(f) f$u)

param.constrOptim <- constrOptim(
  theta = params,
  optimize1FnCond,
  grad = NULL,
  u = u,
  ui = rbind(
    c( 1,  0,  0,  0),
    c(-1,  0,  0,  0),
    c( 0,  0,  1,  0),
    c( 0,  0,  0,  1),
    c( 0,  0, -1, -1)
  ),
  ci = rbind(
    4,
    -20,
    0,
    0,
    -0.9999
  ),
  control = list(
    trace = 5,
    maxit = 10
  )
)

# Get Correlation ----
# get.series.Correlation <- function(eta, df, skew, alpha, beta) {
#   zstar <- get.zstar(eta, df, skew)
#   zstar_bar <- get.zstar_bar(zstar, alpha, beta)
#   Omega <- get.Omega(zstar_bar)
#   Gamma <- get.Gamma(zstar_bar, Omega, alpha, beta)
#   Correlation <- normalize.Gamma(Gamma)
#   Correlation
# }
# 
# corr <- get.series.Correlation(
#   eta,
#   df = df,
#   skew = rep(param.nlm$estimate[1], ncol(eta)),
#   alpha = param.nlm$estimate[2],
#   beta = param.nlm$estimate[3]
# )
# plot(corr[1,2,], type = 'l')

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
