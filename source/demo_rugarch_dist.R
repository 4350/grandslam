#' Demo of our parametrization round-trip
#' 

library(rugarch)
library(ggplot2)
library(ghyp)
library(devtools)
load_all('wimbledon')

rm(list = ls())

# rugarch parametrization (location and scale invariant)
skew <- 0.50
shape <- 8

p <- seq(0.01, 0.99, by = 0.01)

x.ruga <- rugarch:::qsghst(
  p,
  shape = shape,
  skew = skew
)
x.ghyp <- garch.qghyp.rugarch(
  p,
  shape = shape,
  skew = skew
)

# x.ghyp and x.ruga are the same!
