#'
#'


# Libraries ----
library(rugarch)
library(parallel)
library(devtools)
library(dplyr)
load_all('wimbledon')

# Data ----
rm(list = ls())
load('data/derived/model_GARCH.RData')
load('data/derived/weekly-full.RData')
load('data/derived/weekly-estim.RData')

# Get out of sample length
T = nrow(df.estim)
H = nrow(df) - T

rm(df.estim)

df <- df %>% select(-Date)

# Best models, as determined before
kGARCHModels <- list(
  Mkt.RF = garch.specgen(0, 0),
  HML = garch.specgen(1, 1),
  SMB = garch.specgen(1, 1),
  Mom = garch.specgen(1, 0),
  RMW = garch.specgen(1, 0),
  CMA = garch.specgen(1, 0) 
)

#' Function to estimate GARCH models with full data set
estimate.garch <- function(df, nOOS) {
  garch <- lapply(names(kGARCHModels), function(name, nOOS) {
    model <- kGARCHModels[[name]]
    data <- as.data.frame(df[, name])
    ugarchfit(model, data, out.sample = nOOS, solver = 'hybrid')
  },
  nOOS = nOOS)
  names(garch) <- names(kGARCHModels)
  
  garch
}

# Get GARCH fits that preserve the OOS data
garch.fits <- estimate.garch(df, H)

# Make the rolling 1-step ahead forecasts of sigma and the series
n.roll = 30

fcList = list()

fcList <- lapply(garch.fits, 
                 function(fFit, n.roll) {
                   ugarchforecast(
                     fitORspec = fFit,
                     n.ahead = 1,
                     n.roll = n.roll
                   )
                 },
                 n.roll = n.roll)

# Consolidate series of standardized returns and sigma in and out of sample ----
std.ret <- rbind(
  # Get conditional standardized returns in-sample
  sapply(garch.fits,
         function(fFit) {
           std.ret = (fFit@fit$z/100) / fFit@fit$sigma 
         }),
  # Get conditional standardized returns "out-of-sample"
  sapply(fcList,
         function(ffc) {
           std.ret = ffc@forecast$seriesFor / ffc@forecast$sigmaFor
         })
)

sigma.fc <- rbind(
  # Get conditional sigma in-sample
  sapply(garch.fits,
         function(fFit) {
           sigma.fc = fFit@fit$sigma 
         }),
  # Get conditional sigma "out-of-sample"
  sapply(fcList,
         function(ffc) {
           sigma.fc = ffc@forecast$sigmaFor
         })
)


#' Function to update the standardized returns with latest forecast of vol
#' See Hull & White (1998)
#' @param h Index in the forecast period, h<H
#' @param std.ret Consolidated matrix of standardized returns (T+H)xN, both in and out of sample
#' @param sigma.fc Consolidated matrix of conditional volatilities (T+H)*N, both in and out of sample
#' @param T End of estimation period
#' 
#' @return r.star Volatility updated returns (T+h-1)xN
volatilityupdate <- function(h, std.ret, sigma.fc, T) {
  # Standard returns up to last period
  std.ret <- std.ret[1:(T+h-1),]
  r.star <- apply(std.ret, 1, function(r, h) {
    # Multiply each time t standardized rets by next forecast of vol
    r * sigma.fc[T+h,]
  },
  h = h)
  r.star <- t(r.star)
}

# Set probaiblities for VaR
kQs <- c(0.01, 0.05)

#' Historical simulation using volatility updated returns
#' 
#' @param q Lower CDF probability of interest
#' @param std.ret Consolidated matrix of standardized returns (T+H)xN, both in and out of sample
#' @param sigma.fc Consolidated matrix of conditional volatilities (T+H)*N, both in and out of sample
#' @param T End of estimation period
#' 
#' @return HS.VaR Historically simulated VaR
hist.sim <- function(q, std.ret, sigma.fc, T) {
  H = T - nrow(std.ret)
  
  # Get the vol updated return matrix
  r.star = volatilityupdate(1, std.ret, sigma.fc, T)
  # Calculate the VaR from empirical CDF of vol updated data
  var.star <- apply(r.star, 2, function(r, q) {
    quantile(r, probs = q)
  },
  q = q
  )
}
