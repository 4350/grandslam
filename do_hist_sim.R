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
n.roll = H

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


#' Function to calculate VaR, using volatility updated 
#' standardized returns, see Hull & White (1998)
#' 
#' @param h Index in the forecast period, h<H
#' @param q Lower tail probability of interest
#' @param std.ret Consolidated matrix of standardized returns (T+H)xN, both in and out of sample
#' @param sigma.fc Consolidated matrix of conditional volatilities (T+H)*N, both in and out of sample
#' @param T End of estimation period
#' 
#' @return var.star Vector length N. 1-step ahead VaR for forecast index h
volatilityupdate <- function(h, q, std.ret, sigma.fc, T) {
  
  # Standard returns up to last period
  std.ret <- std.ret[1:(T+h-1),]
  r.star <- apply(std.ret, 1, function(r, h) {
    # Multiply each time t standardized rets by next forecast of vol
    r * sigma.fc[(T+h),]
  },
  h = h)
  
  # Turn it right
  r.star <- t(r.star)
  
  # Calculate VaRs
  var.star <- apply(r.star, 2, function(r, q) {
    quantile(r, probs = q)
  },
  q = q
  )
  
  var.star
}

#' Historical simulation using volatility updated returns
#' 
#' @param q Lower CDF probability of interest
#' @param std.ret Consolidated matrix of standardized returns (T+H)xN, both in and out of sample
#' @param sigma.fc Consolidated matrix of conditional volatilities (T+H)*N, both in and out of sample
#' @param T End of estimation period
#' 
#' @return HS.VaR Historically simulated VaR
hist.sim <- function(q, std.ret, sigma.fc, T) {
  H = nrow(std.ret) - T - 1
  HS.VaR = sapply(seq(1:H), function(h, q, std.ret, sigma.fc, T) {
    volatilityupdate(h, q, std.ret, sigma.fc, T)
  },
  q = q,
  std.ret = std.ret,
  sigma.fc = sigma.fc, 
  T = T
  )
  HS.VaR = t(HS.VaR)
}

library(tidyr)
library(ggplot2)

# ugly stuff to try it quickly :)
a <- hist.sim(0.01, std.ret, sigma.fc, T)
a <- as.data.frame(a)

a <- gather(a, 'factor','value',1:6)
a$h <- rep(seq(1,H), 6)
a$order <- factor(a$factor, names(kGARCHModels))

b <- hist.sim(0.05, std.ret, sigma.fc, T)
b <- as.data.frame(b)

b <- gather(b, 'factor','value',1:6)
b$h <- rep(seq(1,H), 6)
b$order <- factor(b$factor, names(kGARCHModels))

emp <- gather(df[(T+1):nrow(df), ], 'factor','value', 1:6)
emp$h <- rep(seq(1,H), 6)
emp$order <- factor(emp$factor, names(kGARCHModels))

ggplot(a, aes(x = h, y = value, group = factor))+
  geom_line()+
  geom_line(aes(color = 'blue'), data = b)+
  geom_line(aes(color = factor), data = emp)+
  facet_grid(order ~ .)

