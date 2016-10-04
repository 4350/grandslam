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
n.roll = H - 1

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

# # Consolidate series of sigma in and out of sample ----
sigma <- rbind(
  # Get conditional sigma in-sample
  sapply(garch.fits,
         function(fFit) {
           sigma = fFit@fit$sigma
         }),
  # Get conditional sigma "out-of-sample"
  sapply(fcList,
         function(ffc) {
           sigma = ffc@forecast$sigmaFor
         })
)


#' Function to calculate VaR, using volatility updated 
#' standardized returns, see Hull & White (1998)
#' 
#' @param h Index in the forecast period, h<=H
#' @param q Lower tail probability of interest
#' @param df 
#' @param sigma Consolidated matrix of conditional volatilities (T+H)*N, both in and out of sample
#' @param T End of estimation period
#' 
#' @return var.star Vector length N. 1-step ahead VaR for forecast index h
volatilityupdate <- function(h, q, df, sigma, T) {
  
  # Standard returns up to last period
  std.ret <- df[1:(T+h-1), ] / sigma[1:(T+h-1), ]
  r.star <- apply(std.ret, 1, function(r, h) {
    # Multiply each time t standardized rets by next forecast of vol
    r * sigma[(T+h),]
  },
  h = h)
  
  # Turn it right
  r.star <- t(r.star)
  
  #Calculate VaRs
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
#' @return HS.VaR df of historically simulated VaR
hist.sim <- function(q, df, sigma.fc, T) {
  H = nrow(df) - T
  HS.VaR = sapply(seq(1:H), function(h, q, df, sigma, T) {
    volatilityupdate(h, q, df, sigma, T)
  },
  q = q,
  df = df,
  sigma = sigma, 
  T = T
  )
  HS.VaR = data.frame(t(HS.VaR))
}

# Libraries for plot ----
library(tidyr)
library(ggplot2)
library(extrafont)

# Tidy data frame return for plot

hist.sim.tidy <- function(q, df, sigma, T, kGARCHModels) {
  H = nrow(df) - T
  N = ncol(df)
  # Run simulation
  out.df <- hist.sim(q, df, sigma, T)
  # Long format data
  out.df <- gather(out.df, 'factor', 'value', 1:N)
  # Add horizon "x values"
  out.df$h <- rep(seq(1,H), N)
  # Order data set for plot
  out.df$order <- factor(out.df$factor, names(kGARCHModels))
  out.df
}

# Get the stuff we want for plot, including empirical

var05 <- hist.sim.tidy(.05, df, sigma, T, kGARCHModels)
var01 <- hist.sim.tidy(.01, df, sigma, T, kGARCHModels)

emp <- gather(df[(T+1):nrow(df), ], 'factor','value', 1:6)
emp$h <- rep(seq(1,H), ncol(df))
emp$order <- factor(emp$factor, names(kGARCHModels))

# Plot that
ggplot(emp, aes(x = h, y = value, group = factor))+
  geom_line(aes(color = 'Realized return'))+
  geom_line(aes(color = '5% HS-HW VaR'), data = var05)+
  geom_line(aes(color = '1% HS-HW VaR'), data = var01)+
  theme_Publication()+
  facet_grid(. ~ order)

# Evalute with tests

var01$violate <- emp$value < var01$value
var05$violate <- emp$value < var05$value
