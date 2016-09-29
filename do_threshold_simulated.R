#' Estimates threshold corrs from simulated data.
#' Uses weekly data of empirical series previously estimated
#' in do_correlation.R code and stored in derived data.
#' 
#' This file is waiting for data from distributions

# Load libraries and functions ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(zoo)
library(scales)
library(devtools)
load_all('wimbledon')

# Reset workspace and load data for each distribution ----
rm(list = ls())
load('data/derived/weekly-estim.RData')

# DATA SHOULD BE LONG FORMAT, E.G. (2479*BOOTRUNS)xN
# load(normaldata)
# load(tdata)
# load(skewtdata)

df.norm <- rbind(df.estim, df.estim) %>% select(-Date)

# Get the correlation series for each distribution
normThCorrelationList <- sim_th_corr(df.norm)
#tThCorrelationListst <- sim_th_corr(df.t)
#skewtThCorrelationListst <- sim_th_corr(df.skewt)

# Save these
# save(normThCorrelationList, file = 'data/derived/simulation/normThCorrelationList.Rdata')
# save(tThCorrelationList, file = 'data/derived/simulation/tThCorrelationList.Rdata')
# save(skewtThCorrelationList, file = 'data/derived/simulation/skewtThCorrelationList.Rdata')

# Load factor lists from threshold and simulations separately
load('data/derived/simulation/normThCorrelationList.Rdata')
#load('data/derived/simulation/tThCorrelationList.Rdata')
#load('data/derived/simulation/skewtThCorrelationList.Rdata')
load('data/derived/weekly-thCorrelationList.Rdata')

# Set up factor groups ----
factors.value = c("HML", "RMW", "CMA")
factors.all   = c("Mkt.RF", "HML", "SMB", "Mom", "RMW", "CMA")

# Set up factor list that consolidates
simThCorrList <- list()

for (values in factors.value) {
  out.df <- data.frame(
    norm = normThCorrelationList[[values]]$simvalues * 1.2,
    #t = tThCorrelationList[[values]]$simvalues,
    #skewt = skewtThCorrelationList[[values]]$simvalues,
    thCorrelationList[[values]]
  )
  simThCorrList[[values]] <- out.df
}

# Bind to one df for plot
plotdf <- bind_rows(simThCorrList$HML, simThCorrList$RMW, simThCorrList$CMA)

plots <- correlations.plot.threshold.sim(plotdf)

# Arrange and save plots
g <- grid.arrange(
  plots,
    top = textGrob("Threshold correlations with simulated distributions", gp = gpar(fontsize = 15))
)

ggsave(file = 'output/ThresholdCorrelationsSimulated.jpeg', g, width = 16.6, height = 11.7, units = 'in')
