#' Estimates threshold and rolling corrs
#' Saves graph grids
#' 

# Libraries ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(zoo)
library(scales)
library(devtools)
# Reset workspace and load return data ----
rm(list = ls())
load('data/derived/daily-estim.RData')

# Load functions ----
load_all('wimbledon')

# Set up factor groups ----
factors.value = c("HML", "RMW", "CMA")
factors.all   = c("Mkt.RF", "HML", "SMB", "Mom", "RMW", "CMA")

# Threshold correlations ----

# Create list to hold correlation sets
thCorrelationList = list()

# Loop for each value factors to get threshold correlation for all factors
# and save to the list created above

for (value in factors.value) {
  
  # Get three statistics sets
correlations <- correlations.factor.apply(factors.all, value, 'coef', df.estim, 1)
  lb <- correlations.factor.apply(factors.all, value, 'lb', df.estim, 1)
  ub <- correlations.factor.apply(factors.all, value, 'ub', df.estim, 1)
  
  # Consolidate and order data frame
  out.df <- data.frame(qs = seq(0.10, 0.90, by = 0.01), correlations, lb = lb[, 'value'], ub = ub[, 'value'])
  out.df$order <- factor(out.df$factor, levels = factors.all)
  
  # Save to common object that holds all factors
  thCorrelationList[[value]] <- out.df
}

# Do plots for every factor's data frame of threshold correlations vs all other factors
plots <- lapply(thCorrelationList, correlations.plot.threshold)

# Arrange and save plots
g <- grid.arrange(
  plots$HML +
    ylab('HML')+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.title.x=element_blank(), axis.ticks.x = element_blank()),
  plots$RMW +
    ylab('RMW')+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.title.x=element_blank(), axis.ticks.x = element_blank()),
  plots$CMA +
    ylab('CMA') +
    xlab('Quantiles'),
  top = textGrob("Threshold correlations on daily data (95% confidence bounds)", gp = gpar(fontsize = 15))
)

ggsave(file = 'output/ThresholdCorrelationsDaily.jpeg', g, width = 8.3, height = 11.7, units = 'in')

# Rolling correlations ----

# Create list to hold correlation sets
rollingCorrelationList = list()

for (value in factors.value) {
  
  # Get three statistics sets
  correlations <- correlations.factor.apply(factors.all, value, 'coef', df.estim, 0)
  lb <- correlations.factor.apply(factors.all, value, 'lb', df.estim, 0)
  ub <- correlations.factor.apply(factors.all, value, 'ub', df.estim, 0)
  
  # Consolidate and order data frame
  out.df <- data.frame(Date = tail(df.estim$Date, (dim(df.estim)[1]-250+1)), 
                       correlations,
                       lb = lb[, 'value'],
                       ub = ub[, 'value']
  )
  out.df$order <- factor(out.df$factor, levels = factors.all)
  
  # Save to common object that holds all factors
  rollingCorrelationList[[value]] <- out.df
  
}

# Do plots for every factor's data frame of threshold correlations vs all other factors
plots <- lapply(rollingCorrelationList, correlations.plot.rolling)

# Arrange and save plots
g <- grid.arrange(
  plots$HML +
    ylab('HML')+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.title.x=element_blank(), axis.ticks.x = element_blank()),
  plots$RMW +
    ylab('RMW')+
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.title.x=element_blank(), axis.ticks.x = element_blank()),
  plots$CMA +
    ylab('CMA') +
    xlab(''),
  top = textGrob("Rolling 250-day correlations (95% confidence bounds)", gp = gpar(fontsize = 15))
)

ggsave(
  file = 'output/rollingCorrelations/250daily.jpeg',
  g,
  width = 16.6,
  height = 11.7,
  units = 'in'
)