#' Estimates threshold and rolling corrs.
#' Currently implementing daily return series. Has previously run weekly
#' series for threshold correlation. The thCorrelationList is saved in derived
#' data for further use in simulated threshold graphs.
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
thCorrList = list()

# Loop for each value factors to get threshold correlation for all factors
# and save to the list created above

thCorrList <- th_corr(df.estim %>% select(-Date), 0)

# Bind to one df for plot
plotdf <- bind_rows(thCorrList$HML, thCorrList$RMW, thCorrList$CMA)

# Do plots for every factor's data frame of threshold correlations vs all other factors
plots <- correlations.plot.threshold(plotdf)+
  ggtitle("Threshold correlations of daily log returns (95% confidence bounds)")

# Arrange and save plots
g <- grid.arrange(
  plots
)

ggsave(file = 'output/ThresholdCorrelationsDaily2.jpeg', g, width = 16.6, height = 11.7, units = 'in')


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

# Bind to one df for plot
plotdf <- bind_rows(rollingCorrelationList$HML, rollingCorrelationList$RMW, rollingCorrelationList$CMA)

plots <- correlations.plot.rolling(plotdf)+
  ggtitle("Rolling 250-day correlations (95% confidence bounds)")

# Arrange and save plots
g <- grid.arrange(
  plots
)

ggsave(
  file = 'output/rollingCorrelations/250daily2.jpeg',
  g,
  width = 16.6,
  height = 11.7,
  units = 'in'
)