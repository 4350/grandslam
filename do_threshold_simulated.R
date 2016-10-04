#' Estimates threshold corrs from simulated data and plots
#' this together with estimated threshold corrs on empirical data
#' 
#' All threshold correlation applied to (GARCH) residuals
#' 

# Load libraries and functions ----
library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(extrafont)
library(zoo)
library(scales)
library(devtools)
load_all('wimbledon')

# Reset workspace and load data for each distribution ----
rm(list = ls())
load('data/derived/garch_res.RData')

# Set up factor groups ----
factors.value = c("HML", "RMW", "CMA")
factors.all   = c("Mkt.RF", "HML", "SMB", "Mom", "RMW", "CMA")

# DATA SHOULD BE LONG FORMAT, E.G. (2479*BOOTRUNS)xN
df.empirical <- df.res %>% select(-Date)
rm(df.res)

# load(normdata)
# load(tdata)

# load(skewtdata)
df.skewt <- read.table('data/derived/simulation/ghskt/1b813702-897c-11e6-bb9c-9dfbf344bd3b.csv', header = TRUE, sep = ',')
df.skewt <- df.skewt[,14:19]
colnames(df.skewt) <- factors.all

# Get the correlation series for each distribution ----

empiricalThCorrList <- th_corr(df.empirical, 0)
#normThCorrList <- th_corr(df.norm)
#tThCorrList <- th_corr(df.t)
skewtThCorrList <- th_corr(df.skewt, 1)

# Set up factor list and consolidate
allThCorrList <- list()

for (values in factors.value) {
  out.df <- data.frame(
    # The empirical imports coef, lb, ub, quantiles and ordering of factors
    empiricalThCorrList[[values]],
    
    # All other distributions only imports coefficients
    #norm = normThCorrList[[values]]$simcoef,
    #t = tThCorrList[[values]]$simcoef,
    skewt = skewtThCorrList[[values]]$simcoef
  )
  allThCorrList[[values]] <- out.df
}

# Bind to one df for plot
plotdf <- bind_rows(allThCorrList$HML, allThCorrList$RMW, allThCorrList$CMA)

plots <- correlations.plot.threshold.sim(plotdf)+
  ggtitle("Threshold correlations of weekly data sets (95% confidence bounds)")

# Arrange and save plots
g <- grid.arrange(
  plots
)

ggsave(file = 'output/ThresholdCorrelationsSimulatedResiduals.jpeg', g, width = 16.6, height = 11.7, units = 'in')
