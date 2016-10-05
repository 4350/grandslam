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
load('data/derived/garch_stdres.RData')

# Set up factor groups ----
factors.value = c("HML", "RMW", "CMA")
factors.all   = c("Mkt.RF", "HML", "SMB", "Mom", "RMW", "CMA")

# DATA SHOULD BE LONG FORMAT, E.G. (2479*BOOTRUNS)xN
df.empirical <- df.stdres %>% select(-Date)
rm(df.res)

# load(normdata)
# load(tdata)

# load(skewtdata)
# Load me some constant simulated data
get.simulation.results <- function(dir) {
  filenames <- file.path(dir, list.files(dir))
  do.call('rbind', lapply(filenames, read.csv))
}

# df.skewt <- get.simulation.results('data/derived/simulation/constant_gauss//')
# raw <- df.skewt
# df.skewt <- df.skewt[, 14:19] / df.skewt[, 8:13]
# colnames(df.skewt) <- factors.all

# Threshold correlations for simulated data ----
do.th.sim <- function(file) {
  df <- get.simulation.results(file)
  stdres <- df[, 14:19] / df[, 8:13]
  colnames(stdres) <- factors.all

  th_corr(stdres, 1)
}

# df.skewt <- get.simulation.results('data/derived/simulation/constant_gauss//')
# raw <- df.skewt
# df.skewt <- df.skewt[, 14:19] / df.skewt[, 8:13]
# colnames(df.skewt) <- factors.all

th.corr.gauss <- do.th.sim('data/derived/simulation/constant_gauss/')
th.corr.ght   <- do.th.sim('data/derived/simulation/constant_ght/')
th.corr.ghskt <- do.th.sim('data/derived/simulation/constant_ghskt/')

th.corr.empirical <- th_corr(df.empirical, 0)

# Draw Threshold Correlations ----

allThCorrList <- list()

for (values in factors.value) {
  out.df <- data.frame(
    # The empirical imports coef, lb, ub, quantiles and ordering of factors
    th.corr.empirical[[values]],

    # All other distributions only imports coefficients
    gauss = th.corr.gauss[[values]]$simcoef,
    ght = th.corr.ght[[values]]$simcoef,
    ghskt = th.corr.ghskt[[values]]$simcoef
  )
  allThCorrList[[values]] <- out.df
}

plotdf <- bind_rows(allThCorrList$HML, allThCorrList$RMW, allThCorrList$CMA)
plots <- correlations.plot.threshold.sim(plotdf)+
  ggtitle("Threshold correlations of weekly data sets (95% confidence bounds)")

# Arrange and save plots
g <- grid.arrange(
  plots +
    coord_cartesian(ylim = c(-0.2, 0.5))
)


# OLD LIKE GUSTAF, STILL YOUNG AND FRESH ----

# df.skewt <- read.table('data/derived/simulation/ghskt/1b813702-897c-11e6-bb9c-9dfbf344bd3b.csv', header = TRUE, sep = ',')
# df.skewt <- df.skewt[,14:19]
# colnames(df.skewt) <- factors.all

# Get the correlation series for each distribution ----


#normThCorrList <- th_corr(df.norm)
#tThCorrList <- th_corr(df.t)
skewtThCorrList <- th_corr(df.skewt, 1)

# load('data/derived/model_copula_constant_ghskt.RData')
# model <- model.copula.constant.ghskt
#
# dist <- ghyp::student.t(
#   nu = model$params$dist.params$df,
#   gamma = model$params$dist.params$skew,
#   mu = rep(0, 6),
#   sigma = cor(model$Correlation)
# )
# df.skewt <- ghyp::rghyp(100000, dist)
# colnames(df.skewt) <- factors.all
# skewtThCorrList <- th_corr(df.skewt, 1)
# cor(df.skewt)

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
