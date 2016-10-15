#' Create marginal summary statistics table and plots

# Libraries ----
library(dplyr)
library(fBasics)
library(tidyr)
library(ggplot2)
library(scales)
library(ggfortify)
library(gridExtra)
library(devtools)
library(extrafont)
library(stargazer)

# Reset workspace and load return data ----
rm(list = ls())
load('data/derived/weekly-estim.RData')

# Load functions ----
load_all('wimbledon')

# Descriptive statistics tables ---
sum.table <- summary.table(df.estim, 'Estim')
#stargazer(sum.table, summary = FALSE)

# Marginal plots ----
varlist = list('Mkt.RF','HML','SMB','Mom','RMW','CMA')

lapply(varlist,
       function(varlist) summary.plots(df.estim, 'Estim', varlist))

# Cumulative plots
summary.cumretplots(df.estim, 'Estim')