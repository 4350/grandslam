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
library(WeightedPortTest)

# Reset workspace and load return data ----
rm(list = ls())
load('data/derived/weekly-estim.RData')

# Load functions ----
load_all('wimbledon')


# LB, ARCH LM -------------------------------------------------------------

ret_LB_5 <- apply(df.estim[,-1], 2, function(ret) Weighted.Box.test(ret, lag = 5, type = "Ljung-Box", 
                                                                    weighted = TRUE)$p.value)
ret_LB_10 <- apply(df.estim[,-1], 2, function(ret) Weighted.Box.test(ret, lag = 10, type = "Ljung-Box", 
                                                                     weighted = TRUE)$p.value)

sqr_ret_LB_5 <- apply(df.estim[,-1], 2, function(ret_sqr) Weighted.Box.test(ret_sqr, lag = 5, type = "Ljung-Box",
                                                                            weighted = TRUE,
                                                                            sqrd.res = TRUE)$p.value)
sqr_ret_LB_10 <- apply(df.estim[,-1], 2, function(ret_sqr) Weighted.Box.test(ret_sqr, lag = 10, type = "Ljung-Box",
                                                                             weighted = TRUE,
                                                                             sqrd.res = TRUE)$p.value)

#' Summary statistics table

summary_table <- df.estim %>%
    dplyr::select(-Date) %>%
    basicStats() %>%
    .[c('nobs','Maximum','Minimum','Mean','Median','Stdev','Skewness','Kurtosis'),] %>%
    round(., digits = 4)

summary_table <- rbind(summary_table,
                       'Return LB [5] p-value' = ret_LB_5,
                       'Return LB [10] p-value' = ret_LB_10,
                       'Squared return LB [5] p-value' = sqr_ret_LB_5,
                       'Squared return LB [10] p-value' = sqr_ret_LB_10) %>%
  round(., digits = 4)

write.table(summary_table, file = 'output/MarginalStats/summaryTable.Estim.csv')
#stargazer(sum.table, summary = FALSE)

# Marginal plots ----
varlist = list('Mkt.RF','HML','SMB','Mom','RMW','CMA')

lapply(varlist,
       function(varlist) summary.plots(df.estim, 'Estim', varlist))

# Cumulative plots
summary.cumretplots(df.estim, 'Estim')