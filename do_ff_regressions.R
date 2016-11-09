# Using four factors to explain the fifth --------------------------------


# Data loading and setup --------------------------------------------------

library(sandwich)
library(lmtest)
library(stargazer)
library(dplyr)
rm(list = ls())
load('data/derived/weekly-estim.RData')

# Simplify returns
df <- exp(df.estim[,-1]) - 1
rm(df.estim)


# FF regressions ----------------------------------------------------------

#' Returns linear regression result with robust standard errors
#' @param y String, dependent variable in df
#' @param x String vector, independent variables in df
#' @param df Data frame with x, y
#' 
#' @return coeftest class with parameter estimates and 
#' robust errors
ff_reg <- function(y, x, df) {
  # Get the regression formula
  formula = as.formula(paste(y, " ~ ", 
                             paste(x, collapse = '+')))
  # Run regression
  reg_fit <- lm(formula, df)
  coeftest(reg_fit, vcov = sandwich)
}

ff_reg_rsq <- function(y, x, df) {
  # Get the regression formula
  formula = as.formula(paste(y, " ~ ", 
                             paste(x, collapse = '+')))
  # Run regression
  reg_fit <- lm(formula, df)
  summary(reg_fit)$r.squared
}

#' Returns latex/ascii table of regression results
#' on a data frame where each column takes turns
#' as the dependent variable, the rest being indep.
#' @param df Data frame with x, y
#' 
#' @return latex and ascii tables from stargazer
#' with robust standard errors
ff_reg_exclude_one <- function(df) {
  factors <- colnames(df)
  
  # Regressions
  out_list = lapply(seq(length(factors)), function(f) {
    y = factors[f]
    x = factors[-f]
    ff_reg(y, x, df)
  })
  
  # Rsquareds
  out_rsq = lapply(seq(length(factors)), function(f) {
    y = factors[f]
    x = factors[-f]
    ff_reg_rsq(y, x, df)
  })
  
  names(out_list) = factors
  names(out_rsq) = factors
  stargazer(out_list, order = c('Constant', factors), column.labels = factors, digits = 4)
  stargazer(out_list, order = c('Constant', factors), type = 'text', column.labels = factors, digits = 4)
  out_rsq
}


# Do ----------------------------------------------------------------------

df %>% select(Mkt.RF, SMB, HML, CMA, RMW) %>% ff_reg_exclude_one()

df %>% select(Mkt.RF, SMB, HML, CMA, RMW, Mom) %>% ff_reg_exclude_one()


# Get rsquared ------------------------------------------------------------


