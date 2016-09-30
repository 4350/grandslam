# Output relevant statistics based on the chosen few GARCH models
# 
# Prerequisites:
# 
# * Need to run `do_garch_fit` first to estimate all GARCH models

# Libraries ----
library(devtools)
load_all('wimbledon')

# Reset Worksspace ----
rm(list = ls())
load('data/derived/weekly-estim.RData')
load('data/derived/garch_fits.RData')

# Calculate BIC ----
BIC <- sapply(fits, function(fit) {
  sapply(fit, function(spec) {
    BIC = tryCatch(infocriteria(spec)['Bayes',], error = function(err) NA)
  })
})

rankBIC <- apply(as.data.frame(BIC), 2, min_rank)
write.table(rankBIC, file = 'output/garch_diagnostics/bic_table.csv', sep = ',')

# Pick models with highest BIC ----
# Consolidate these fits in a list to extract all plots

# NOTE: Hard-coded! May want to verify this if changing sample
bestfits <- list(
  Mkt.RF = fits$Mkt.RF$ARMA00,
  HML = fits$HML$ARMA11,
  SMB = fits$SMB$ARMA11,
  Mom = fits$Mom$ARMA10,
  RMW = fits$RMW$ARMA10,  
  CMA = fits$CMA$ARMA10
)
model.GARCH <- bestfits
save(model.GARCH, file = 'data/derived/model_GARCH.RData')


# Get and save standardized residuals ----
df.stdres <- data.frame(
  Date = df.estim$Date,
  as.data.frame(
    lapply(bestfits,
           function(factor) factor@fit$residuals/factor@fit$sigma)
  )
)
save(df.stdres, file = 'data/derived/garch_stdres.RData')

# Get and save uniform residuals ----
df.u <- data.frame(
  Date = df.estim$Date,
  as.data.frame(
    sapply(
      bestfits,
      function(factor) {
        rugarch:::psghst(
          (factor@fit$residuals/factor@fit$sigma),
          shape = factor@fit$coef[['shape']],
          skew = factor@fit$coef[['skew']]
        )
      }
    )
  )
)
save(df.u, file = 'data/derived/garch_unires_model.RData')

# Get and save uniform residuals, empirical CDF ----
df.u.e <- data.frame(
  Date = df.estim$Date,
  
  sapply(
    bestfits,
    function(factor) {
      f = ecdf(factor@fit$residuals)
      f(factor@fit$residuals)
    }
  )
)
save(df.u.e, file = 'data/derived/garch_unires_empirical.RData')

# Output various other things ----

# GARCH parameters
varnames <- list('Mkt.RF', 'HML', 'SMB', 'Mom', 'RMW', 'CMA')
lapply(
  varnames,
  function(vars) {
    write.table(
      bestfits[[vars]]@fit$robust.matcoef,
      file = paste('output/garch_diagnostics/garch_params_', vars, '.csv', sep = ''),
      sep = ','
    )
  }
)

newsimp <- lapply(bestfits, function(fit) {
  newsimpactlist <- rugarch::newsimpact(fit)
  data.frame(x = newsimpactlist$zx, y = newsimpactlist$zy)
})

# Get list of empirical density data for all fits
empdens <- lapply(bestfits, function(fit) garch.empirical.density(fit))

# Make and save the nice diagnostic plots ----
library(ggplot2)
library(ggfortify)
library(gridExtra)
source('func/garch_diagplots.R')
lapply(
  list('Mkt.RF','HML','SMB','Mom','CMA','RMW'),
  function(varlist) garch.diagplots(df.stdres, varlist, newsimp, empdens)
)