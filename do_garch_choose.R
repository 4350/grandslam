# Output relevant statistics based on the chosen few GARCH models
# 
# Prerequisites:
# 
# * Need to run `do_garch_fitting.R` first to estimate all GARCH models

# Libraries ----
library(devtools)
library(dplyr)
library(rugarch)
load_all('wimbledon')

# Reset Worksspace ----
rm(list = ls())
load('data/derived/weekly-estim.RData')


# BIC bestfits and save out stuff -----------------------------------------



do_garch_BIC_bestfits <- function(df.estim, submodel, dist) {
  # Load data
  load(sprintf('data/derived/garch_fits_%s_%s.RData', submodel, dist))
  
  # Create folders for output
  SAVEPATH = 'output/garch_diagnostics'
  save.directory <- file.path(SAVEPATH, submodel, dist)
  dir.create(save.directory, showWarnings = FALSE, recursive = TRUE)
  
  # Calculate BIC
  BIC <- sapply(fits, function(fit) {
    sapply(fit, function(spec) {
      BIC = tryCatch(infocriteria(spec)['Bayes',], error = function(err) NA)
    })
  })
  
  # Rank BIC and write table of BICs
  rankBIC <- apply(as.data.frame(BIC), 2, min_rank)
  write.table(rankBIC, 
              file = sprintf('%s/bic_table_%s_%s.csv',
                                      file.path(SAVEPATH, submodel, dist),
                                      submodel,
                                      dist),
              sep = ','
              )
  
  # Pick models with highest BIC. Consolidate these fits in bestfits list
  best_list_name <- rownames(rankBIC)[apply(rankBIC, 2, function(f) match(1,f))]
  bestfits <- mapply(
    function(f, best_spec_name) {
      f[[best_spec_name]]
    }, 
    fits,
    best_list_name
  )
  
  # Save GARCH parameters of best fits to csv files
  lapply(
    colnames(df.estim[,-1]),
    function(vars, SAVEPATH, submodel, dist) {
      write.table(
        bestfits[[vars]]@fit$robust.matcoef,
        file = sprintf('%s/garch_params_%s_%s_%s.csv',
                       file.path(SAVEPATH, submodel, dist),
                       vars,
                       submodel,
                       dist),
        sep = ','
      )
    },
    SAVEPATH = SAVEPATH,
    submodel = submodel,
    dist = dist
  )
  
  # Name the bestfits model.GARCH
  model.GARCH <- bestfits
  
  # Do the QQ plots
  
  g_qq <- do_qq_plot(do_qq_data(model.GARCH))
  ggsave(file = sprintf('output/garch_diagnostics/%s/%s/qqplot.png', submodel, dist),
         g_qq, width = 14.0, height = 10, units = 'cm', limitsize = F
         ) 
  
  
  # Save model.GARCH to file
  save(model.GARCH, file = sprintf('data/derived/model_GARCH_%s_%s.Rdata', submodel, dist))
}

do_garch_BIC_bestfits(df.estim, 'GJRGARCH', 'ghst')
do_garch_BIC_bestfits(df.estim, 'GJRGARCH', 'std')
do_garch_BIC_bestfits(df.estim, 'GJRGARCH', 'norm')

do_garch_BIC_bestfits(df.estim, 'GARCH', 'ghst')
do_garch_BIC_bestfits(df.estim, 'GARCH', 'std')
do_garch_BIC_bestfits(df.estim, 'GARCH', 'norm')

# Below only for the chosen -----------------------------------------------



# Get and save residuals ----
df.res <- data.frame(
  Date = df.estim$Date,
  as.data.frame(
    lapply(bestfits,
           function(factor) factor@fit$residuals)
  )
)
save(df.res, file = 'data/derived/garch_res.RData')

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

newsimp <- lapply(bestfits, function(fit) {
  newsimpactlist <- rugarch::newsimpact(fit)
  data.frame(x = newsimpactlist$zx, y = newsimpactlist$zy)
})

# Get list of empirical density data for all fits
empdens <- lapply(bestfits, function(fit) garch.empirical.density(fit))



# GARCH diagnostics -------------------------------------------------------

# QQ PLOT
# Sign bias




# Make and save the nice diagnostic plots ----
library(ggplot2)
library(ggfortify)
library(gridExtra)
library(extrafont)
lapply(
  list('Mkt.RF','HML','SMB','Mom','RMW','CMA'),
  function(varlist) garch.diagplots(df.stdres, varlist, newsimp, empdens)
)