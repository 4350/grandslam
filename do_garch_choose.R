# Output relevant statistics based on the chosen few GARCH models
# 
# Prerequisites:
# 
# * Need to run `do_garch_fitting.R` first to estimate all GARCH models

# Libraries ----
library(devtools)
library(dplyr)
library(rugarch)
library(WeightedPortTest)
library(tidyr)
library(ggplot2)
library(scales)
library(ggfortify)
library(gridExtra)
library(devtools)
library(extrafont)

load_all('wimbledon')

# Reset Worksspace ----
rm(list = ls())
load('data/derived/weekly-estim.RData')


# BIC bestfits and save out stuff -----------------------------------------



do_garch_BIC_bestfits <- function(df.estim, submodel, dist) {
  # Load data
  load(sprintf('data/derived/garch/garch_fits_%s_%s.RData', submodel, dist))
  
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
  
  model.GARCH <- bestfits
  rm(bestfits)
  
  # Save GARCH parameters of best fits to csv files
  lapply(
    colnames(df.estim[,-1]),
    function(vars, SAVEPATH, submodel, dist) {
      write.table(
        model.GARCH[[vars]]@fit$robust.matcoef,
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
  
  # Get and save standardized residuals ----
  df.stdres <- data.frame(
    Date = df.estim$Date,
    as.data.frame(
      lapply(model.GARCH,
             function(factor) factor@fit$residuals/factor@fit$sigma)
    )
  )
  save(df.stdres, file = sprintf('data/derived/garch/garch_stdres_%s_%s.Rdata', submodel, dist))
  
  # Do the QQ plots
  
  g_qq <- do_qq_plot(do_qq_data(model.GARCH))
  ggsave(file = sprintf('output/garch_diagnostics/%s/%s/qqplot.png', submodel, dist),
         g_qq, width = 14.0, height = 10, units = 'cm', limitsize = F
         ) 
  
  
  # Save end table of test results etc
  test.table <- sapply(model.GARCH, do_garch_end_table)
  write.table(test.table, file = sprintf('output/garch_diagnostics/%s/%s/end_table_%s_%s.csv',
                                         submodel, dist, submodel, dist), sep = ',')
  
  # Save model.GARCH to file
  save(model.GARCH, file = sprintf('data/derived/garch/model_GARCH_%s_%s.Rdata', submodel, dist))
}

do_garch_BIC_bestfits(df.estim, 'GJRGARCH', 'ghst')
do_garch_BIC_bestfits(df.estim, 'GJRGARCH', 'std')
do_garch_BIC_bestfits(df.estim, 'GJRGARCH', 'norm')

do_garch_BIC_bestfits(df.estim, 'GARCH', 'ghst')
do_garch_BIC_bestfits(df.estim, 'GARCH', 'std')
do_garch_BIC_bestfits(df.estim, 'GARCH', 'norm')

# Below only for the chosen -----------------------------------------------
load('data/derived/garch/model_GARCH_GARCH_ghst.Rdata')
model.GARCH.GARCH <- model.GARCH
rm(model.GARCH)
load('data/derived/garch/model_GARCH_GJRGARCH_ghst.Rdata')
model.GARCH.GJRGARCH <- model.GARCH
rm(model.GARCH)

bestfits <- 
  list(
    Mkt.RF = model.GARCH.GJRGARCH$Mkt.RF,
    HML = model.GARCH.GARCH$HML,
    SMB = model.GARCH.GARCH$SMB,
    Mom = model.GARCH.GARCH$Mom,
    RMW = model.GARCH.GARCH$RMW,
    CMA = model.GARCH.GJRGARCH$CMA
  )

model.GARCH <- bestfits
rm(bestfits)
save(model.GARCH, file = 'data/derived/garch/model_GARCH_chosen.RData')
rm(model.GARCH.GJRGARCH, model.GARCH.GARCH)

# Get and save residuals ----
df.res <- data.frame(
  Date = df.estim$Date,
  as.data.frame(
    lapply(model.GARCH,
           function(factor) factor@fit$residuals)
  )
)
save(df.res, file = 'data/derived/garch/model_GARCH_chosen_res.RData')

# Get and save uniform residuals ----
df.u <- data.frame(
  Date = df.estim$Date,
  as.data.frame(
    sapply(
      model.GARCH,
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
save(df.u, file = 'data/derived/garch/model_GARCH_chosen_u.RData')

# Get and save standardized residuals ----
df.stdres <- data.frame(
  Date = df.estim$Date,
  as.data.frame(
    lapply(model.GARCH,
           function(factor) factor@fit$residuals/factor@fit$sigma)
  )
)
save(df.stdres, file = 'data/derived/garch/model_GARCH_chosen_stdres.RData')

# Below only for the chosen  OOS -----------------------------------------------
rm(list = ls())
load('data/derived/weekly-estim.RData')
load('data/derived/garch/oos_garch_fits_GARCH_ghst.RData')
garch_fits_GARCH <- fits
rm(fits)
load('data/derived/garch/oos_garch_fits_GJRGARCH_ghst.RData')
garch_fits_GJRGARCH <- fits
rm(fits)

oos_bestfits <- 
  list(
    Mkt.RF = garch_fits_GJRGARCH$Mkt.RF$ARMA00,
    HML = garch_fits_GARCH$HML$ARMA11,
    SMB = garch_fits_GARCH$SMB$ARMA11,
    Mom = garch_fits_GARCH$Mom$ARMA10,
    RMW = garch_fits_GARCH$RMW$ARMA11,
    CMA = garch_fits_GJRGARCH$CMA$ARMA11
  )

model.GARCH <- oos_bestfits
rm(oos_bestfits)
save(model.GARCH, file = 'data/derived/garch/oos_model_GARCH_chosen.RData')
rm(garch_fits_GARCH, garch_fits_GJRGARCH)

nOOS <- length(model.GARCH$Mkt.RF@fit$residuals)

# Get and save residuals ----
df.res <- data.frame(
  Date = tail(df.estim$Date, nOOS),
  as.data.frame(
    lapply(model.GARCH,
           function(factor) factor@fit$residuals)
  )
)
save(df.res, file = 'data/derived/garch/oos_model_GARCH_chosen_res.RData')

# Get and save uniform residuals ----
df.u <- data.frame(
  Date = tail(df.estim$Date, nOOS),
  as.data.frame(
    sapply(
      model.GARCH,
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
save(df.u, file = 'data/derived/garch/oos_model_GARCH_chosen_u.RData')

# Get and save standardized residuals ----
df.stdres <- data.frame(
  Date = tail(df.estim$Date, nOOS),
  as.data.frame(
    lapply(model.GARCH,
           function(factor) factor@fit$residuals/factor@fit$sigma)
  )
)
save(df.stdres, file = 'data/derived/garch/oos_model_GARCH_chosen_stdres.RData')