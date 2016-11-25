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

do_garch_BIC_bestfits <- function(df.estim, model, dist, nIS = 0) {
  prefix = ""
  if(nIS != 0) {
    prefix = "oos_"
    df.estim <- head(df.estim, nIS)
  }
  # Load data
  load(sprintf('data/derived/garch/%sgarch_fits_%s_%s.RData', prefix, model, dist))
  
  # Create folders for output
  SAVEPATH = 'output/garch_diagnostics'
  save.directory <- file.path(SAVEPATH, model, dist)
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
              file = sprintf('%s/%sbic_table_%s_%s.csv',
                                      file.path(SAVEPATH, model, dist),
                             prefix,         
                             model,
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
    function(vars, SAVEPATH, prefix, model, dist) {
      write.table(
        model.GARCH[[vars]]@fit$robust.matcoef,
        file = sprintf('%s/%sgarch_params_%s_%s_%s.csv',
                       file.path(SAVEPATH, model, dist),
                       prefix,
                       vars,
                       model,
                       dist),
        sep = ','
      )
    },
    SAVEPATH = SAVEPATH,
    prefix = prefix,
    model = model,
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
  save(df.stdres, file = sprintf('data/derived/garch/%sgarch_stdres_%s_%s.Rdata', prefix, model, dist))
  
  # Do the QQ plots
  
  g_qq <- do_qq_plot(do_qq_data(model.GARCH))
  ggsave(file = sprintf('output/garch_diagnostics/%s/%s/%sqqplot_%s_%s.png', model, dist, prefix, model, dist),
         g_qq, width = 9, height = 5.8, units = 'cm', limitsize = F
         ) 
  
  
  # Save end table of test results etc
  test.table <- sapply(model.GARCH, do_garch_end_table)
  write.table(test.table, file = sprintf('output/garch_diagnostics/%s/%s/%send_table_%s_%s.csv',
                                         model, dist, prefix, model, dist), sep = ',')
  
  # Save model.GARCH to file
  save(model.GARCH, file = sprintf('data/derived/garch/%smodel_GARCH_%s_%s.Rdata', prefix, model, dist))
}

# Total
do_garch_BIC_bestfits(df.estim, 'gjrGARCH', 'ghst')
do_garch_BIC_bestfits(df.estim, 'gjrGARCH', 'std')
do_garch_BIC_bestfits(df.estim, 'gjrGARCH', 'norm')

do_garch_BIC_bestfits(df.estim, 'sGARCH', 'ghst')
do_garch_BIC_bestfits(df.estim, 'sGARCH', 'std')
do_garch_BIC_bestfits(df.estim, 'sGARCH', 'norm')

# OOS
do_garch_BIC_bestfits(df.estim, 'gjrGARCH', 'ghst', 1852)
do_garch_BIC_bestfits(df.estim, 'sGARCH', 'ghst', 1852)


# Below only for the chosen -----------------------------------------------
rm(list = ls())
load('data/derived/weekly-estim.RData')
load('data/derived/garch/model_GARCH_sGARCH_ghst.Rdata')
model.GARCH.sGARCH <- model.GARCH
rm(model.GARCH)
load('data/derived/garch/model_GARCH_gjrGARCH_ghst.Rdata')
model.GARCH.gjrGARCH <- model.GARCH
rm(model.GARCH)

bestfits <- 
  list(
    Mkt.RF = model.GARCH.gjrGARCH$Mkt.RF,
    HML = model.GARCH.sGARCH$HML,
    SMB = model.GARCH.sGARCH$SMB,
    Mom = model.GARCH.sGARCH$Mom,
    RMW = model.GARCH.sGARCH$RMW,
    CMA = model.GARCH.sGARCH$CMA
  )

model.GARCH <- bestfits
rm(bestfits)
save(model.GARCH, file = 'data/derived/garch/model_GARCH_chosen.RData')
rm(model.GARCH.sGARCH, model.GARCH.gjrGARCH)


# QQ plots for best fits --------------------------------------------------

# Do the QQ plots

g_qq_bestfits <- do_qq_plot(do_qq_data(model.GARCH))
ggsave(file = 'output/garch_diagnostics/qqplot_bestfits.png',
       g_qq_bestfits, width = 14.0, height = 8, units = 'cm', limitsize = F
) 


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

# Filter ----
specs <- garch.fit2spec(model.GARCH)

df <- data.frame(df.estim)[, -1]

filtered <- lapply(seq_along(specs), function(i) {
  filtered <- ugarchfilter(specs[[i]], c(df[, i]))
  
  stopifnot(head(filtered@filter$z) == head(model.GARCH[[i]]@fit$z))
  filtered
})
names(filtered) <- names(model.GARCH)

save(filtered, file = 'data/derived/garch/model_GARCH_chosen_filtered.RData')

# Below only for the chosen (hyperparameters from full sample) OOS -----------------------------------------------
rm(list = ls())
load('data/derived/weekly-estim.RData')
load('data/derived/garch/oos_garch_fits_sGARCH_ghst.RData')
oos.sGARCH.fits <- fits
rm(fits)
load('data/derived/garch/oos_garch_fits_gjrGARCH_ghst.RData')
oos.gjrGARCH.fits <- fits
rm(fits)

# NOTE HARDCODED IMPORTANT TO CHANGE IF SPECS CHANGE
####################################################
oos_bestfits <- 
  list(
    Mkt.RF = oos.gjrGARCH.fits$Mkt.RF$ARMA00,
    HML = oos.sGARCH.fits$HML$ARMA11,
    SMB = oos.sGARCH.fits$SMB$ARMA11,
    Mom = oos.sGARCH.fits$Mom$ARMA10,
    RMW = oos.sGARCH.fits$RMW$ARMA11,
    CMA = oos.sGARCH.fits$CMA$ARMA10
  )

model.GARCH <- oos_bestfits
rm(oos_bestfits)
save(model.GARCH, file = 'data/derived/garch/oos_model_GARCH_chosen.RData')
rm(oos.sGARCH.fits, oos.gjrGARCH.fits)

nIS <- length(model.GARCH$Mkt.RF@fit$residuals)

# Get and save residuals ----
df.res <- data.frame(
  Date = head(df.estim$Date, nIS),
  as.data.frame(
    lapply(model.GARCH,
           function(factor) factor@fit$residuals)
  )
)
save(df.res, file = 'data/derived/garch/oos_model_GARCH_chosen_res.RData')

# Get and save uniform residuals ----
df.u <- data.frame(
  Date = head(df.estim$Date, nIS),
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
  Date = head(df.estim$Date, nIS),
  as.data.frame(
    lapply(model.GARCH,
           function(factor) factor@fit$residuals/factor@fit$sigma)
  )
)
save(df.stdres, file = 'data/derived/garch/oos_model_GARCH_chosen_stdres.RData')

# Get and save filtered data for full sample ----
specs <- garch.fit2spec(model.GARCH)
n.old <- length(model.GARCH$Mkt.RF@fit$z)

df <- data.frame(df.estim)[, -1]

filtered <- lapply(seq_along(specs), function(i) {
  filtered <- ugarchfilter(specs[[i]], c(df[, i]), n.old = n.old)
  
  stopifnot(head(filtered@filter$z) == head(model.GARCH[[i]]@fit$z))
  filtered
})
names(filtered) <- names(model.GARCH)

save(filtered, file = 'data/derived/garch/oos_model_GARCH_chosen_filtered.RData')
