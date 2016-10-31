#' ARCHIVED BECAUSE THIS GARCH MODEL IS NOT CURRENT NOR IS THE
#' COPULA INTERFACE CURRENT

#' Filter models
#' 
#' Give me estimated models and I will give you their filtered output against
#' data. I don't care what data you fitted your models against.
#' 
#' I expect you to give me the following parameters:
#' 
#' kSeries: A zoo object of the data you're filtering me against
#' kGARCHModels: List of ugarchfit objects
#' kCopulaModel: A parametrization of a copula model
#' kOutputFile: Where do I save my results?

# Demo ----
library(zoo)
rm(list = ls())

# Get series
load('data/derived/weekly-full.RData')
kSeries <- zoo(df[, -1], order.by = df$Date)
rm(df)

# Load some nice models
load('data/derived/model_GARCH.RData')
load('data/derived/model_copula_dynamic_ghskt.RData')

kGARCHModels <- model.GARCH
kCopulaModel <- model.copula.dynamic.ghskt
kOutputPath = 'data/derived/filtered_dynamic_ghskt.RData'

# Setup ------------------------------------------------------------------

library(devtools)
library(rugarch)
library(zoo)
load_all('wimbledon')

# Build GARCH specifications from garchfit objects to ensure that the fitted
# data does not pollute it 
kGARCHSpecs <- garch.fit2spec(kGARCHModels)

# Filter GARCH models ----------------------------------------------------

filtered <- lapply(seq(kGARCHSpecs), function(i) {
  ugarchfilter(kGARCHSpecs[[i]], kSeries[, i])
})
names(filtered) <- colnames(kSeries)

kResid <- zoo(sapply(filtered, rugarch::residuals), order.by = index(kSeries))
kSigma <- zoo(sapply(filtered, rugarch::sigma), order.by = index(kSeries))
rm(filtered)

# Extract uniform residuals
# Takes relatively long time
stdresid <- kResid / kSigma
u <- sapply(seq(kGARCHSpecs), function(i) {
  pars <- kGARCHSpecs[[i]]@model$fixed.pars
  
  rugarch:::psghst(
    stdresid[, i],
    shape = pars['shape'],
    skew = pars['skew']
  )
})
colnames(u) <- colnames(kSeries)
u <- zoo(u, order.by = index(kSeries))

rm(stdresid)

# Filter copula ----------------------------------------------------------

copula <- dc.run.model(
  u,
  dist.params = kCopulaModel$params$dist.params,
  alpha = kCopulaModel$params$alpha,
  beta = kCopulaModel$params$beta,
  
  # We *DON'T* want to re-estimate Omega but rather we take it as given
  # from before (it *is* a model parameter, however, estimation is by
  # method of moments eqn)
  Omega = kCopulaModel$Omega
)

# Save results -----------------------------------------------------------

filtered.series <- list(
  filtered = list(
    copula = copula,
    series = kSeries,
    resid = kResid,
    sigma = kSigma
  ),
  copula = kCopulaModel,
  garch = kGARCHSpecs
)

save(filtered.series, file = kOutputPath)