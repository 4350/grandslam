#' Configure constants and required functions and source this file to
#' run an appropriate bootstrap.

library(devtools)
load_all('wimbledon')

library(tictoc)
library(parallel)
library(rugarch)
library(ghyp)

load('data/derived/weekly-estim.RData')
kT <- dim(df.estim)[1]

kBSName <- sprintf(
  "RNG %05d - BLL %04d - REP %04d",
  kRandomSeed,
  kBSBlockLength,
  kBSRepetitions
)

# Functions --------------------------------------------------------------

#' Estimate GARCH models on bootstrap sample
estimate.garch <- function(b.sample) {
  garch <- lapply(names(kGARCHModels), function(name) {
    model <- kGARCHModels[[name]]
    data <- as.data.frame(b.sample[, name])
    ugarchfit(model, data, solver = 'hybrid')
  })
  names(garch) <- names(kGARCHModels)
  
  garch
}

#' Uniform residuals from GARCH models (assume GHSKT for marginal)
extract.uniforms <- function(b.garch) {
  sapply(b.garch, function(model) {
    stdres <- model@fit$residuals / model@fit$sigma
    
    rugarch:::psghst(
      stdres,
      shape = model@fit$coef[['shape']],
      skew = model@fit$coef[['skew']]
    )
  })
}

#' Bootstrap output for GARCH models (parameters and diags)
build.bootstrap.output.garch <- function(b.garch) {
  unlist(lapply(names(b.garch), function(name) {
    model <- b.garch[[name]]
    out <- c(
      model@fit$coef,
      model@fit$LLH,
      model@fit$convergence
    )
    names(out) <- paste0(
      'garch.', name, '.',
      c(names(model@fit$coef), 'll', 'convergence')
    )
    
    out
  }))
}

build.bootstrap.output <- function(b, b.garch, b.copula) {
  # Build a vector of GARCH output, where each element is named
  # "garch.FACTOR.ELEMENT_NAME", e.g. "garch.CMA.beta1"
  #
  # In addition to parameter estimates, we save log likelihood and convergence
  # indicator
  b.out.garch <- build.bootstrap.output.garch(b.garch)
  
  # Build a similar vector from copula estimates
  b.out.copula <- build.bootstrap.output.copula(b.copula)
  
  # Combine into one big vector and prepend the iteration number
  out <- c(b, b.out.garch, b.out.copula)
  names(out)[1] <- 'b'
  
  out
}

# Initialization ---------------------------------------------------------

set.seed(kRandomSeed)

# Load or restore index from file for this bootstrap configuration
path <- file.path(kBSOutputPath, 'index', paste0(kBSName, '.RData'))
if (!file.exists(path)) {
  bs.index <- bs.stationary.index(kBSRepetitions, kT, kBSBlockLength)
  save(bs.index, file = path)
} else {
  load(path)
}
rm(path)

# Bootstrap Loop ----

# Create output folder; each run is saved to a separate file
# It's too much work to save into one big file -- we can do that later
run.directory <- file.path(kBSOutputPath, 'runs', kBSName)
dir.create(run.directory, showWarnings = FALSE)

for (b in seq(kBSIteration, kBSRepetitions)) {
  cat(noquote(strrep('-', 79)))
  cat(noquote(sprintf("%s - ITERATION %04d - (%s)", kBSName, b, Sys.time())))
  
  # Bootstrapped sample
  b.sample <- df.estim[bs.index[, b], ]
  
  # Estimate GARCH models for this sample.
  # Takes roughly 15 seconds on Victor's MacBook -- 4 seconds if parallelized,
  # which isn't worth the complexity since the copula estimation takes
  # ~30 minutes
  cat(noquote("GARCH ESTIMATION"))
  tic()
    b.garch <- estimate.garch(b.sample)
  toc()
  
  # Extract uniform residuals
  # Takes roughly 20 seconds on Victor's MacBook
  cat(noquote("EXTRACT UNIFORM RESIDUALS"))
  tic()
    b.u <- extract.uniforms(b.garch)
  toc()
  
  cat(noquote("COPULA ESTIMATION"))
  tic()
    b.copula.param <- estimate.copula(b.u)
  toc()
  
  # Save this iteration to a data.frame
  b.out <- build.bootstrap.output(b, b.garch, b.copula.param)
  write.csv(
    t(b.out),
    file = file.path(run.directory, paste0(sprintf("%04d", b), '.csv')),
    row.names = F
  )
}