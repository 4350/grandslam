#' Simulate from constant copulas

# Setup Workspace --------------------------------------------------------

library(devtools)
library(parallel)
library(tictoc)
library(uuid)
load_all('wimbledon')

rm(list = ls())

copula.do.sim <- function(copula, garch, T, output.dir, block.size = 8) {
  cluster <- prepare.cluster()

  simulate <- function(i, T, garch, copula) {
    # Simulate...
    result <- sim.c(T, garch, copula)

    # Pretty results
    data.frame(
      series = result$series,
      sigma = result$sigma,
      resid = result$resid
    )
  }

  n <- 0
  while (TRUE) {
    tic()
    results <- parLapplyLB(
      cluster,
      seq(block.size),
      simulate,

      T = T,
      garch = garch,
      copula = copula
    )
    toc()

    lapply(results, function(result) {
      file <- file.path(output.dir, paste0(UUIDgenerate(), '.csv'))
      write.csv(result, file = file)
    })

    n <- n + length(results)
    cat("TOTAL: ", n, "\n")
  }

  stopCluster()
}

# Constant Asymmetric Student's t ---------------------------------------

load('data/derived/model_GARCH.RData')
load('data/derived/model_copula_constant_gauss.RData')

kCopulaModel <- model.copula.constant.gauss
kGARCHModels <- model.GARCH
kOutputPath <- 'data/derived/simulation/constant_gauss'

dir.create(kOutputPath, recursive = T, showWarnings = F)

copula.do.sim(kCopulaModel, kGARCHModels, 1000, kOutputPath)

# Debugging --------------------------------------------------------------

copula.sim <- copula.sim <- sim.c.copula(10000, kCopulaModel)
colnames(copula.sim) <- c("Mkt.RF", "HML", "SMB", "Mom", "RMW", "CMA")

th.corr.sim <- correlations.threshold(copula.sim, "Mkt.RF", "HML")
series.sim <- sim.c.GARCH(copula.sim, kGARCHModels)

colnames(series.sim$resid) <- c("Mkt.RF", "HML", "SMB", "Mom", "RMW", "CMA")
colnames(series.sim$series) <- c("Mkt.RF", "HML", "SMB", "Mom", "RMW", "CMA")
colnames(series.sim$sigma) <- c("Mkt.RF", "HML", "SMB", "Mom", "RMW", "CMA")

stdresid <- series.sim$resid / series.sim$sigma
th.corr.series <- correlations.threshold(series.sim$series, "Mkt.RF", "HML")
th.corr.resid <- correlations.threshold(series.sim$resid, "Mkt.RF", "HML")
th.corr.stdresid <- correlations.threshold(stdresid, "Mkt.RF", "HML")

plot(th.corr.sim['coef', ], type = 'l')
lines(th.corr.resid['coef', ], col = 'red')
lines(th.corr.series['coef', ], col = 'green')
lines(th.corr.stdresid['coef', ], col = 'gray')
