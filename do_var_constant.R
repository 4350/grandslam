# Setup ------------------------------------------------------------------

library(ghyp)
library(rugarch)
library(zoo)
library(tictoc)

rm(list = ls())
MODEL_NAME <- 'constant_ghskt'
N <- 1000

load(sprintf('data/derived/model_copula_%s.RData', MODEL_NAME))
load('data/derived/model_GARCH.RData')
load('data/derived/weekly-full.RData')
model <- get(sprintf('model.copula.%s', gsub('_', '.', MODEL_NAME)))

garch <- garch.fit2spec(model.GARCH)

# Filter series ----
series <- zoo(df[, -1], order.by = df$Date)
filtered <- lapply(seq_along(garch),
                   function(i) ugarchfilter(garch[[i]], series[, i]))
names(filtered) <- colnames(series)

resid <- zoo(sapply(filtered, rugarch::residuals), order.by = index(series))
sigma <- zoo(sapply(filtered, rugarch::sigma), order.by = index(series))

filtered <- list(
  series = series,
  resid = resid,
  sigma = sigma
)

# Simulate standardized residuals ----
# Distribution same for all T with constant
mv_distribution <- dc.mv.dist(model$params$dist.params, diag(1, 6))
uv_distributions <- dc.uv.dists(ncol(model$Omega), model$params$dist.params)

shocks <- ghyp::rghyp(N, mv_distribution)
stdresid <- shocks2stdresid(shocks, uv_distributions, garch)

# Distributions ----

time <- 2479:2766
series <- array(NA, dim = c(N, 6, length(time)))
colnames(series) <- names(garch)

tic("Simulate t+1 paths for all factors")
for (t in seq(time)) {
  for (i in 1:6) {
    path <- garch.path.it(i, time[t], garch, filtered, stdresid)
    series[, i, t] <- t(path@path$seriesSim)
  }
}
toc()

# VAR and ES ----

var.distribution <- function(series, limit) {
  endT <- dim(series)[3]
  
  data.frame(t(sapply(seq(endT), function(t) {
    apply(series[,, t], 2, function(x) stats::quantile(x, probs = limit))
  })))
}

var.expected.shortfall <- function(dist, vars) {
  es <- sapply(seq(ncol(vars)), function(i) {
    sapply(seq(nrow(vars)), function(t) {
      dist.it <- dist[, i, t]
      
      # Express losses as positive
      -1 * mean(dist.it[dist.it < vars[t, i]])
    })
  })
  
  colnames(es) <- colnames(vars)
  data.frame(es)
}

df.var01.model <- var.distribution(series, 0.01)
df.var05.model <- var.distribution(series, 0.05)
df.es01.model <- var.expected.shortfall(series, df.var01.model)
df.es05.model <- var.expected.shortfall(series, df.var05.model)

var_es <- list(
  var01 = df.var01.model,
  var05 = df.var05.model,
  es01 = df.es01.model,
  es05 = df.es05.model
)

save(var_es, file = sprintf('data/derived/var_es_%s.RData', MODEL_NAME))
