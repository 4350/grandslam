#' Use simulated distributions of returns to get VaR estimates
#' 

# Load Data ----

library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)

rm(list = ls())

# Load distribution 
kModelName <- 'dynamic_ghskt'
load(sprintf('data/derived/distribution_%s.RData', kModelName))
factors <- c("Mkt.RF", "HML", "SMB", "Mom", "RMW", "CMA")
colnames(distribution) <- factors

# Load empirical data
load('data/derived/weekly-full.RData')
df.date <- df %>% select(Date) %>% as.vector()
df <- df %>% select(-Date)

# Run the hist_sim first to get this
load('data/derived/var_hshw_01.RData')
load('data/derived/var_hshw_05.RData')

# Do it ----
var.distribution <- function(series, limit) {
  endT <- dim(series)[3]
  
  data.frame(t(sapply(seq(endT), function(t) {
    apply(series[,, t], 2, function(x) stats::quantile(x, probs = limit))
  })))
}

H <- dim(distribution)[3]
T <- 2479

df.var01.model <- var.distribution(distribution[,, 1:H], limit = c(0.01))
df.var05.model <- var.distribution(distribution[,, 1:H], limit = c(0.05))

save(df.var01.model, file = sprintf('data/derived/var_%s_01.RData', kModelName))
save(df.var05.model, file = sprintf('data/derived/var_%s_05.RData', kModelName))

# Plotting time ----

gather.var <- function(var) {
  gathered <- gather(var, 'factor', 'VaR')
  
  gathered$h <- rep(seq(nrow(var)), ncol(var))
  gathered$order <- factor(gathered$factor, colnames(var))
  gathered
}
  
# Prettify empirical data
df.emp <- gather(df[(T + 1):(T + H), ], 'factor','value', 1:6)
df.emp$h <- rep(seq(1, H), ncol(df))
df.emp$order <- factor(df.emp$factor, factors)

ggplot(df.emp, aes(x = h, y = value, group = factor)) +
  geom_line(aes(color = 'Realized Return')) +
  geom_line(aes(x = h, y = VaR, color = '1% GARCH-Copula VaR'), data = gather.var(df.var01.model)) +
  geom_line(aes(x = h, y = VaR, color = '5% GARCH-Copula VaR'), data = gather.var(df.var05.model)) +
  geom_line(aes(x = h, y = VaR, color = '1% HS-HW VaR'), data = df.var01) +
  geom_line(aes(x = h, y = VaR, color = '5% HS-HW VaR'), data = df.var05) +
  theme_Publication() +
  ylab('') +
  xlab('Horizon') +
  scale_y_continuous(labels = percent) +
  facet_grid(. ~ order)

# Difference ----
df.diff01 <- data.frame(df.var01)
df.diff05 <- data.frame(df.var05)

# If difference is POSITIVE, the historical VAR is HIGHER (= GREATER LOSS)
df.diff01[, 'VaR'] <- abs(df.var01[, 'VaR']) - abs(gather.var(df.var01.model)[, 'VaR'])
df.diff05[, 'VaR'] <- abs(df.var05[, 'VaR']) - abs(gather.var(df.var05.model)[, 'VaR'])

ggplot() +
  geom_line(aes(x = h, y = VaR, color = '1% Difference'), data = df.diff01) +
  geom_line(aes(x = h, y = VaR, color = '5% Difference'), data = df.diff05) +
  theme_Publication() +
  ylab('') +
  xlab('Horizon') +
  scale_y_continuous(labels = percent) +
  facet_grid(. ~ order)


# Expected shortfall ----

var.expected.shortfall <- function(dist, vars) {
  es <- sapply(seq(ncol(vars)), function(i) {
    sapply(seq(nrow(vars)), function(t) {
      dist.it <- dist[, i, t]
      
      # Express losses as positive
      -1 * mean(dist.it[dist.it < vars[t, i]])
    })
  })
  
  colnames(es) <- colnames(vars)
  es
}

# Express losses as positive
df.es01.model <- var.expected.shortfall(distribution, df.var01.model)
df.es05.model <- var.expected.shortfall(distribution, df.var05.model)
