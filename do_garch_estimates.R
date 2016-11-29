rm(list = ls())

library(tidyr)
library(stargazer)

factors <- list(
  "Mkt.RF",
  "SMB",
  "HML",
  "CMA",
  "RMW",
  "Mom"
)
estimates <- c(
  'mu',
  'ar1',
  'ma1',
  'alpha1',
  'beta1',
  'gamma1',
  'shape',
  'skew',
  'omega'
)

load('data/derived/garch/model_GARCH_chosen.RData')

df2 <- lapply(factors, function(factor) {
  coef <- model.GARCH[[factor]]@fit$robust.matcoef
  colnames(coef) <- c('coef', 'se', 't', 'p')
  coef['mu', c('coef', 'se')] <- 100 * coef['mu', c('coef', 'se')]
  
  coef[, 'p'] <- 100 * coef[, 'p']
  data.frame(coef)
})
names(df2) <- factors

selfie <- function(dfs, column, digits = 2) {
  table <- bind_rows(lapply(dfs, function(df) {
    points <- df[[column]]
    names(points) <- rownames(df)
    data.frame(as.list(points))
  }))
  rownames(table) <- names(dfs)
  stargazer(
    t(table[estimates]),
    type = 'text',
    summary = F,
    digits = digits,
    digits.extra = 0
  )
}

selfie(df2, 'coef', digits = 2)
selfie(df2, 'se', digits = 2)
selfie(df2, 'p', digits = 0)
# lapply(df2, selfie('coef'))