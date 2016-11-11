# Main ----
library(dplyr)
library(ggplot2)
library(tidyr)
library(devtools)
library(zoo)

load_all('wimbledon')

rm(list = ls())

load('data/derived/weekly-estim.RData')
MODEL_NAME <- 'full_dynamic_std_10000'

factor_models <- list(
  "5F" = list(
    '5F',
    '5F_EXCL_CMA',
    '5F_EXCL_HML'
  ),
  "6F" = list(
    '6F',
    '6F_EXCL_CMA',
    '6F_EXCL_HML'
  )
)

get_cdb <- function(name) {
  load(sprintf('data/derived/cdb/%s.RData', name))
  cdb_results
}

SMOOTH <- 13

cdb_cdb <- lapply(names(factor_models), function(factor_name) {
  factor_strategy_names <- factor_models[[factor_name]]
  factor_cdb <- lapply(factor_strategy_names, function(name) {
    cdb_name <- sprintf('%s_%s', MODEL_NAME, name)
    data.frame(
      Week = df.estim$Date[(SMOOTH + 1):length(df.estim$Date)],
      CDB = rollmeanr(get_cdb(cdb_name)$cdb, SMOOTH)
    )
  })
  
  names(factor_cdb) <- c('ALL', 'EXCL_CMA', 'EXCL_HML')
  bind_rows(factor_cdb, .id = 'Strategy')
})
names(cdb_cdb) <- c('5 factors', '6 factors')
cdb_cdb <- bind_rows(cdb_cdb, .id = 'Factors')

cdb_cdb$Strategy <- factor(
  cdb_cdb$Strategy,
  c('ALL', 'EXCL_CMA', 'EXCL_HML'),
  c('All', 'Excl. CMA', 'Excl. HML')
)

g <- ggplot(cdb_cdb, aes(x = Week, y = CDB, color = Strategy)) +
  geom_line() +
  facet_grid(Factors ~ .) +
  coord_cartesian(ylim = c(0.70, 1.00)) +
  scale_colour_Publication() +
  theme_Publication()

g

ggsave(
  'output/cdb/cdb_5F_6F.png',
  g,
  width = 21.0,
  height = 21,
  units = 'cm',
  dpi = 144
)

# Means and hypothesis testing ----

cdb <- lapply(factor_models, function(factor) {
  cdb <- lapply(factor, function(strategy) get_cdb(sprintf('%s_%s', MODEL_NAME, strategy))$cdb)
  names(cdb) <- factor
  cdb
})

cdb <- data.frame(cdb)

# H0: X == Y
pairs <- list(
  list(cdb$X5F.5F, cdb$X5F.5F_EXCL_CMA),
  list(cdb$X5F.5F, cdb$X5F.5F_EXCL_HML),
  list(cdb$X5F.5F_EXCL_CMA, cdb$X5F.5F_EXCL_HML),
  
  list(cdb$X6F.6F, cdb$X6F.6F_EXCL_CMA),
  list(cdb$X6F.6F, cdb$X6F.6F_EXCL_HML),
  list(cdb$X6F.6F_EXCL_CMA, cdb$X6F.6F_EXCL_HML),
  
  list(cdb$X6F.6F, cdb$X5F.5F)
)

tested <- bind_rows(lapply(pairs, function(pair) {
  test <- t.test(100 * pair[[1]], 100 * pair[[2]], paired = T)
  list(
    statistic = test$statistic,
    estimate = test$estimate,
    se = test$estimate / test$statistic,
    p = test$p.value
  )
}))

stargazer::stargazer(100 * cdb, type = 'text')
stargazer::stargazer(data.frame(tested), type = 'text', summary = F, digits = 3, digits.extra = 0)
