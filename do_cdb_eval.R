library(ggplot2)
library(tidyr)

rm(list = ls())

get_cdb <- function(name) {
  load(sprintf('data/derived/cdb/%s.RData', name))
  cdb_results
}

cdb_results <- list(
  CMA = get_cdb('CMA_full_dynamic_ghst_10000'),
  HML = get_cdb('HML_full_dynamic_ghst_10000')
)

cdb <- data.frame(
  t = 1:length(cdb_results[[1]]$cdb),
  lapply(cdb_results, function(c) c$cdb)
)

ggplot(
  gather(cdb, Strategy, CDB, -1),
  aes(x = t, y = CDB, color = Strategy)
) + geom_line()