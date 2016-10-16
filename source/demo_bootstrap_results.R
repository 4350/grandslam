library(ggplot2)
library(grid)
library(tidyr)

rm(list = ls())
path <- 'data/derived/bootstrap/ghskt/runs/RNG 00403 - BLL 0045 - REP 1000'

filenames <- file.path(path, list.files(path))
bootstraps <- do.call('rbind', lapply(filenames, read.csv))

# Copula parameters
gathered.copula <- gather(bootstraps, 'param', 'value', 63:71)

ggplot(gathered.copula, aes(value)) +
  geom_histogram(bins = 20) +
  facet_grid(. ~ param, scales = 'free')