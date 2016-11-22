# Setup -------------------------------------------------------------------

rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(gtable)
library(grid)
library(gridExtra)
library(extrafont)
library(zoo)
library(scales)
library(devtools)

load_all('wimbledon')

load('data/derived/rolling/roll_corr_stdres.RData')
load('data/derived/rolling/roll_corr_returns.RData')
load('data/derived/rolling/roll_corr_simulated.RData')
load('data/derived/garch/model_GARCH_chosen_stdres.RData')

# Set ordering for plots

MODEL_ORDER <- c('stdres','returns','norm','std','ghst')

ORDER_PAGE_1 <- list(
    c('Mkt.RF', 'HML'),
    c('Mkt.RF', 'CMA'),
    c('Mkt.RF', 'RMW'),
    c('Mom', 'HML'),
    c('Mom', 'CMA'),
    c('Mom', 'RMW')
)

ORDER_PAGE_2 <- list(
    c('SMB', 'HML'),
    c('SMB', 'CMA'),
    c('SMB', 'RMW'),
    c('HML', 'CMA'),
    c('HML', 'RMW'),
    c('CMA', 'RMW')
  )

WIDTH = 16
HEIGHT = 20

# Print the pages with plots ----------------------------------------------

# Main section
g <- generate_page(ORDER_PAGE_1, plot_roll_main)
ggsave('output/rollingCorrelations/rolling1.png', g, width = WIDTH, height = HEIGHT, units = 'cm')

g <- generate_page(ORDER_PAGE_2, plot_roll_main)
ggsave('output/rollingCorrelations/rolling2.png', g, width = WIDTH, height = HEIGHT, units = 'cm')

# Appendix section
g <- generate_page(ORDER_PAGE_1, plot_roll_appendix, c('ARMA-GARCH standardized residuals', 'Returns'))
ggsave('output/rollingCorrelations/appendix_rolling1.png', g, width = WIDTH, height = HEIGHT, units = 'cm')

g <- generate_page(ORDER_PAGE_2, plot_roll_appendix, c('ARMA-GARCH standardized residuals', 'Returns'))
ggsave('output/rollingCorrelations/appendix_rolling2.png', g, width = WIDTH, height = HEIGHT, units = 'cm')

# Simulated section
g <- generate_page(ORDER_PAGE_1, plot_roll_simulated, c('Standardized residuals', 'Symmetric t copula standardized residuals'))
ggsave('output/rollingCorrelations/rolling_simulated1.png', g, width = WIDTH, height = HEIGHT, units = 'cm')

g <- generate_page(ORDER_PAGE_2, plot_roll_simulated, c('Standardized residuals', 'Symmetric t copula standardized residuals'))
ggsave('output/rollingCorrelations/rolling_simulated2.png', g, width = WIDTH, height = HEIGHT, units = 'cm')