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

load('data/derived/threshold/th_corr_returns.RData')
load('data/derived/threshold/th_corr_stdres.RData')
load('data/derived/threshold/th_corr_constant_simulated.RData')
load('data/derived/garch/model_GARCH_chosen_stdres.RData')

# Set ordering for plots

MODEL_ORDER <- c('stdres','returns','norm','std','ghst')

ORDER_PAGE_1 <- list(
  c('Mkt.RF', 'HML'),
  c('Mkt.RF', 'CMA'),
  c('Mkt.RF', 'RMW'),
  c('SMB', 'HML'),
  c('SMB', 'CMA'),
  c('SMB', 'RMW')
)

ORDER_PAGE_2 <- list(
  c('Mom', 'HML'),
  c('Mom', 'CMA'),
  c('Mom', 'RMW'),
  c('HML', 'CMA'),
  c('HML', 'RMW'),
  c('CMA', 'RMW')
)

WIDTH = 16
HEIGHT = 20

# Print the pages with plots ----------------------------------------------

# Main section
g <- generate_page(ORDER_PAGE_1, plot_th_main)
ggsave('output/thresholdCorrelations/threshold1.png', g, width = WIDTH, height = HEIGHT, units = 'cm')

g <- generate_page(ORDER_PAGE_2, plot_th_main)
ggsave('output/thresholdCorrelations/threshold2.png', g, width = WIDTH, height = HEIGHT, units = 'cm')

# Appendix section
g <- generate_page(ORDER_PAGE_1, plot_th_appendix, c('ARMA-GARCH standardized residuals', 'Returns'))
ggsave('output/thresholdCorrelations/appendix_threshold_1.png', g, width = WIDTH, height = HEIGHT, units = 'cm')

g <- generate_page(ORDER_PAGE_2, plot_th_appendix, c('ARMA-GARCH standardized residuals', 'Returns'))
ggsave('output/thresholdCorrelations/appendix_threshold_2.png', g, width = WIDTH, height = HEIGHT, units = 'cm')

# Simulated section
g <- generate_page(ORDER_PAGE_1, plot_th_simulated, c('Standardized residuals', 'Normal','Symmetric t','Skewed t'))
ggsave('output/thresholdCorrelations/threshold_simulated_1.png', g, width = WIDTH, height = HEIGHT, units = 'cm')

g <- generate_page(ORDER_PAGE_2, plot_th_simulated, c('Standardized residuals', 'Normal','Symmetric t','Skewed t'))
ggsave('output/thresholdCorrelations/threshold_simulated_2.png', g, width = WIDTH, height = HEIGHT, units = 'cm')  

# Explanatory plot
g <- plot_th_explain(c('Mkt.RF','HML'))

  