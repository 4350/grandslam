#'
#'


# Libraries ----
library(rugarch)
library(parallel)
library(devtools)
library(dplyr)
load_all('wimbledon')

# Data ----
rm(list = ls())
load('data/derived/model_GARCH.RData')
load('data/derived/weekly-full.RData')
load('data/derived/weekly-estim.RData')

# Get out of sample length
T = nrow(df.estim)
H = nrow(df) - T
rm(df.estim)

df.date <- df %>% select(Date) %>% as.vector()
df <- df %>% select(-Date)

# Best models, as determined before
kGARCHModels <- list(
  Mkt.RF = garch.specgen(0, 0),
  HML = garch.specgen(1, 1),
  SMB = garch.specgen(1, 1),
  Mom = garch.specgen(1, 0),
  RMW = garch.specgen(1, 0),
  CMA = garch.specgen(1, 0) 
)

# Get GARCH fits that preserve the OOS data
garch.fits <- estimate.garch(df, H)

# Make the rolling 1-step ahead forecasts of sigma and the series
fc <- fc_List(garch.fits, H)

# Consolidate series of sigma in and out of sample ----
sigma <- rbind(
  # Get conditional sigma in-sample
  sapply(garch.fits,
         function(fFit) {
           sigma = fFit@fit$sigma
         }),
  # Get conditional sigma "out-of-sample"
  sapply(fc,
         function(ffc) {
           sigma = ffc@forecast$sigmaFor
         })
)

# TESTING SECTION ----------------------------------------------------------------------------------

# Kupiec and Christoffersen tests using rugarch for simplicity
kupiec.05 <- kupiec_test.var(.05, df, sigma, T, kGARCHModels)
kupiec.01 <- kupiec_test.var(.01, df, sigma, T, kGARCHModels)

# PLOT SECTION -------------------------------------------------------------------------------------
# Libraries for plot ----
library(tidyr)
library(ggplot2)
library(extrafont)
library(scales)

# Make tidy empirical data for plot ----
df.emp <- gather(df[(T+1):nrow(df), ], 'factor','value', 1:6)
df.emp$h <- rep(seq(1,H), ncol(df))
df.emp$order <- factor(df.emp$factor, names(kGARCHModels))


# Get tidy data sets for plots
df.var05 <- hist.sim.tidy.var(.05, df, sigma, T, kGARCHModels)
df.var01 <- hist.sim.tidy.var(.01, df, sigma, T, kGARCHModels)

# Plot that
g <- ggplot(df.emp, aes(x = h, y = value, group = factor))+
  geom_line(aes(color = 'Realized return'))+
  geom_line(aes(x = h, y = VaR, color = '5% HS-HW VaR'), data = df.var05)+
  geom_line(aes(x = h, y = VaR, color = '1% HS-HW VaR'), data = df.var01)+
  theme_Publication()+
  ylab('')+
  xlab('Horizon')+
  scale_y_continuous(labels = percent)+
  facet_grid(. ~ order)+
  ggtitle("Value-at-Risk out-of-sample using HS-HW method")
ggsave(file = paste('output/VaR/VaRHSHW', 'T', T, 'jpeg', sep = '.'), g, width = 16.6, height = 11.7, units = 'in')



