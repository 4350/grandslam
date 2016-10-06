#' Estimates threshold corrs from simulated data and plots
#' this together with estimated threshold corrs on empirical data
#'
#' All threshold correlation applied to (GARCH) residuals
#'

# Load libraries, functions, data -----------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(extrafont)
library(zoo)
library(scales)
library(devtools)
load_all('wimbledon')

# Reset workspace and load data for each distribution
rm(list = ls())
load('data/derived/garch_stdres.RData')

# Empircal data
df.empirical <- df.stdres %>% select(-Date)
rm(df.stdres)

# Calculate threshold correlation data frame lists ------------------------

# Get simulated df lists
th.corr.gauss <- do.th.sim('data/derived/simulation/constant_gauss/')
th.corr.ght   <- do.th.sim('data/derived/simulation/constant_ght/')
th.corr.ghskt <- do.th.sim('data/derived/simulation/constant_ghskt/')

# Get empirical df list
th.corr.empirical <- th_corr(df.empirical, 0)

# Consolidate all dfs to one list
allThCorrList <- list()

# Set up factor groups ----
factors.value <- c("HML", "RMW", "CMA")
factors.all <- c("Mkt.RF", "HML", "SMB", "Mom", "RMW", "CMA")

for (values in factors.value) {
  out.df <- data.frame(
    # The empirical imports coef, lb, ub, quantiles and ordering of factors
    th.corr.empirical[[values]],

    # All other distributions only import coefficients
    gauss = th.corr.gauss[[values]]$simcoef,
    ght = th.corr.ght[[values]]$simcoef,
    ghskt = th.corr.ghskt[[values]]$simcoef
  )
  allThCorrList[[values]] <- out.df
}

# Bind all value dfs to one df for plot
plotdf <- bind_rows(allThCorrList$HML, allThCorrList$RMW, allThCorrList$CMA)


# Do threshold plot and save ----------------------------------------------
g <- ggplot(plotdf, aes(x = qs, y = value)
            ) +
  geom_ribbon(aes(ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
              fill = 'grey10',
              alpha = 0.1
              ) +
  geom_line(aes(color = "Empirical distribution")) +
  geom_line(aes(x = qs, y = gauss, color = "Simulated Gaussian"), linetype = 2) +
  geom_line(aes(x = qs, y = ght, color = "Simulated Student-t"), linetype = 3) +
  geom_line(aes(x = qs, y = ghskt, color = "Simulated Skewed Student-t"), linetype = 4) +
  theme_Publication() +
  ylab('') +
  xlab('Quantiles') +
  coord_cartesian(xlim = c(0, 1), ylim = c(-0.5, 1)) +
  scale_x_continuous(labels = scales::percent) +
  facet_grid(order2 ~ order) +
  ggtitle('Threshold correlations of weekly data (95% confidence bounds)')

ggsave(file = 'output/thresholdCorrelations/WeeklySimulated.jpeg', g, width = 16.6, height = 11.7, units = 'in')
