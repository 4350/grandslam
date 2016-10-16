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

FACTORS <- c("Mkt.RF", "HML", "SMB", "Mom", "RMW", "CMA")

# Empircal data
df.empirical <- df.stdres %>% select(-Date)
rm(df.stdres)

# Calculate threshold correlation data frame lists ------------------------

# Get simulated df lists
#sim <- get.simulation.results('data/derived/stdresid/constant_gauss')

th.corr.gauss <- do.th.sim('data/derived/stdresid/constant_gauss')
th.corr.ght   <- do.th.sim('data/derived/stdresid/constant_ght')
th.corr.ghskt <- do.th.sim('data/derived/stdresid/constant_ghskt')

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


# Unconditional ----------------------------------------------------------

stdresid <- read.csv('data/derived/stdresid/constant_gauss/100000_1.csv')
stdresid <- data.frame(stdresid[, -1])
colnames(stdresid) <- FACTORS

# Do Good Threshold Plot -------------------------------------------------

series <- c('HML', 'RMW', 'CMA')
# row_series <- c('HML', 'RMW', 'CMA')
# 
# unconditional <- data.frame(factor = FACTORS, unconditional)
# unconditional <- gather(unconditional, factor2, value, -factor, factor_key = TRUE)

plotdf <- dplyr::filter(plotdf, order %in% series)

ggplot(plotdf,
       aes(x = qs, y = value)) +
  geom_ribbon(aes(ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
              fill = 'grey10',
              alpha = 0.1
  ) +
  geom_line(aes(color = "Standardized Residuals")) +
  geom_line(aes(y = gauss, color = 'Simulated Gaussian')) +
  geom_line(aes(y = ght, color = "Simulated Student's t")) +
  geom_line(aes(y = ghskt, color = "Simulated Skewed Student's t")) +
  theme_Publication() +
  scale_colour_Publication() +
  ylab('Correlation') +
  xlab('Quantiles') +
  coord_cartesian(xlim = c(0,1), ylim = c(-0.5, 1)) + 
  scale_x_continuous(labels = scales::percent) +
  facet_grid(order2 ~ order)

