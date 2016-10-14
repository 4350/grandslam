#' Estimates threshold and rolling corrs.
#' Currently implementing daily return series. Has previously run weekly
#' series for threshold correlation. The thCorrelationList is saved in derived
#' data for further use in simulated threshold graphs.
#' 


# Libraries ---------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(extrafont)
library(zoo)
library(scales)
library(devtools)

# Reset workspace and load return data
rm(list = ls())
load('data/derived/weekly-estim.RData')
load('data/derived/garch_stdres.RData')


# Load functions
load_all('wimbledon')

# Threshold correlations --------------------------------------------------

# Loop for each value factors to get threshold correlation for all factors
# and save to the list created above

thCorrList.ret <- th_corr(df = df.estim %>% select(-Date),
                          simulatetoggle = 0
                          )

thCorrList.res <- th_corr(df = df.stdres %>% select(-Date),
                          simulatetoggle = 0
                          )


# Bind to one df for plot
plotdf.ret <- bind_rows(thCorrList.ret$HML, thCorrList.ret$RMW, thCorrList.ret$CMA)
plotdf.res <- bind_rows(thCorrList.res$HML, thCorrList.res$RMW, thCorrList.res$CMA)


# Do threshold plot and save ----------------------------------------------

.plot_th_corr <- function(plotdf.ret, plotdf.res, COLFACTORS, OUTNAME) {
  # Select the column factors for plot this plot
  plotdf.ret <- plotdf.ret[plotdf.ret$order == COLFACTORS,]
  plotdf.res <- plotdf.res[plotdf.res$order == COLFACTORS,]
  
  # Then do plot
  g <- ggplot(plotdf.ret, aes(x = qs, y = value)
  ) +
    geom_ribbon(aes(ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
                fill = 'grey10',
                alpha = 0.1
    ) +
    geom_ribbon(aes(ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
                data = plotdf.res,
                fill = 'grey10',
                alpha = 0.1
    ) +
    geom_line(aes(color = 'Return series')) +
    geom_line(aes(x = qs, y = value, color = 'Residual series'), data = plotdf.res) +
    theme_Publication() +
    scale_colour_Publication() +
    ylab('Correlation') +
    xlab('Quantiles') +
    coord_cartesian(xlim = c(0,1), ylim = c(-0.5, 1)) + 
    scale_x_continuous(labels = scales::percent) +
    facet_grid(order ~ order2)
    #ggtitle('Threshold correlations of weekly data (95% confidence bounds)') + 
    #theme(axis.text = element_text(size = rel(0.6), colour = "grey30")) 
  
  # Save plot
  OUTPATH <- 'output/thresholdCorrelations/threshold_%s.png'
  ggsave(sprintf(OUTPATH, OUTNAME), 
    g, device = 'png', width = 14, height = 16, units = 'cm'
    )

}

# Do the plots ------------------------------------------------------------


.plot_th_corr(plotdf.ret, plotdf.res, c('Mkt.RF', 'SMB','Mom'), 'Nonvalue')
.plot_th_corr(plotdf.ret, plotdf.res, c('HML','RMW','CMA'), 'Value')


# Rolling correlations ----------------------------------------------------

# Create list to hold correlation sets
rollCorrList.ret = roll_corr(df = df.estim %>% select(-Date), 
                                   df.date = df.estim %>% select(Date), 
                                   window = 45
                                   )


# Bind to one df for plot
plotdf.ret <- bind_rows(rollCorrList.ret$HML, rollCorrList.ret$RMW, rollCorrList.ret$CMA)

# Do roll plot and save ----------------------------------------------
g <- ggplot(plotdf.ret, aes(x = Date, y = value)
) +
  geom_ribbon(aes(ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
              fill = 'grey10',
              alpha = 0.1
  ) +
  geom_line(aes(color = 'Return series')) +
  theme_Publication() +
  scale_colour_Publication() +
  ylab('Correlation') +
  xlab('Quantiles') +
  scale_x_date(date_labels = "%y") +
  coord_cartesian(ylim = c(-1, 1)) + 
  facet_grid(order ~ order2)
  #ggtitle('Rolling 45-week correlations (95% confidence bounds)')#+
  #theme(axis.text = element_text(size = rel(0.6), colour = "grey30")) 

ggsave('output/rollingCorrelations/rolling45.png', g, device = 'png', width = 14, height = 18, units = 'cm', limitsize = F)