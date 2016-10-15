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


# Function ----------------------------------------------------------------

# Quick fix to add standard correlations to a graph data frame

.append_standard_corr <- function(plot_df, corr_df) {
  # Check which unique pairs to do correlations between
  label_pairs <- unique(select(plot_df, c(order2, order)))
  
  # Get these correlations
  standard_corr <- apply(label_pairs, 1, function(var_pair, corr_df) {
    cor(corr_df[,var_pair])[1,2]
  },
  corr_df = corr_df
  )
  
  # Make data frame from the labels and corrs
  data.frame(
    standard_corr = round(standard_corr,2),
    order2 = label_pairs[,1],
    order = label_pairs[,2]
  )
  
}

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

.plot_th_corr <- function(plotdf.ret, plotdf.res, df.labels, ROWFACTORS, OUTNAME) {
  # Select the column factors for plot this plot
  plotdf.ret <- plotdf.ret[plotdf.ret$order == ROWFACTORS,]
  plotdf.res <- plotdf.res[plotdf.res$order == ROWFACTORS,]
  
  df.labels <- df.labels %>% dplyr::filter(order %in% ROWFACTORS)
  
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
    #annotate("rect", xmin = 0.375, xmax = 0.625, ymin = -0.50, ymax = -0.25, alpha = 0.8, fill = 'grey80')+
    geom_text(data = df.labels, aes(x = 0.50, y = -0.45, label = paste("r = ", standard_corr)), family = 'Minion Pro', size = 3, parse = F)+
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

# Quick fix to append standard correlations to graphs
df.labels <- .append_standard_corr(plotdf.ret, df.estim)

.plot_th_corr(plotdf.ret, plotdf.res, df.labels, c('Mkt.RF', 'SMB','Mom'), 'Nonvalue')
.plot_th_corr(plotdf.ret, plotdf.res, df.labels, c('HML','RMW','CMA'), 'Value')


# Rolling correlations ----------------------------------------------------

# Create list to hold correlation sets
rollCorrList.ret = roll_corr(df = df.estim %>% select(-Date), 
                                   df.date = df.estim %>% select(Date), 
                                   window = 45
                                   )


# Bind to one df for plot
plotdf.ret <- bind_rows(rollCorrList.ret$HML, rollCorrList.ret$RMW, rollCorrList.ret$CMA)

# Quick fix to append standard correlations to graphs
df.labels <- .append_standard_corr(plotdf.ret, df.estim)

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
  xlab('Year') +
  scale_x_date(date_labels = "%y") +
  coord_cartesian(ylim = c(-1, 1)) + 
  #annotate("rect", xmin = as.Date('1986-01-01'), xmax = as.Date('1994-01-01'), ymin = -0.95, ymax = -0.5, alpha = 0.8, fill = 'grey80')+
  #geom_text(data = df.labels, aes(x = as.Date('1990-01-01'), y = -0.725, label = paste('r = ',standard_corr)), family = 'Minion Pro', size = 3, parse = FALSE)+
  geom_text(data = df.labels, aes(x = as.Date('2010-01-01'), y = -0.90, label = paste('r = ',standard_corr)), family = 'Minion Pro', size = 3, parse = FALSE)+
  facet_grid(order ~ order2)
  #ggtitle('Rolling 45-week correlations (95% confidence bounds)')#+
  #theme(axis.text = element_text(size = rel(0.6), colour = "grey30")) 

ggsave('output/rollingCorrelations/rolling45.png', g, device = 'png', width = 14, height = 18, units = 'cm', limitsize = F)