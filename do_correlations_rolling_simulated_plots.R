rm(list = ls())

library(ggplot2)
library(foreach)
library(devtools)
library(dplyr)
library(gridExtra)
load_all('wimbledon')

load('data/derived/garch/model_GARCH_chosen_stdres.RData')
load('data/derived/correlations_rolling_simulated_std.RData')

Date <- df.stdres$Date

# Clean up empirical
stdresid_roll_corr = roll_corr(
  df = df.stdres %>% select(-Date), 
  df.date = df.stdres %>% select(Date), 
  window = 52
)
stdresid_roll_corr <- bind_rows(stdresid_roll_corr) %>%
  select(Date, order, order2, value, lb, ub)
colnames(stdresid_roll_corr) <- c('Date', 'factor1', 'factor2', 'value', 'lb', 'ub')

plot_pairs <- function(pairs) {
  # Get the rolling simulated correlations for these pairs
  simulated_roll_corr_df <- bind_rows(lapply(pairs, function(pair) {
    corr <- sapply(simulated_roll_corr, function(c) c[pair[1], pair[2]])
    data.frame(
      Date = tail(Date, length(simulated_roll_corr)),
      factor1 = pair[1],
      factor2 = pair[2],
      value = corr
    )
  }))
  
  # Generate plots for each pair
  plots <- lapply(pairs, function(pair) {
    ggplot(
      data = stdresid_roll_corr %>% filter(factor1 == pair[1], factor2 == pair[2]),
      aes(x = Date, y = value)
    ) +
      # Uncertainty, empirical stdresid
      geom_ribbon(
        aes(ymin = lb, ymax = ub),
        fill = 'grey10',
        alpha = 0.1
      ) +
      # Rolling correlations, empirical stdresid
      geom_line(
        aes(color = 'Rolling Correlations, Empirical')
      ) +
      
      # Rolling correlations, simulated stdresid
      geom_line(
        data = simulated_roll_corr_df %>% filter(factor1 == pair[1], factor2 == pair[2]),
        aes(x = Date, y = value, color = 'Rolling Correlations, Simulated')
      ) +
      
      # Theming
      theme_Publication() +
      scale_colour_Publication() +
      theme(legend.position = 'none') +
      ylab('Correlation') +
      xlab('Year') +
      scale_x_date(date_labels = "%y") +
      coord_cartesian(
        ylim = c(-1, 1),
        xlim = c(Date[1], tail(Date, 1))
      ) +
      annotate(
        "segment",
        x = tail(Date, 1),
        xend = Date[1],
        y = Inf,
        yend = Inf,
        color = "black",
        lwd = 1
      )+
      ggtitle(sprintf("%s - %s", pair[1], pair[2]))
  })
  
  grid.arrange(
    grobs = plots,
    ncol = 2,
    nrow = 3,
    as.table = FALSE
  )
}

# Plot Pairs

g1 <- plot_pairs(
  list(
    c('Mkt.RF', 'HML'),
    c('Mkt.RF', 'CMA'),
    c('Mkt.RF', 'RMW'),
    c('Mom', 'HML'),
    c('Mom', 'CMA'),
    c('Mom', 'RMW')
  )
)
ggsave(
  'output/rollingCorrelations/rolling_simulated1.png',
  g1,
  width = 16,
  height = 17,
  limitsize = FALSE,
  units = 'cm'
)

#
g2 <- plot_pairs(
  list(
    c('SMB', 'HML'),
    c('SMB', 'CMA'),
    c('SMB', 'RMW'),
    c('HML', 'CMA'),
    c('HML', 'RMW'),
    c('CMA', 'RMW')
  )
)
ggsave(
  'output/rollingCorrelations/rolling_simulated2.png',
  g2,
  width = 16,
  height = 17,
  limitsize = FALSE,
  units = 'cm'
)

# Blue -- empirical
# Yellow -- Simulated


