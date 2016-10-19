
# Libraries and setup -----------------------------------------------------

library(tidyr)
library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(extrafont)
library(zoo)
library(scales)
library(devtools)

load_all('wimbledon')

rm(list = ls())

load('data/derived/usrec-weekly.RData')
load('data/derived/weekly-estim.RData')
df_dates <- as.data.frame(df.estim[,1])
df <- as.data.frame(df.estim[,-1])
rm(df.estim)

#  ------------------------------------------------------------------------

below_above_quantile <- function(df, df_dates, q){
  # Calculate a list of return, below and above index for each column
  df_list <- lapply(df, function(f, df_dates, q) {
    data.frame(
      Date = df_dates,
      ret = f,
      below = as.numeric(f < quantile(f, probs = q)),
      above = as.numeric(f > quantile(f, probs = (1-q)))
    )
  }, df_dates = df_dates, q = q)
}

list_below_above <- below_above_quantile(df, df_dates, 0.10)

# Improve that df with dates and the factor name and collect all
for(f in 1:6){
  list_below_above[[f]]$factor = colnames(df)[f]
}
plot_below_above <- bind_rows(list_below_above)
plot_below_above$order <- factor(plot_below_above$factor, colnames(df))

# Do plot with NBER dummy --------------------------------------------

g <- ggplot(plot_below_above) +
  geom_bar(stat = 'identity', aes(x = Date, y = ret*above, color = 'Top 10% percentile returns')) +
  geom_bar(stat = 'identity', aes(x = Date, y = ret*below, color = 'Bottom 10% percentile returns')) +
  geom_ribbon(data = usrec, 
              mapping = aes(x = Date, ymin = -1, 
                            ymax = -1 + 2 * recdummy, 
                            linetype = NA, 
                            fill = 'sienna2'
              ),
              fill = 'sienna2',
              alpha = 0.3
  ) +
  theme_Publication() +
  scale_colour_Publication() +
  scale_y_continuous(labels = scales::percent)+
  ylab('Return') +
  xlab('Year') +
  scale_x_date(date_labels = "%y") +
  coord_cartesian(ylim = c(-.20, .20)) +
  facet_grid(order ~ .)


ggsave('output/below_above_quantile.png', g, device = 'png', width = 14, height = 18, units = 'cm', limitsize = F)

