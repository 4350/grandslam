library(mvtnorm)
library(tidyr)
library(dplyr)
library(ggplot2)
rm(list = ls())

load('data/derived/weekly-full.RData')
load('data/derived/garch_stdres.RData')

df <- df.stdres
X <- df
Date <- X$Date
X <- df[, -1]

U <- apply(X, 2, function(x) ecdf(x)(x))
distance_to_median <- U - 0.50

df_plot <- data.frame(
  Date = Date,
  distance_to_median
)
df_plot <- gather(df_plot, Factor, Distance, -1)

df_pos <- filter(df_plot, Distance >= 0)
df_neg <- filter(df_plot, Distance <  0)

ggplot(df_pos) +
  geom_area(data = df_pos, aes(x = Date, y = Distance, color = Factor))
  # geom_area(data = df_neg, aes(x = Date, y = Distance, color = Factor))