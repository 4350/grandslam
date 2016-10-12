#' Look at density plots of GARCH residuals - is there dependency left?
#' Returns a grid of density plots

# Libraries ----
library(ggplot2)
library(GGally)
library(dplyr)
library(ellipse)
library(mvtnorm)
library(ggthemes)
library(devtools)
library(extrafont)
load_all('wimbledon')

# Reset workspace and load residuals from GARCH ----
rm(list = ls())
load('data/derived/garch_stdres.RData')

# Do and save density plots ----
jpeg(filename = "output/densityGARCHresiduals.jpeg", height = 12, width = 14, units = 'cm', res = 300)
g <- df.stdres %>% select(-Date) %>%
  ggpairs(lower = list(
    continuous = function(data, mapping, ...) {
      ggplot(data = data, mapping = mapping)+
        geom_density2d(..., size = 0.3)+
        theme_Publication()+
        coord_cartesian(xlim =c(-2.5, 2.5), ylim =c(-2.5,2.5))+
        theme(axis.text = element_text(size = rel(0.6), colour = "grey30")) +
        scale_colour_Publication()
      }
    ),
    diag = 'blank',
    upper = list(
      continuous = function(data, mapping, ...) {
        ggplot(data = data, mapping = mapping)+
          geom_bin2d(..., bins = 45)+
          stat_ellipse(..., type = 'norm', level = 0.9995, size = 0.3)+
          theme_Publication()+
          scale_colour_Publication()
      }
    ),
    columnLabels = rep("", 6),
    showStrips = NULL
  )
print(g)
dev.off()