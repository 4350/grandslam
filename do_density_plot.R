#' Look at density plots of GARCH residuals - is there dependency left?
#' Returns a grid of density plots

# Libraries ----
library(ggplot2)
library(GGally)
library(dplyr)
library(ellipse)
library(mvtnorm)


# Reset workspace and load residuals from GARCH ----
rm(list = ls())
load('data/derived/garch_stdres.RData')

# Set colors and seed
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
set.seed(17)

# Do and save density plots ----
jpeg(filename = "output/densityGARCHresiduals.jpeg", height = 8.3, width = 11.7, units = 'in', res = 300)
g <- df.stdres %>% select(-Date) %>%
  ggpairs(lower = list(
    continuous = function(data, mapping, ...) {
      ggplot(data = data, mapping = mapping)+
        geom_density2d(...)
      }
    ),
    diag = 'blank',
    upper = list(
      continuous = function(data, mapping, ...) {
        ggplot(data = data, mapping = mapping)+
          geom_bin2d(..., bins = 45)+
          stat_ellipse(..., type = 'norm', level = 0.9995, color = cbbPalette[2])
      }
    )
)
print(g)
dev.off()