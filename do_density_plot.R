#' Look at density plots of GARCH residuals - is there dependency left?
#' Returns a grid of density plots

# Libraries ----
library(ggplot2)
library(GGally)
library(dplyr)

# Reset workspace and load residuals from GARCH ----
rm(list = ls())
load('data/derived/garch_stdres.RData')

# Do and save density plots ----
jpeg(filename = "output/densityGARCHresiduals.jpeg", height = 8.3, width = 11.7, units = 'in', res = 300)
g <- df.stdres %>% select(-Date) %>%
  ggpairs(lower = list(
    continuous = function(data, mapping, ...) {
      ggplot(data = data, mapping = mapping)+
        geom_density2d(...)
      }
    ),
    diag = 'blank'
)
print(g)
dev.off()
