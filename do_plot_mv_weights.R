
# Setup --------------------------------------------------------------------
rm(list = ls())
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)
library(ggfortify)
library(zoo)
library(gridExtra)
library(devtools)
library(extrafont)

load_all('wimbledon')

MODEL = 'results_full_dynamic_std_10000'
MODEL_NAME_1 = '5F_EXCL_HML'
MODEL_NAME_2 = '5F'

load(sprintf('data/derived/mv/%s_%s.RData', MODEL, MODEL_NAME_1))
results1 <- results

load(sprintf('data/derived/mv/%s_%s.RData', MODEL, MODEL_NAME_2))
results2 <- results

rm(results)

# weightplot  -----------------------------------------------------------------

tutti <- bind_rows(results1 = gather(data.frame(t = 1:2765, results1$weights), 'Factor', 'Weight', -1),
                   results2 = gather(data.frame(t = 1:2765, results2$weights), 'Factor', 'Weight', -1), 
                   .id = 'Model'
                   )
tutti$Factor <- factor(tutti$Factor, levels = c('HML','CMA','Mkt.RF','SMB','RMW'))
tutti$Model[tutti$Model == 'results1'] = MODEL_NAME_1
tutti$Model[tutti$Model == 'results2'] = MODEL_NAME_2

tutti <- tutti %>% group_by(Model, Factor) %>% mutate(ma = rollapply(Weight, 52, mean, align = 'right', fill = NA))

grid.arrange(
  ggplot(tutti, aes(x = t, y = Weight, color = Model)) +
    facet_grid(Factor ~ ., switch = 'y') +
    geom_line()+
    theme_Publication()+
    theme(strip.background = element_blank())+
    theme(legend.position = 'none')+
    scale_colour_Publication(),
  
  ggplot(tutti, aes(x = t, y = ma, color = Model)) +
    facet_grid(Factor ~ ., switch = 'y') +
    geom_line()+
    theme_Publication()+
    theme(strip.background = element_blank())+
    theme(legend.position = 'none')+
    scale_colour_Publication()+
    ylab('Smoothed weight (1-year moving average)')
    ,
  ncol =2, nrow = 1
)

