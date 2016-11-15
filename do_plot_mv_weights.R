
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

SAMPLE_MODEL = 'results_full_sample'

load(sprintf('data/derived/mv/%s_%s.RData', SAMPLE_MODEL, MODEL_NAME_1))
sample_results1 <- results

load(sprintf('data/derived/mv/%s_%s.RData', SAMPLE_MODEL, MODEL_NAME_2))
sample_results2 <- results

rm(results)


# date --------------------------------------------------------------------

load('data/derived/weekly-estim.RData')


# weightplot  -----------------------------------------------------------------

tutti <- bind_rows(results1 = gather(data.frame(Date = df.estim$Date[2:2766], results1$weights), 'Factor', 'Weight', -1),
                   results2 = gather(data.frame(Date = df.estim$Date[2:2766], results2$weights), 'Factor', 'Weight', -1),
                   sample_results1 = gather(data.frame(Date = df.estim$Date[2:2766], sample_results1$weights), 'Factor', 'Weight', -1),
                   sample_results2 = gather(data.frame(Date = df.estim$Date[2:2766], sample_results2$weights), 'Factor', 'Weight', -1),
                   .id = 'Model'
                   )
tutti$Factor <- factor(tutti$Factor, levels = c('HML','CMA','Mkt.RF','SMB','RMW', 'Mom'))
tutti$Model[tutti$Model == 'results1'] = MODEL_NAME_1
tutti$Model[tutti$Model == 'results2'] = MODEL_NAME_2
tutti$Model[tutti$Model == 'sample_results1'] = sprintf('Sample %s', MODEL_NAME_1)
tutti$Model[tutti$Model == 'sample_results2'] = sprintf('Sample %s', MODEL_NAME_2)

tutti <- tutti %>% 
  group_by(Model, Factor) %>% 
  mutate(ma = rollapply(Weight, 52, mean, align = 'right', fill = NA))

ggplot(tutti, aes(x = Date, y = ma, color = Model)) +
    facet_grid(Factor ~ ., switch = 'y') +
    geom_line()+
    theme_Publication()+
    theme(strip.background = element_blank())+
    theme(panel.margin = unit(1, "lines"))+
    theme(legend.direction = 'vertical')+
    theme(legend.key.size = unit(0.75, 'lines'))+
    scale_colour_Publication()+
    coord_cartesian(ylim = c(0, 0.60))+  
    scale_y_continuous(labels = scales::percent, breaks = c(0, 0.30, 0.60))+
    ylab('Smoothed weight (1-year moving average)')


ggsave(sprintf('output/mv/Weights_%s_%s.png', MODEL_NAME_1, MODEL_NAME_2),
       width = 8,
       height = 18,
       units = 'cm',
       limitsize = FALSE)
