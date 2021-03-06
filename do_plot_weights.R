
# Setup --------------------------------------------------------------------
rm(list = ls())
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)
library(gtable)
library(ggfortify)
library(zoo)
library(gridExtra)
library(devtools)
library(extrafont)

load_all('wimbledon')


# Inputs ------------------------------------------------------------------

# COMPARE
# NAME = 'CDB_MV'
# MODEL = 'mv/results_full_dynamic_std_10000'
# SAMPLE_MODEL = 'cdb/constrOptim_q5_full_dynamic_std_10000'

MV
NAME = 'MV'
MODEL = 'mv/results_full_dynamic_std_10000'
SAMPLE_MODEL = 'mv/results_sample'

#CDB

# NAME = 'CDB'
# MODEL = 'cdb/constrOptim_q5_full_dynamic_std_10000'
# SAMPLE_MODEL = 'cdb/constrOptim_q5_full_dynamic_std_10000' # is not active when NAME = 'CDB'

MODEL_NAME_1 = '5F'
MODEL_NAME_2 = '5F_EXCL_RMW' # not active when NAME = 'CDB_MV'

LABELS = c("Five-factor (model)", #model 1
           "Five-factor excl. RMW (model)", #model-2
           "Five-factor (sample)", #sample-1
           "Five-factor excl. RMW (sample)" #sample-2
           )


# Loading -----------------------------------------------------------------

load(sprintf('data/derived/%s_%s.RData', MODEL, MODEL_NAME_1))
results1 <- results

load(sprintf('data/derived/%s_%s.RData', MODEL, MODEL_NAME_2))
results2 <- results

load(sprintf('data/derived/%s_%s.RData', SAMPLE_MODEL, MODEL_NAME_1))
sample_results1 <- results

load(sprintf('data/derived/%s_%s.RData', SAMPLE_MODEL, MODEL_NAME_2))
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

if(NAME == 'CDB') {
  tutti <- tutti %>%
    filter(Model %in% c(MODEL_NAME_1, MODEL_NAME_2))
}

if(NAME == 'CDB_MV') {
  tutti <- tutti %>%
    filter(Model %in% c(MODEL_NAME_1, sprintf('Sample %s', MODEL_NAME_1)))
}

g <- ggplot(tutti, aes(x = Date, y = ma, color = Model, linetype = Model)) +
    facet_grid(Factor ~ ., switch = 'y') +
    geom_line()+
    theme_Publication()+
    theme(strip.background = element_blank())+
    theme(panel.margin = unit(1, "lines"))+
    theme(legend.direction = 'vertical')+
    theme(legend.position = 'none')+
    scale_colour_Publication()+
    coord_cartesian(ylim = c(0, 0.60))+  
    scale_y_continuous(labels = scales::percent, breaks = c(0, 0.30, 0.60))+
    scale_colour_manual(labels = LABELS, values = c("#386cb0","#fdb462","#386cb0","#fdb462"))+
    scale_linetype_manual(labels = LABELS, values = c("solid","solid","dashed","dashed"))+
    theme(legend.key.size = unit(0.75, 'lines'))+
    theme(legend.key.width = unit(0.6, 'cm'))+
    ylab('Smoothed weight (1-year moving average)')


# And draw for appendix part

# First legend
g_legend <- g+theme(legend.position = 'bottom')+
  scale_colour_manual(labels = LABELS, values = c("#386cb0","#fdb462","#386cb0","#fdb462"))

g_legend = gtable_filter(ggplotGrob(g_legend), "guide-box") 

# Then graphs

g <- 
  grid.arrange(
    g,
    g_legend,
    nrow = 2,
    heights = c(19,2)
  )

ggsave(sprintf('output/weights/Weights_%s_%s_%s.png', NAME, MODEL_NAME_1, MODEL_NAME_2),
       g,
       width = 7.9,
       height = 20,
       units = 'cm',
       limitsize = FALSE)
