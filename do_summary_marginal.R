#' Create marginal summary statistics table and plots

# Libraries ----
library(dplyr)
library(fBasics)
library(tidyr)
library(ggplot2)
library(scales)
library(ggfortify)
library(gridExtra)
library(devtools)
library(extrafont)
library(stargazer)
library(WeightedPortTest)

# Reset workspace and load return data ----
rm(list = ls())
load('data/derived/weekly-estim.RData')

# Load functions ----
load_all('wimbledon')


# LB, ARCH LM -------------------------------------------------------------

ret_LB_5 <- apply(df.estim[,-1], 2, function(ret) Weighted.Box.test(ret, lag = 5, type = "Ljung-Box", 
                                                                    weighted = TRUE)$p.value)
ret_LB_10 <- apply(df.estim[,-1], 2, function(ret) Weighted.Box.test(ret, lag = 10, type = "Ljung-Box", 
                                                                     weighted = TRUE)$p.value)

sqr_ret_LB_5 <- apply(df.estim[,-1], 2, function(ret_sqr) Weighted.Box.test(ret_sqr, lag = 5, type = "Ljung-Box",
                                                                            weighted = TRUE,
                                                                            sqrd.res = TRUE)$p.value)
sqr_ret_LB_10 <- apply(df.estim[,-1], 2, function(ret_sqr) Weighted.Box.test(ret_sqr, lag = 10, type = "Ljung-Box",
                                                                             weighted = TRUE,
                                                                             sqrd.res = TRUE)$p.value)

#' Summary statistics table

summary_table <- df.estim %>%
    dplyr::select(-Date) %>%
    select(Mkt.RF, SMB, Mom, HML, CMA, RMW) %>%
    basicStats() %>%
    .[c('nobs','Maximum','Minimum','Mean','Median','Stdev','Skewness','Kurtosis'),]

summary_table <- rbind(summary_table,
                       'Return LB [5] p-value' = ret_LB_5,
                       'Return LB [10] p-value' = ret_LB_10,
                       'Squared return LB [5] p-value' = sqr_ret_LB_5,
                       'Squared return LB [10] p-value' = sqr_ret_LB_10)

write.table(summary_table, file = 'output/MarginalStats/summaryTable.Estim.csv')
#stargazer(summary_table, summary = FALSE)


# QQ Plots ----------------------------------------------------------------

plot_df <- df.estim[,-1] %>%
  gather('factor','value')

plot_df$factor <- factor(plot_df$factor, levels = c('Mkt.RF','SMB','HML','CMA','RMW','Mom'))

line_df <- plot_df %>%
  group_by(factor) %>%
  summarize(q25    = quantile(value,0.25),
            q75    = quantile(value,0.75),
            norm25 = qnorm(0.25),
            norm75 = qnorm(0.75),
            slope  = (q25 - q75) / (norm25 - norm75),
            intercept    = q25 - slope * norm25
            ) %>%
  select(factor,slope,intercept)

g <- ggplot(plot_df, 
       aes(
         sample = value
         )
       )+
  stat_qq(size = 0.5)+
  geom_abline(
    aes(intercept = intercept, slope = slope),
    linetype = 2,
    data = line_df
  )+
  coord_cartesian(xlim = c(-4,4))+
  annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=0.25)+
  facet_wrap( ~ factor, nrow = 2, ncol = 3)+
  theme_Publication() +
  theme(strip.background = element_blank())+
  scale_colour_Publication() +
  ylab('Sample returns') +
  xlab('Theoretical quantiles') +
  theme(legend.position = 'none')+
  scale_y_continuous(labels = scales::percent)

ggsave('output/MarginalStats/qq_returns.png', g, device = 'png', width = 9, height = 5.8, units = 'cm', limitsize = F)
  

# Corr matrix data --------------------------------------------------------

cor_matrix <- cor(df.estim[,-1] %>% select(Mkt.RF, SMB, HML, CMA, RMW, Mom))
cor_matrix[upper.tri(cor_matrix)] <- NA
stargazer(cor_matrix)
stargazer(cor_matrix, type = 'text', digits = 2)

# Marginal plots ----
varlist = list('Mkt.RF','HML','SMB','Mom','RMW','CMA')

lapply(varlist,
       function(varlist) summary.plots(df.estim, 'Estim', varlist))

# Cumulative plots
summary.cumretplots(df.estim, 'Estim')