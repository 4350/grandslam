#' Estimates threshold and rolling corrs.
#' Currently implementing daily return series. Has previously run weekly
#' series for threshold correlation. The thCorrelationList is saved in derived
#' data for further use in simulated threshold graphs.
#' 



# Setup -------------------------------------------------------------------


library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(extrafont)
library(zoo)
library(scales)
library(gtable)
library(devtools)
library(mvtnorm)

# Reset workspace and load return data
rm(list = ls())
load('data/derived/weekly-estim.RData')
ID = 'model_GARCH_chosen_stdres.RData'
load(sprintf('data/derived/garch/%s', ID))


# Load functions
load_all('wimbledon')

# Quick fix to add standard correlations to a graph data frame

.append_standard_corr <- function(plot_df, corr_df) {
  # Check which unique pairs to do correlations between
  label_pairs <- unique(select(plot_df, c(order2, order)))
  
  # Get these correlations
  standard_corr <- apply(label_pairs, 1, function(var_pair, corr_df) {
    cor(corr_df[,var_pair])[1,2]
  },
  corr_df = corr_df
  )
  
  # Make data frame from the labels and corrs
  data.frame(
    standard_corr = round(standard_corr,2),
    order2 = label_pairs[,1],
    order = label_pairs[,2]
  )
  
}

# Threshold correlations --------------------------------------------------

# Loop for each value factors to get threshold correlation for all factors
# and save to the list created above

thCorrList.ret <- th_corr(df = df.estim %>% select(-Date),
                          simulatetoggle = 0
)

thCorrList.res <- th_corr(df = df.stdres %>% select(-Date),
                          simulatetoggle = 0
)


# Bind to one df for plot
plotdf.ret <- bind_rows(thCorrList.ret$HML, thCorrList.ret$RMW, thCorrList.ret$CMA)
plotdf.res <- bind_rows(thCorrList.res$HML, thCorrList.res$RMW, thCorrList.res$CMA)


# Prepare scatter data ----------------------------------------------------

# Return scatter data frame

do_scatter_df <- function(df) {
  factors_value = c('HML','RMW','CMA')
  out_list <- lapply(
    factors_value, 
    function(value, df) {
      
      df <- df[,-1]
      factors_all = c('Mkt.RF','HML','SMB','Mom','RMW','CMA')
      
      bind_rows(
        lapply(
          factors_all, 
          function(other_factor, df) {
            out.df <- data.frame(
              value1 = as.vector(df[,other_factor]),
              value2 = as.vector(df[,value]),
              order = factor(other_factor, levels = factors_all),
              order2 = factor(value, levels = factors_all),
              uniq = paste0(other_factor, value)
            )
            colnames(out.df) <- c('value1','value2','order','order2', 'uniq')
            return(out.df)
          },
          df = df
        )
      )
    }, 
    df = df
  )
  out_df <- bind_rows(out_list)
  # Get breakpoints for x y quantiles
  out_df <- out_df %>% 
    group_by(uniq) %>%
    mutate(
      px10 = quantile(value1, probs = .10),
      py10 = quantile(value2, probs = .10),
      px49 = quantile(value1, probs = .49),
      py49 = quantile(value2, probs = .49),
      px51 = quantile(value1, probs = .51),
      py51 = quantile(value2, probs = .51),
      px90 = quantile(value1, probs = .90),
      py90 = quantile(value2, probs = .90)
    ) %>% ungroup()
}

plotdf.scatter.ret <- do_scatter_df(df.estim)
plotdf.scatter.res <- do_scatter_df(df.stdres)

# Do threshold plot and save for main text ----------------------------------------------

.plot_th_corr <- function(plotdf = plotdf.res,
                          df.labels = df.labels.res, COLFACTORS, ROWFACTORS) {
  # Select the column factors for plot this plot
  plotdf <- plotdf %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  df.labels <- df.labels %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  
  # Then do threshold plot
  g <- ggplot(data = plotdf) +
    geom_ribbon(aes(x = qs, ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
                fill = 'grey10',
                alpha = 0.1
    ) +
    geom_line(aes(x = qs, y = value, colour = 'Threshold correlation')) +
    geom_abline(aes(colour = 'Correlation', intercept = standard_corr, slope = 0),
                colour = 'grey20', size = 0.25, linetype = 2, data = df.labels
    )+
    theme_Publication() +
    scale_colour_Publication()+
    theme(legend.position = 'none')+
    ylab('Correlation') +
    xlab('Quantiles') +
    coord_cartesian(xlim = c(0.10,0.90), ylim = c(-0.50, .75)) + 
    scale_x_continuous(labels = scales::percent, breaks = c(0.10, 0.50, 0.90)) +
    annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=1)+
    ggtitle(sprintf("%s - %s", COLFACTORS, ROWFACTORS))
  
  g
  
}


# Threshold plot including returns ----------------------------------------

.plot_th_corr_incl_ret <- function(plotdf = plotdf.res, plotdf2 = plotdf.ret,
                                   df.labels = df.labels.res, df.labels2 = df.labels.ret,
                                   COLFACTORS, ROWFACTORS) {
  # Select the column factors for plot this plot
  plotdf <- plotdf %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  plotdf2 <- plotdf2 %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  df.labels <- df.labels %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  df.labels2 <- df.labels2 %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  
  # Then do threshold plot
  g <- ggplot(data = plotdf) +
    geom_ribbon(aes(x = qs, ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
                fill = 'grey10',
                alpha = 0.1
    ) +
    geom_line(aes(x = qs, y = value, color = 'Threshold correlation (residuals)')) +
    geom_ribbon(aes(x = qs, ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
                fill = 'grey10',
                alpha = 0.1,
                data = plotdf2
    ) +
    geom_line(aes(x = qs, y = value, color = 'Threshold correlation (returns)'),
              data = plotdf2) +
    geom_abline(aes(slope = 0, intercept = standard_corr,
                    color = 'Correlation (residuals)'), 
                colour = 'grey20', size = 0.25, linetype = 2, data = df.labels
    )+
    geom_abline(aes(slope = 0, intercept = standard_corr,
                    color = 'Correlation (returns)'), 
                colour = 'grey20', size = 0.25, linetype = 3, data = df.labels2,
                show.legend = TRUE)+
    theme_Publication() +
    scale_colour_Publication() +
    theme(legend.position = 'none')+
    ylab('Correlation') +
    xlab('Quantiles') +
    coord_cartesian(xlim = c(0.10,0.90), ylim = c(-0.50, .75)) + 
    scale_x_continuous(labels = scales::percent, breaks = c(0.10, 0.50, 0.90)) +
    annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=1)+
    ggtitle(sprintf("%s - %s", COLFACTORS, ROWFACTORS))
  
  g
  
}

# Do the threshold plots for threshold part ------------------------------------------------------------

# Quick fix to append standard correlations to graphs
df.labels.ret <- .append_standard_corr(plotdf.ret, df.estim)
df.labels.res <- .append_standard_corr(plotdf.res, df.stdres)
  
# Reorder variables
plotdf.res$order <- factor(plotdf.res$order, levels = c('Mkt.RF','SMB','Mom','HML','CMA','RMW'))
plotdf.res$order2 <- factor(plotdf.res$order2, levels = c('Mkt.RF','SMB','Mom','HML','CMA','RMW'))
plotdf.ret$order <- factor(plotdf.ret$order, levels = c('Mkt.RF','SMB','Mom','HML','CMA','RMW'))
plotdf.ret$order2 <- factor(plotdf.ret$order2, levels = c('Mkt.RF','SMB','Mom','HML','CMA','RMW'))

# Draw the pages for results main part

threshold1 <- grid.arrange(
  .plot_th_corr(COLFACTORS = 'Mkt.RF', ROWFACTORS = 'HML'),
  .plot_th_corr(COLFACTORS = 'Mkt.RF', ROWFACTORS = 'CMA'),
  .plot_th_corr(COLFACTORS = 'Mkt.RF', ROWFACTORS = 'RMW'),
  .plot_th_corr(COLFACTORS = 'Mom', ROWFACTORS = 'HML'),
  .plot_th_corr(COLFACTORS = 'Mom', ROWFACTORS = 'CMA'),
  .plot_th_corr(COLFACTORS = 'Mom', ROWFACTORS = 'RMW'),
  ncol = 2,
  nrow = 3,
  as.table = FALSE
)
ggsave('output/thresholdCorrelations/threshold1.png', threshold1, width = 16, height = 17, limitsize = FALSE, units = 'cm')

threshold2 <- grid.arrange(
  .plot_th_corr(COLFACTORS = 'SMB', ROWFACTORS = 'HML'),
  .plot_th_corr(COLFACTORS = 'SMB', ROWFACTORS = 'CMA'),
  .plot_th_corr(COLFACTORS = 'SMB', ROWFACTORS = 'RMW'),
  .plot_th_corr(COLFACTORS = 'HML', ROWFACTORS = 'CMA'),
  .plot_th_corr(COLFACTORS = 'HML', ROWFACTORS = 'RMW'),
  .plot_th_corr(COLFACTORS = 'CMA', ROWFACTORS = 'RMW'),
  ncol = 2,
  nrow = 3,
  as.table = FALSE
)
ggsave('output/thresholdCorrelations/threshold2.png', threshold2, width = 16, height = 17, limitsize = FALSE, units = 'cm')

# And draw for appendix part

# First legend
g_appendix_legend <- .plot_th_corr_incl_ret(COLFACTORS = 'Mkt.RF', ROWFACTORS = 'HML')+theme(legend.position = 'bottom')+
  scale_colour_manual(labels = c("ARMA-GARCH residual series","Return series"), values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33"))

appendix_legend = gtable_filter(ggplotGrob(g_appendix_legend), "guide-box") 

# Then graphs

appendix_threshold1 <- 
  grid.arrange(grid.arrange(
    .plot_th_corr_incl_ret(COLFACTORS = 'Mkt.RF', ROWFACTORS = 'HML'),
    .plot_th_corr_incl_ret(COLFACTORS = 'Mkt.RF', ROWFACTORS = 'CMA'),
    .plot_th_corr_incl_ret(COLFACTORS = 'Mkt.RF', ROWFACTORS = 'RMW'),
    .plot_th_corr_incl_ret(COLFACTORS = 'Mom', ROWFACTORS = 'HML'),
    .plot_th_corr_incl_ret(COLFACTORS = 'Mom', ROWFACTORS = 'CMA'),
    .plot_th_corr_incl_ret(COLFACTORS = 'Mom', ROWFACTORS = 'RMW'),
    ncol = 2,
    nrow = 3,
    as.table = FALSE
  ),
  appendix_legend,
  nrow = 2,
  heights = c(13,1)
  )
  
ggsave('output/thresholdCorrelations/appendix_threshold_1.png', appendix_threshold1, width = 16, height = 18, limitsize = FALSE, units = 'cm')

appendix_threshold2 <- 
  grid.arrange(
    grid.arrange(
      .plot_th_corr_incl_ret(COLFACTORS = 'SMB', ROWFACTORS = 'HML'),
      .plot_th_corr_incl_ret(COLFACTORS = 'SMB', ROWFACTORS = 'CMA'),
      .plot_th_corr_incl_ret(COLFACTORS = 'SMB', ROWFACTORS = 'RMW'),
      .plot_th_corr_incl_ret(COLFACTORS = 'HML', ROWFACTORS = 'CMA'),
      .plot_th_corr_incl_ret(COLFACTORS = 'HML', ROWFACTORS = 'RMW'),
      .plot_th_corr_incl_ret(COLFACTORS = 'CMA', ROWFACTORS = 'RMW'),
      ncol = 2,
      nrow = 3,
      as.table = FALSE
    ),
    appendix_legend,
    nrow = 2,
    heights = c(13,1)
  )
  
    
ggsave('output/thresholdCorrelations/appendix_threshold_2.png', appendix_threshold2, width = 16, height = 17, limitsize = FALSE, units = 'cm')

# Threshold correlation with scatter to  explain before main results, for MKT-HML pair ------------------------------
# One function for scatter with residuals and one for scatter with returns, axes differ etc...

# Function for graph

.plot_th_corr_scatter <- function(plotdf, plotdf.scatter,
                                  df.labels, COLFACTORS, ROWFACTORS, OUTNAME,
                                  width, height) {
  # Select the column factors for plot this plot
  plotdf <- plotdf %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  plotdf.scatter <- plotdf.scatter %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  df.labels <- df.labels %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  
  # Create the text df, to not print 1000 values, looks bad
  plotdf.scatter.text <- plotdf.scatter %>% select(-value1, -value2) %>% distinct()
  
  # Then do threshold plot
  g <- ggplot(data = plotdf) +
    geom_ribbon(aes(x = qs, ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
                fill = 'grey10',
                alpha = 0.1
    ) +
    geom_line(aes(x = qs, y = value, color = 'Threshold correlation')) +
    geom_abline(aes(slope = 0, intercept = standard_corr,
                    color = 'Correlation'), colour = 'grey20', size = 0.25, linetype = 2, data = df.labels)+
    theme_Publication() +
    scale_colour_Publication() +
    ylab('Correlation') +
    xlab('Quantiles') +
    coord_cartesian(xlim = c(0.10,0.90), ylim = c(-0.5, 1)) +
    theme(legend.position = 'none')+
    scale_x_continuous(labels = scales::percent, breaks = c(0.10, 0.50, 0.90))+
    annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=1)+
    ggtitle("Correlations plot")
  
  # Then do scatter plot(s)
  g_scatter <- ggplot() +
    geom_bin2d(
      data = plotdf.scatter,
      mapping = aes(
        x = value1, y = value2, color = 'GARCH residuals'
      ),
      binwidth = c(0.10, 0.10)
    )+
    # Boxes
    # 10%
    geom_rect(
      data = plotdf.scatter.text,
      mapping = aes(
        xmin = -5, xmax = px10, ymin = -5, ymax = py10
      ),
      fill = 'grey80', alpha = 0.4, linetype = 1, colour = 'black', size = 0.25
    ) +
    # 49%
    geom_rect(
      data = plotdf.scatter.text,
      mapping = aes(
        xmin = -5, xmax = px49, ymin = -5, ymax = py49
      ),
      fill = 'grey80', alpha = 0.5, linetype = 2, colour = 'black', size = 0.25
    )+
    # 51%
    geom_rect(
      data = plotdf.scatter.text,
      mapping = aes(
        xmin = px51, xmax = 5, ymin = py51, ymax = 5
      ),
      fill = 'grey80', alpha = 0.5, linetype = 2, colour = 'black', size = 0.25
    )+
    # 90%
    geom_rect(
      data = plotdf.scatter.text,
      mapping = aes(
        xmin = px90, xmax = 5, ymin = py90, ymax = 5
      ),
      fill = 'grey80', alpha = 0.4, linetype = 1, colour = 'black', size = 0.25
    )+
    # Text
    # 10 %
    geom_text(data = plotdf.scatter.text,
              aes(x = -3.75, y = py10-0.3, label = 'ThCorr 10%'), family = 'Minion Pro', size = 2, parse = F)+
    # 49%
    geom_text(data = plotdf.scatter.text,
              aes(x = -3.75, y = py49-0.3, label = 'ThCorr 49%'), family = 'Minion Pro', size = 2, parse = F)+
    # 51%
    geom_text(data = plotdf.scatter.text,
              aes(x = 3.75, y = py51+0.3, label = 'ThCorr 51%'), family = 'Minion Pro', size = 2, parse = F)+
    # 90%
    geom_text(data = plotdf.scatter.text,
              aes(x = 3.75, y = py90+0.3, label = paste('ThCorr 90%')), family = 'Minion Pro', size = 2, parse = F)+
    
    # Theme options
    theme_Publication() +
    scale_colour_Publication()+
    theme(legend.position = 'none')+
    # Axes
    ylab('HML') +
    xlab('Mkt.RF') +
    coord_cartesian(xlim = c(-5,5), ylim = c(-5, 5))+
    annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=1)+
    ggtitle("Scatter plot")
  
  # Combine the two in grid
  out.graph <- arrangeGrob(g_scatter, g, ncol = 2)
  # Save plot
  
  OUTPATH <- 'output/thresholdCorrelations/threshold_scatter_%s.png'
  ggsave(sprintf(OUTPATH, OUTNAME),
         out.graph, device = 'png', width = width, height = height, units = 'cm'
  )
  
}

.plot_th_corr_scatter_ret <- function(plotdf, plotdf.scatter,
                                      df.labels, COLFACTORS, ROWFACTORS, OUTNAME,
                                      width, height) {
  # Select the column factors for plot this plot
  plotdf <- plotdf %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  plotdf.scatter <- plotdf.scatter %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  df.labels <- df.labels %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  
  # Create the text df, to not print 1000 values, looks bad
  plotdf.scatter.text <- plotdf.scatter %>% select(-value1, -value2) %>% distinct()
  
  # Then do threshold plot
  g <- ggplot(data = plotdf) +
    geom_ribbon(aes(x = qs, ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
                fill = 'grey10',
                alpha = 0.1
    ) +
    geom_line(aes(x = qs, y = value, color = 'Threshold correlation')) +
    geom_abline(aes(slope = 0, intercept = standard_corr,
                    color = 'Correlation'), colour = 'grey20', size = 0.25, linetype = 2, data = df.labels)+
    theme_Publication() +
    scale_colour_Publication() +
    ylab('Correlation') +
    xlab('Quantiles') +
    coord_cartesian(xlim = c(0.10,0.90), ylim = c(-0.5, 1)) +
    theme(legend.position = 'none')+
    scale_x_continuous(labels = scales::percent, breaks = c(0.10, 0.50, 0.90))+
    annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=1)+
    ggtitle("Correlations plot")
  
  # Then do scatter plot(s)
  g_scatter <- ggplot() +
    geom_bin2d(
      data = plotdf.scatter,
      mapping = aes(
        x = value1, y = value2, color = 'GARCH residuals'
      ),
      binwidth = c(0.002, 0.002)
    )+
    # Boxes
    # 10%
    geom_rect(
      data = plotdf.scatter.text,
      mapping = aes(
        xmin = -.1, xmax = px10, ymin = -.1, ymax = py10
      ),
      fill = 'grey80', alpha = 0.4, linetype = 1, colour = 'black', size = 0.25
    ) +
    # 49%
    geom_rect(
      data = plotdf.scatter.text,
      mapping = aes(
        xmin = -.1, xmax = px49, ymin = -.1, ymax = py49
      ),
      fill = 'grey80', alpha = 0.5, linetype = 2, colour = 'black', size = 0.25
    )+
    # 51%
    geom_rect(
      data = plotdf.scatter.text,
      mapping = aes(
        xmin = px51, xmax = .1, ymin = py51, ymax = .1
      ),
      fill = 'grey80', alpha = 0.5, linetype = 2, colour = 'black', size = 0.25
    )+
    # 90%
    geom_rect(
      data = plotdf.scatter.text,
      mapping = aes(
        xmin = px90, xmax = .1, ymin = py90, ymax = .1
      ),
      fill = 'grey80', alpha = 0.4, linetype = 1, colour = 'black', size = 0.25
    )+
    # Text
    # 10 %
    geom_text(data = plotdf.scatter.text,
              aes(x = -.075, y = py10-0.006, label = 'ThCorr 10%'), family = 'Minion Pro', size = 2, parse = F)+
    # 49%
    geom_text(data = plotdf.scatter.text,
              aes(x = -.075, y = py49-0.006, label = 'ThCorr 49%'), family = 'Minion Pro', size = 2, parse = F)+
    # 51%
    geom_text(data = plotdf.scatter.text,
              aes(x = .075, y = py51+0.006, label = 'ThCorr 51%'), family = 'Minion Pro', size = 2, parse = F)+
    # 90%
    geom_text(data = plotdf.scatter.text,
              aes(x = .075, y = py90+0.006, label = paste('ThCorr 90%')), family = 'Minion Pro', size = 2, parse = F)+
    
    # Theme options
    theme_Publication() +
    scale_colour_Publication()+
    scale_x_continuous(labels = scales::percent, breaks = c(-.10, -.05, 0, .05, .10)) +
    scale_y_continuous(labels = scales::percent, breaks = c(-.10, -.05, 0, .05, .10)) +
    theme(legend.position = 'none')+
    # Axes
    ylab('HML') +
    xlab('Mkt.RF') +
    coord_cartesian(xlim = c(-.1,.1), ylim = c(-.1, .1))+
    # Correlation coefficient
    #geom_text(data = df.labels, aes(x = 0, y = -4.5, label = paste("r = ", standard_corr)), family = 'Minion Pro', size = 3, parse = F)+
    # Facet
    annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=1)+
    ggtitle("Scatter plot")
  
  # Combine the two in grid
  out.graph <- arrangeGrob(g_scatter, g, ncol = 2)
  # Save plot
  
  OUTPATH <- 'output/thresholdCorrelations/threshold_scatter_ret_%s.png'
  ggsave(sprintf(OUTPATH, OUTNAME),
         out.graph, device = 'png', width = width, height = height, units = 'cm'
  )
  
}

# Do graph with residuals and returns
.plot_th_corr_scatter(plotdf = plotdf.res, plotdf.scatter = plotdf.scatter.res, df.labels.res,
                      COLFACTORS = c('Mkt.RF'), ROWFACTORS = c('HML'), sprintf('%s_MKT_HML', ID),
                      16, 9)

.plot_th_corr_scatter_ret(plotdf = plotdf.ret, plotdf.scatter = plotdf.scatter.ret, df.labels.ret,
                          COLFACTORS = c('Mkt.RF'), ROWFACTORS = c('HML'), sprintf('%s_MKT_HML', ID),
                          16, 9)

# Threshold correlations simulated - for comparing compula performance to real world ----------------------------------------

# Function for plot

.plot_th_corr_simulated <- function(plotdf = models, ribbondf = models_empirical, df.labels = df.labels.res,
                                    COLFACTORS, ROWFACTORS) {
  # Select the column factors for plot this plot
  plotdf <- plotdf %>% filter(factor1 %in% COLFACTORS, factor2 %in% ROWFACTORS)
  ribbondf <- ribbondf %>% filter(factor1 %in% COLFACTORS, factor2 %in% ROWFACTORS)
  df.labels$factor1 <- df.labels$order
  df.labels$factor2 <- df.labels$order2
  df.labels <- df.labels %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  
  # Then do threshold plot
  g <- ggplot(data = plotdf) +
    geom_line(aes(x = qs, y = coef, color = model)) +
    geom_ribbon(data = ribbondf,
                mapping = aes(x = qs, ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
                fill = 'grey10',
                alpha = 0.1
    ) +
    # geom_abline(aes(slope = 0, intercept = standard_corr,
    #                 color = 'Correlation'), colour = 'grey20', linetype = 2, size = 0.25, data = df.labels)+
    theme_Publication() +
    scale_colour_Publication() +
    ylab('Correlation') +
    xlab('Quantiles') +
    coord_cartesian(xlim = c(0.10,0.90), ylim = c(-0.20, .60)) + 
    theme(legend.position = 'none')+
    scale_x_continuous(labels = scales::percent, breaks = c(0.10, 0.50, 0.90)) +
    annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=1)+
    ggtitle(sprintf("%s - %s", COLFACTORS, ROWFACTORS))
  
  g
  
}

# Consolidate simulated and empirical data

load('data/derived/correlations_threshold_copula.RData')
models_simulated <- models
models_empirical <- plotdf.res %>% 
  select(factor1 = order,
         factor2 = order2,
         qs = qs,
         coef = value,
         lb = lb,
         ub = ub
  ) %>%
  mutate(model = 'empirical')

models <- bind_rows(models, models_empirical) %>% select(-ub, -lb)
models$model <- factor(models$model, levels = c('empirical','norm','std','ghst'))
models$factor1 <- factor(models$factor1, levels = c('Mkt.RF','SMB','Mom','HML','CMA','RMW'))
models$factor2 <- factor(models$factor2, levels = c('Mkt.RF','SMB','Mom','HML','CMA','RMW'))

# Save legend
g_simulated_legend <- .plot_th_corr_simulated(COLFACTORS = 'Mkt.RF', ROWFACTORS = 'HML')+theme(legend.position = 'bottom')+
  scale_colour_manual(labels = c("Empirical distribution","Normal","Symmetric t","Asymmetric t"), values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33"))

simulated_legend = gtable_filter(ggplotGrob(g_simulated_legend), "guide-box") 

# Plot simulated graphs

threshold_simulated1 <- 
  grid.arrange(
    grid.arrange(
      .plot_th_corr_simulated(COLFACTORS = 'Mkt.RF', ROWFACTORS = 'HML'),
      .plot_th_corr_simulated(COLFACTORS = 'Mkt.RF', ROWFACTORS = 'CMA'),
      .plot_th_corr_simulated(COLFACTORS = 'Mkt.RF', ROWFACTORS = 'RMW'),
      .plot_th_corr_simulated(COLFACTORS = 'Mom', ROWFACTORS = 'HML'),
      .plot_th_corr_simulated(COLFACTORS = 'Mom', ROWFACTORS = 'CMA'),
      .plot_th_corr_simulated(COLFACTORS = 'Mom', ROWFACTORS = 'RMW'),
      ncol = 2,
      nrow = 3,
      as.table = FALSE
    ),
    simulated_legend,
    nrow = 2,
    heights = c(13,1)
  )
  
ggsave('output/thresholdCorrelations/threshold_simulated_1.png', threshold_simulated1, width = 16, height = 18, limitsize = FALSE, units = 'cm')

threshold_simulated2 <- 
  grid.arrange(
    grid.arrange(
      .plot_th_corr_simulated(COLFACTORS = 'SMB', ROWFACTORS = 'HML'),
      .plot_th_corr_simulated(COLFACTORS = 'SMB', ROWFACTORS = 'CMA'),
      .plot_th_corr_simulated(COLFACTORS = 'SMB', ROWFACTORS = 'RMW'),
      .plot_th_corr_simulated(COLFACTORS = 'HML', ROWFACTORS = 'CMA'),
      .plot_th_corr_simulated(COLFACTORS = 'HML', ROWFACTORS = 'RMW'),
      .plot_th_corr_simulated(COLFACTORS = 'CMA', ROWFACTORS = 'RMW'),
      ncol = 2,
      nrow = 3,
      as.table = FALSE
    ),
    simulated_legend,
    nrow = 2,
    heights = c(13,1)
  )
ggsave('output/thresholdCorrelations/threshold_simulated_2.png', threshold_simulated2, width = 16, height = 17, limitsize = FALSE, units = 'cm')


#  ------------------------------------------------------------------------


# Rolling correlations ----------------------------------------------------

# Create list to hold correlation sets
rollCorrList.ret = roll_corr(df = df.estim %>% select(-Date), 
                             df.date = df.estim %>% select(Date), 
                             window = 52
)

rollCorrList.res = roll_corr(df = df.stdres %>% select(-Date), 
                             df.date = df.stdres %>% select(Date), 
                             window = 52
)



# Bind to one df for plot
plotdf.ret <- bind_rows(rollCorrList.ret$HML, rollCorrList.ret$RMW, rollCorrList.ret$CMA)
plotdf.res <- bind_rows(rollCorrList.res$HML, rollCorrList.res$RMW, rollCorrList.res$CMA)

# Quick fix to append standard correlations to graphs
df.labels.ret <- .append_standard_corr(plotdf.ret, df.estim)
df.labels.res <- .append_standard_corr(plotdf.res, df.stdres)

# Order the variables
plotdf.ret$order <- factor(plotdf.ret$order, levels = c('Mkt.RF','SMB','Mom','HML','CMA','RMW'))
plotdf.ret$order2 <- factor(plotdf.ret$order2, levels = c('Mkt.RF','SMB','Mom','HML','CMA','RMW'))
plotdf.res$order <- factor(plotdf.res$order, levels = c('Mkt.RF','SMB','Mom','HML','CMA','RMW'))
plotdf.res$order2 <- factor(plotdf.res$order2, levels = c('Mkt.RF','SMB','Mom','HML','CMA','RMW'))

# Roll plot functions ----------------------------------------------

.plot_roll <- function(plot_df = plotdf.res, 
                       df_labels = df.labels.res, 
                       df_stdres = df.stdres, COLFACTORS, ROWFACTORS) {
  # For COLFACTORS ROWFACTORS only
  plot_df <- plot_df %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  df_labels <- df_labels %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  
  g <- ggplot(plot_df, aes(x = Date, y = value)
  ) +
    geom_ribbon(aes(ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
                fill = 'grey10',
                alpha = 0.1
    ) +
    geom_line(aes(color = 'Rolling correlation')) +
    geom_abline(aes(slope = 0, intercept = standard_corr,
                    color = 'Unconditional correlation'), colour = 'grey20', size = 0.25, linetype = 2, data = df_labels)+
    theme_Publication() +
    scale_colour_Publication() +
    ylab('Correlation') +
    xlab('Year') +
    scale_x_date(date_labels = "%y") +
    theme(legend.position = 'none')+
    coord_cartesian(ylim = c(-1, 1), xlim = c(df_stdres$Date[1], df_stdres$Date[length(df_stdres$Date)]))+
    annotate("segment",x=tail(df_stdres$Date,1),xend=head(df_stdres$Date,1),y=Inf,yend=Inf,color="black",lwd=1)+
    ggtitle(sprintf("%s - %s", COLFACTORS, ROWFACTORS))
  
  g
}

.plot_roll_incl_ret <- function(plot_df = plotdf.res, plot_df2 = plotdf.ret, 
                                df_labels = df.labels.res, df_labels2 = df.labels.ret, 
                                df_stdres = df.stdres,
                                COLFACTORS, ROWFACTORS) {
  
  # For COLFACTORS ROWFACTORS only
  plot_df <- plot_df %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  plot_df2 <- plot_df2 %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  df_labels <- df_labels %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  df_labels2 <- df_labels2 %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  
  g <- ggplot(plot_df, aes(x = Date, y = value)
  ) +
    geom_ribbon(aes(ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
                fill = 'grey10',
                alpha = 0.1
    ) +
    geom_ribbon(aes(ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
                fill = 'grey10',
                alpha = 0.1,
                data = plot_df2
    ) +
    geom_line(aes(color = 'Rolling correlation residuals')) +
    geom_line(aes(color = 'Rolling correlations returns'),
              data = plot_df2) +
    geom_abline(aes(slope = 0, intercept = standard_corr,
                    color = 'Unconditional correlation residuals'), colour = 'grey20', size = 0.25, linetype = 2, data = df_labels)+
    geom_abline(aes(slope = 0, intercept = standard_corr,
                    color = 'Unconditional correlation returns'), colour = 'grey20', size = 0.25, linetype = 3, data = df_labels2)+
    theme_Publication() +
    theme(legend.position = 'none')+
    scale_colour_Publication() +
    ylab('Correlation') +
    xlab('Year') +
    scale_x_date(date_labels = "%y") +
    coord_cartesian(ylim = c(-1, 1), xlim = c(df_stdres$Date[1], df_stdres$Date[length(df_stdres$Date)]))+
    annotate("segment",x=tail(df_stdres$Date,1),xend=head(df_stdres$Date,1),y=Inf,yend=Inf,color="black",lwd=1)+
    ggtitle(sprintf("%s - %s", COLFACTORS, ROWFACTORS))
  
  g
}

# Do roll plots and save

rolling1 <- grid.arrange(
  .plot_roll(COLFACTORS = 'Mkt.RF', ROWFACTORS = 'HML'),
  .plot_roll(COLFACTORS = 'Mkt.RF', ROWFACTORS = 'CMA'),
  .plot_roll(COLFACTORS = 'Mkt.RF', ROWFACTORS = 'RMW'),
  .plot_roll(COLFACTORS = 'Mom', ROWFACTORS = 'HML'),
  .plot_roll(COLFACTORS = 'Mom', ROWFACTORS = 'CMA'),
  .plot_roll(COLFACTORS = 'Mom', ROWFACTORS = 'RMW'),
  ncol = 2,
  nrow = 3,
  as.table = FALSE
)
ggsave('output/rollingCorrelations/rolling1.png', rolling1, width = 16, height = 17, limitsize = FALSE, units = 'cm')


rolling2 <- grid.arrange(
  .plot_roll(COLFACTORS = 'SMB', ROWFACTORS = 'HML'),
  .plot_roll(COLFACTORS = 'SMB', ROWFACTORS = 'CMA'),
  .plot_roll(COLFACTORS = 'SMB', ROWFACTORS = 'RMW'),
  .plot_roll(COLFACTORS = 'HML', ROWFACTORS = 'CMA'),
  .plot_roll(COLFACTORS = 'HML', ROWFACTORS = 'RMW'),
  .plot_roll(COLFACTORS = 'CMA', ROWFACTORS = 'RMW'),
  ncol = 2,
  nrow = 3,
  as.table = FALSE
)
ggsave('output/rollingCorrelations/rolling2.png', rolling2, width = 16, height = 17, limitsize = FALSE, units = 'cm')


appendix_rolling1 <- grid.arrange(
  .plot_roll_incl_ret(COLFACTORS = 'Mkt.RF', ROWFACTORS = 'HML'),
  .plot_roll_incl_ret(COLFACTORS = 'Mkt.RF', ROWFACTORS = 'CMA'),
  .plot_roll_incl_ret(COLFACTORS = 'Mkt.RF', ROWFACTORS = 'RMW'),
  .plot_roll_incl_ret(COLFACTORS = 'Mom', ROWFACTORS = 'HML'),
  .plot_roll_incl_ret(COLFACTORS = 'Mom', ROWFACTORS = 'CMA'),
  .plot_roll_incl_ret(COLFACTORS = 'Mom', ROWFACTORS = 'RMW'),
  ncol = 2,
  nrow = 3,
  as.table = FALSE
)
ggsave('output/rollingCorrelations/appendix_rolling1.png', appendix_rolling1, width = 16, height = 17, limitsize = FALSE, units = 'cm')


appendix_rolling2 <- grid.arrange(
  .plot_roll_incl_ret(COLFACTORS = 'SMB', ROWFACTORS = 'HML'),
  .plot_roll_incl_ret(COLFACTORS = 'SMB', ROWFACTORS = 'CMA'),
  .plot_roll_incl_ret(COLFACTORS = 'SMB', ROWFACTORS = 'RMW'),
  .plot_roll_incl_ret(COLFACTORS = 'HML', ROWFACTORS = 'CMA'),
  .plot_roll_incl_ret(COLFACTORS = 'HML', ROWFACTORS = 'RMW'),
  .plot_roll_incl_ret(COLFACTORS = 'CMA', ROWFACTORS = 'RMW'),
  ncol = 2,
  nrow = 3,
  as.table = FALSE
)
ggsave('output/rollingCorrelations/appendix_rolling2.png', appendix_rolling2, width = 16, height = 17, limitsize = FALSE, units = 'cm')
