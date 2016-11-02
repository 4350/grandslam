#' Estimates threshold and rolling corrs.
#' Currently implementing daily return series. Has previously run weekly
#' series for threshold correlation. The thCorrelationList is saved in derived
#' data for further use in simulated threshold graphs.
#' 


# Libraries ---------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(extrafont)
library(zoo)
library(scales)
library(devtools)
library(mvtnorm)

# Reset workspace and load return data
rm(list = ls())
load('data/derived/weekly-estim.RData')
ID = 'model_GARCH_chosen_stdres.RData'
load(sprintf('data/derived/garch/%s', ID))


# Load functions
load_all('wimbledon')


# Function ----------------------------------------------------------------

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

# Do threshold plot and save ----------------------------------------------

.plot_th_corr <- function(plotdf,
                          df.labels, COLFACTORS, ROWFACTORS, OUTNAME,
                          width, height) {
  # Select the column factors for plot this plot
  plotdf <- plotdf %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  df.labels <- df.labels %>% filter(order %in% COLFACTORS, order2 %in% ROWFACTORS)
  
  # Then do threshold plot
  g <- ggplot(data = plotdf) +
    geom_ribbon(aes(x = qs, ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
                fill = 'grey10',
                alpha = 0.1
    ) +
    geom_line(aes(x = qs, y = value, color = 'Threshold correlation')) +
    geom_abline(aes(slope = 0, intercept = standard_corr,
                    color = 'Correlation'), linetype = 2, data = df.labels)+
    theme_Publication() +
    scale_colour_Publication() +
    ylab('Correlation') +
    xlab('Quantiles') +
    coord_cartesian(xlim = c(0.10,0.90), ylim = c(-0.5, 1)) + 
    theme(legend.position = 'none')+
    scale_x_continuous(labels = scales::percent, breaks = c(0.10, 0.50, 0.90)) +
    facet_grid(order2 ~ order, switch = 'y')
  
  OUTPATH <- 'output/thresholdCorrelations/threshold_%s.png'
  ggsave(sprintf(OUTPATH, OUTNAME),
    g, device = 'png', width = width, height = height, units = 'cm'
    )

}

# Do the plots ------------------------------------------------------------

# Quick fix to append standard correlations to graphs
df.labels <- .append_standard_corr(plotdf.ret, df.estim)

.plot_th_corr(plotdf = plotdf.res, df.labels,
              COLFACTORS = c('Mkt.RF','Mom'), ROWFACTORS = c('HML', 'CMA','RMW'), sprintf('%s_Page1', ID),
              14, 16)
.plot_th_corr(plotdf = plotdf.res, df.labels,
              COLFACTORS = c('RMW','CMA'), ROWFACTORS = c('HML','CMA'), sprintf('%s_Page2', ID),
              14, 12)


# Threshold correlation explain before results, for MKT-HML pair ------------------------------

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
                    color = 'Correlation'), linetype = 2, data = df.labels)+
    theme_Publication() +
    scale_colour_Publication() +
    ylab('Correlation') +
    xlab('Quantiles') +
    coord_cartesian(xlim = c(0.10,0.90), ylim = c(-0.5, 1)) +
    theme(legend.position = 'none')+
    scale_x_continuous(labels = scales::percent, breaks = c(0.10, 0.50, 0.90)) +
    facet_grid(order2 ~ order, switch = 'y')
   
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
    ylab('Standardized residuals (column factor)') +
    xlab('Standardized residuals (row factor)') +
    coord_cartesian(xlim = c(-5,5), ylim = c(-5, 5)) +
    # Correlation coefficient
    #geom_text(data = df.labels, aes(x = 0, y = -4.5, label = paste("r = ", standard_corr)), family = 'Minion Pro', size = 3, parse = F)+
    # Facet
    facet_grid(order2 ~ order, switch = 'y')

  # Combine the two in grid
  out.graph <- arrangeGrob(g_scatter, g, ncol = 2)
  # Save plot

  OUTPATH <- 'output/thresholdCorrelations/threshold_scatter_%s.png'
  ggsave(sprintf(OUTPATH, OUTNAME),
         out.graph, device = 'png', width = width, height = height, units = 'cm'
  )
  
}

# Do graph
.plot_th_corr_scatter(plotdf = plotdf.res, plotdf.scatter = plotdf.scatter.res, df.labels,
              COLFACTORS = c('Mkt.RF'), ROWFACTORS = c('HML'), sprintf('%s_MKT_HML', ID),
              14, 7)


# Threshold correlation fake data scatter for method ------------------------------

# Generate data
sigma = diag(2)
sigma[1,2] = -0.25
sigma[2,1] = -0.25
rand_data <- data.frame(rmvnorm(n = 1000, mean = c(0,0), sigma = sigma))
rand_data$cor = -0.25
colnames(rand_data) = c('y1','y2')

# Plot data
g <- ggplot(rand_data, aes(x = y1, y = y2)
) +
  geom_point() +
  theme_Publication() +
  scale_colour_Publication() +
  ylab('Factor 1 return') +
  xlab('Factor 2 return') +
  coord_cartesian(xlim = c(-4,4), ylim = c(-4, 4)) + 
  annotate("rect", xmin = -4, xmax = 0, ymin = -4, ymax = 0, alpha = 0.4, fill = 'grey80', linetype = 1, color = 'black')+
  annotate("text", x = -3.85, y = -0.2, label = 'A', family = 'Minion Pro')+
  annotate("rect", xmin = -4, xmax = -1, ymin = -4, ymax = -1, alpha = 0.5, fill = 'grey80', linetype = 2, color = 'black')+
  annotate("text", x = -3.85, y = -1.2, label = 'B', family = 'Minion Pro')+
  annotate("rect", xmin = -4, xmax = -2, ymin = -4, ymax = -2, alpha = 0.6, fill = 'grey80', linetype = 3, color = 'black')+
  annotate("text", x = -3.85, y = -2.2, label = 'C', family = 'Minion Pro')+
  annotate("rect", xmin = 0, xmax = 4, ymin = 0, ymax = 4, alpha = 0.4, fill = 'grey80', linetype = 1, color = 'black')+
  annotate("text", x = 3.85, y = 0.2, label = 'D', family = 'Minion Pro')+
  annotate("rect", xmin = 1, xmax = 4, ymin = 1, ymax = 4, alpha = 0.5, fill = 'grey80', linetype = 2, color = 'black')+
  annotate("text", x = 3.85, y = 1.2, label = 'E', family = 'Minion Pro')+
  annotate("rect", xmin = 2, xmax = 4, ymin = 2, ymax = 4, alpha = 0.6, fill = 'grey80', linetype = 3, color = 'black')+
  annotate("text", x = 3.85, y = 2.2, label = 'F', family = 'Minion Pro')

# Save plot
ggsave('output/thresholdCorrelations/threshold_explain.png', g, device = 'png', width = 14, height = 12, units = 'cm')
              

# Threshold correlations simulated ----------------------------------------

# Function for plot


.plot_th_corr_simulated <- function(plotdf, ribbondf,
                          COLFACTORS, ROWFACTORS, OUTNAME,
                          width, height) {
  # Select the column factors for plot this plot
  plotdf <- plotdf %>% filter(factor1 %in% COLFACTORS, factor2 %in% ROWFACTORS)
  ribbondf <- ribbondf %>% filter(factor1 %in% COLFACTORS, factor2 %in% ROWFACTORS)
  
  # Then do threshold plot
  g <- ggplot(data = plotdf) +
    geom_ribbon(data = ribbondf,
                mapping = aes(x = qs, ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
                fill = 'grey10',
                alpha = 0.1
    ) +
    geom_line(aes(x = qs, y = coef, color = model)) +
    theme_Publication() +
    scale_colour_Publication() +
    ylab('Correlation') +
    xlab('Quantiles') +
    coord_cartesian(xlim = c(0.10,0.90), ylim = c(-0.5, 1)) + 
    theme(legend.position = 'none')+
    scale_x_continuous(labels = scales::percent, breaks = c(0.10, 0.50, 0.90)) +
    facet_grid(factor2 ~ factor1, switch = 'y')
  
  OUTPATH <- 'output/thresholdCorrelations/threshold_simulated_%s.png'
  ggsave(sprintf(OUTPATH, OUTNAME),
         g, device = 'png', width = width, height = height, units = 'cm'
  )
  
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


# Plot

.plot_th_corr_simulated(plotdf = models, ribbondf = models_empirical,
              COLFACTORS = c('Mkt.RF','Mom'), ROWFACTORS = c('HML', 'CMA','RMW'), sprintf('%s_Page1', ID),
              14, 16)
.plot_th_corr_simulated(plotdf = models, ribbondf = models_empirical,
              COLFACTORS = c('RMW','CMA'), ROWFACTORS = c('HML','CMA'), sprintf('%s_Page2', ID),
              14, 12)


# Rolling correlations ----------------------------------------------------

# Create list to hold correlation sets
rollCorrList.ret = roll_corr(df = df.estim %>% select(-Date), 
                                   df.date = df.estim %>% select(Date), 
                                   window = 45
                                   )


# Bind to one df for plot
plotdf.ret <- bind_rows(rollCorrList.ret$HML, rollCorrList.ret$RMW, rollCorrList.ret$CMA)

# Quick fix to append standard correlations to graphs
df.labels <- .append_standard_corr(plotdf.ret, df.estim)

# Do roll plot and save ----------------------------------------------
g <- ggplot(plotdf.ret, aes(x = Date, y = value)
) +
  geom_ribbon(aes(ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
              fill = 'grey10',
              alpha = 0.1
  ) +
  geom_line(aes(color = 'Return series')) +
  theme_Publication() +
  scale_colour_Publication() +
  ylab('Correlation') +
  xlab('Year') +
  scale_x_date(date_labels = "%y") +
  coord_cartesian(ylim = c(-1, 1), xlim = c(df.estim$Date[1], df.estim$Date[length(df.estim$Date)])) +
  #annotate("rect", xmin = as.Date('1986-01-01'), xmax = as.Date('1994-01-01'), ymin = -0.95, ymax = -0.5, alpha = 0.8, fill = 'grey80')+
  #geom_text(data = df.labels, aes(x = as.Date('1990-01-01'), y = -0.725, label = paste('r = ',standard_corr)), family = 'Minion Pro', size = 3, parse = FALSE)+
  geom_text(data = df.labels, aes(x = as.Date('2010-01-01'), y = -0.90, label = paste('r = ',standard_corr)), family = 'Minion Pro', size = 3, parse = FALSE)+
  facet_grid(order ~ order2)
  #ggtitle('Rolling 45-week correlations (95% confidence bounds)')#+
  #theme(axis.text = element_text(size = rel(0.6), colour = "grey30")) 

ggsave('output/rollingCorrelations/rolling45.png', g, device = 'png', width = 14, height = 18, units = 'cm', limitsize = F)


# Do roll plot with NBER dummy --------------------------------------------
load('data/derived/usrec-weekly.RData')

g <- ggplot(plotdf.ret) +
  geom_ribbon(aes(x = Date, ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
              fill = 'grey10',
              alpha = 0.1
  ) +
  geom_line(aes(x = Date, y = value, color = 'Return series')) +
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
  ylab('Correlation') +
  xlab('Year') +
  scale_x_date(date_labels = "%y") +
  coord_cartesian(ylim = c(-1, 1), xlim = c(df.estim$Date[1], df.estim$Date[length(df.estim$Date)])) +
  #annotate("rect", xmin = as.Date('1986-01-01'), xmax = as.Date('1994-01-01'), ymin = -0.95, ymax = -0.5, alpha = 0.8, fill = 'grey80')+
  #geom_text(data = df.labels, aes(x = as.Date('1990-01-01'), y = -0.725, label = paste('r = ',standard_corr)), family = 'Minion Pro', size = 3, parse = FALSE)+
  geom_text(data = df.labels, aes(x = as.Date('2010-01-01'), y = -0.90, label = paste('r = ',standard_corr)), family = 'Minion Pro', size = 3, parse = FALSE)+
  facet_grid(order ~ order2)
#ggtitle('Rolling 45-week correlations (95% confidence bounds)')#+
#theme(axis.text = element_text(size = rel(0.6), colour = "grey30")) 

ggsave('output/rollingCorrelations/rolling45NBER.png', g, device = 'png', width = 14, height = 18, units = 'cm', limitsize = F)