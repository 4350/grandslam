#' Functionality related to threshold correlations
#'

#' Threshold correlations
#'
#' Takes a data frame, and returns estimates of threshold corr, lb, ub
#' using the factor names given
#'
#' 
#' @param pair Two strings in character vector
#' @param qs Sequence of quantile values
#' @param df Data frame of return series
#' 
#' @return result GGplot-tidy data frame w coef, lb, ub, unc_coef, qs, factors
#' @export
correlations_threshold <- function(pair, qs, df) {
  
  factor1 <- pair[[1]]
  factor2 <- pair[[2]]
  
  result <- lapply(
    qs, function(q, df, factor1, factor2) {
      # Extract vectors
      x = df[, factor1]
      y = df[, factor2]
      # Find empirical percentiles
      quantile.x <- quantile(x, q)
      quantile.y <- quantile(y, q)
      
      # Check which of x and y to include in correlation estimation
      if (q < 0.50) {
        inc <- x < quantile.x & y < quantile.y
      }
      else {
        inc <- x >= quantile.x & y >= quantile.y
      }
      
      # Estimate correlation
      result <- cor.test(x[inc], y[inc])
      
      # Return coef, lb, ub, unc_coef
      result <- data.frame(coef = result$estimate[[1]],
                           lb = result$conf.int[1],
                           ub = result$conf.int[2],
                           unc_coef = cor(x,y)[[1]])
      
      
      return(result)
    },
    df = df, factor1 = factor1, factor2 = factor2
  )
  
  result <- data.frame(
    bind_rows(result),
    qs = qs,
    factor1 = factor1,
    factor2 = factor2
  )
  
  return(result)
}



#' Get threshold correlation data frame for plotting for a
#' set of factor pairs (in pairs list)
#' @param pairs list of vectors with the two factor strings
#' @param qs sequence of quantiles of interest
#' @param df to run the threshold correlation on
#' @return results rows bound together for all pairs, tidy data
#' @export
threshold_correlations_pairs <- function(pairs, qs, df) {
  # List apply pairs in threshold correlation function
  results <- lapply(pairs, function(pair, qs, df) {
    correlations_threshold(pair, qs, df)
  },
  qs = qs, 
  df = df
  )
  
  bind_rows(results)
  
}

#' Get threshold correlation data frame for plotting for a
#' number of copula models. Loads stdresid data accordingly
#' @param dynamic 'constant' or 'dynamic'
#' @param models list of vectors with the model names norm, std, ghst
#' @param pairs list of vectors with the two factor strings
#' @param qs sequence of quantiles of interest
#' @param df to run the threshold correlation on
#' @return results rows bound together for all model&pairs, tidy data
#' @export
threshold_correlations_models <- function(dynamic, models, pairs, qs) {
  # List apply models (each for pairs) in threshold correlation function
  results <- lapply(models, function(model, dynamic, pairs, qs) {
    load(sprintf('data/derived/stdresid/full_%s_%s.RData', dynamic, model))
    
    threshold_correlations_pairs(pairs, qs, stdresid)
  },
  dynamic = dynamic,
  pairs = pairs,
  qs = qs
  )
  
  names(results) <- models
  
  bind_rows(results, .id = 'model')
}

#' Makes the explanatory plot, comparing standardized residuals to the scatter
#' to better understand the threshold function
#' 
#' @param pair two factor strings in vector
#' @return exports graph
#' @export
plot_th_explain <- function(pair, width = 16, height = 9) {
  
  # Subset data to pair
  plot_df_th <- th_corr_stdres %>% 
    filter(factor1 == pair[1], factor2 == pair[2])
  plot_df_scatter <- data.frame(
    x = df.stdres[,pair[1]],
    y = df.stdres[,pair[2]]
  )
  
  # Do threshold plot
  g_th <- plot_th_main(pair) +
    ggtitle('Correlations plot')
  
  # Then do scatter plot
  g_scatter <- ggplot() +
    geom_bin2d(
      data = plot_df_scatter,
      mapping = aes(x = x, y = y),
      binwidth = c(0.10, 0.10)
    )
  
  # Get quantiles to add boxes
  px10 = quantile(plot_df_scatter$x, probs = .10)
  py10 = quantile(plot_df_scatter$y, probs = .10)
  px49 = quantile(plot_df_scatter$x, probs = .49)
  py49 = quantile(plot_df_scatter$y, probs = .49)
  px51 = quantile(plot_df_scatter$x, probs = .51)
  py51 = quantile(plot_df_scatter$y, probs = .51)
  px90 = quantile(plot_df_scatter$x, probs = .90)
  py90 = quantile(plot_df_scatter$y, probs = .90)
  
  # Add boxes to scatter
  g_scatter <- g_scatter +
    # 49%
    geom_rect(
      mapping = aes(
        xmin = -5, xmax = py49, ymin = -5, ymax = py49
      ),
      fill = 'grey80', alpha = 0.5, linetype = 2, colour = 'black', size = 0.25
    )+
    # 51%
    geom_rect(
      mapping = aes(
        xmin = px51, xmax = 5, ymin = py51, ymax = 5
      ),
      fill = 'grey80', alpha = 0.5, linetype = 2, colour = 'black', size = 0.25
    )+
    # 10%
    geom_rect(
      mapping = aes(
        xmin = -5, xmax = px10, ymin = -5, ymax = py10
      ),
      fill = 'grey80', alpha = 0.4, linetype = 1, colour = 'black', size = 0.25
    ) +
    # 90%
    geom_rect(
      mapping = aes(
        xmin = px90, xmax = 5, ymin = py90, ymax = 5
      ),
      fill = 'grey80', alpha = 0.4, linetype = 1, colour = 'black', size = 0.25
    )
  
  # Add text to scatter
  g_scatter <- g_scatter +
    # 10 %
    geom_text(aes(x = -3.75, y = py10-0.3, label = 'ThCorr 10%'), family = 'Minion Pro', size = 2, parse = F)+
    # 49%
    geom_text(aes(x = -3.75, y = py49-0.3, label = 'ThCorr 49%'), family = 'Minion Pro', size = 2, parse = F)+
    # 51%
    geom_text(aes(x = 3.75, y = py51+0.3, label = 'ThCorr 51%'), family = 'Minion Pro', size = 2, parse = F)+
    # 90%
    geom_text(aes(x = 3.75, y = py90+0.3, label = 'ThCorr 90%'), family = 'Minion Pro', size = 2, parse = F)
  
  # Do theme and other options
  g_scatter <- g_scatter +
    theme_Publication() +
    scale_colour_Publication()+
    theme(legend.position = 'none')+
    # Axes
    xlab('Mkt.RF') +
    ylab('HML') +
    coord_cartesian(xlim = c(-5,5), ylim = c(-5, 5))+
    annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=1)+
    ggtitle("Scatter plot")
  
  # Combine the two in grid
  out.graph <- arrangeGrob(g_scatter, g_th, ncol = 2)
  
  # Save plot
  ggsave('output/thresholdCorrelations/threshold_explain_res.png',
         out.graph, width = width, height = height, units = 'cm'
  )
}