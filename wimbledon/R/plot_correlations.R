# Threshold related -------------------------------------------------------

#' One function for main graph, one for appendix, one for simulated, one to consolidate
#' them into pages

plot_th_main <- function(pair) {
  # Subset data to pair
  plot_df <- th_corr_stdres %>%
    filter(factor1 == pair[1], factor2 == pair[2])
  
  # Do plot
  g <- ggplot(plot_df, aes(x = qs, y = coef)) +
    # Conf. bound
    geom_ribbon(aes(x = qs, ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
                fill = 'grey10',
                alpha = 0.1,
                show.legend = FALSE
    ) +
    
    # Coef. line
    geom_line(aes(x = qs, y = coef, colour = model)) +
    
    # Unc. coef dash line
    geom_line(aes(x = qs, y = unc_coef, colour = model),
              size = 0.5, linetype = 2
    )+
    
    # Theme options 
    theme_Publication() +
    scale_colour_Publication()+
    theme(legend.position = 'none')+
    annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=1)+
    
    # Axes options
    ylab('Correlation') +
    xlab('Quantiles') +
    coord_cartesian(xlim = c(0.10,0.90), ylim = c(-0.50, .75)) + 
    scale_x_continuous(labels = scales::percent, breaks = c(0.10, 0.50, 0.90)) +
    
    ggtitle(sprintf("%s - %s", pair[1], pair[2]))
  
  g
  
}

plot_th_appendix <- function(pair) {
  
  # Bind return and stdres data and subset to pair
  plot_df <- bind_rows(th_corr_stdres, th_corr_returns) %>%
    filter(factor1 == pair[1], factor2 == pair[2])
  
  # Set model ordering
  plot_df$model <- factor(plot_df$model, levels = MODEL_ORDER)
  
  # Do plot
  g <- ggplot(plot_df) +
    # Conf. bounds
    geom_ribbon(aes(x = qs, ymin = lb, ymax = ub, linetype = NA, fill = 'grey40', colour = model),
                fill = 'grey10',
                alpha = 0.1,
                show.legend = FALSE
    ) +
    
    # Coef. lines
    geom_line(aes(x = qs, y = coef, colour = model)) +
    
    # Unc. coef dash line
    geom_line(aes(x = qs, y = unc_coef, colour = model),
              linetype = 2, size = 0.5
    )+
    
    # Theme options 
    theme_Publication() +
    scale_colour_Publication()+
    theme(legend.position = 'none')+
    annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=1)+
    
    # Axes options
    ylab('Correlation') +
    xlab('Quantiles') +
    coord_cartesian(xlim = c(0.10,0.90), ylim = c(-0.50, .75)) + 
    scale_x_continuous(labels = scales::percent, breaks = c(0.10, 0.50, 0.90)) +
    
    ggtitle(sprintf("%s - %s", pair[1], pair[2]))
  
  g
  
}

plot_th_simulated <- function(pair) {
  
  # Bind return and stdres data and subset to pair
  plot_df <- bind_rows(th_corr_stdres, th_corr_simulated) %>%
    filter(factor1 == pair[1], factor2 == pair[2])
  
  # Set model ordering
  plot_df$model <- factor(plot_df$model, levels = MODEL_ORDER)
  
  # Do plot
  g <- ggplot(plot_df, aes(x = qs, y = coef)) +
    # Conf. bounds
    geom_ribbon(aes(x = qs, ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
                data = plot_df %>% filter(model == 'stdres'),
                fill = 'grey10',
                alpha = 0.1,
                show.legend = FALSE
    ) +
    
    # Coef. lines
    geom_line(aes(x = qs, y = coef, colour = model)) +
    
    # Theme options 
    theme_Publication() +
    scale_colour_Publication()+
    theme(legend.position = 'none')+
    annotate("segment",x=Inf,xend=-Inf,y=Inf,yend=Inf,color="black",lwd=1)+
    
    # Axes options
    ylab('Correlation') +
    xlab('Quantiles') +
    coord_cartesian(xlim = c(0.10,0.90), ylim = c(-0.20, .60)) + 
    scale_x_continuous(labels = scales::percent, breaks = c(0.10, 0.50, 0.90)) +
    
    ggtitle(sprintf("%s - %s", pair[1], pair[2]))
  
  g
  
}


# Rolling related ---------------------------------------------------------

plot_roll_main <- function(pair) {
  # Subset data to pair
  plot_df <- roll_corr_stdres %>%
    filter(factor1 == pair[1], factor2 == pair[2])
  
  # Do plot
  g <- ggplot(plot_df, aes(x = Date, y = coef)) +
    # Conf. bound
    geom_ribbon(aes(x = Date, ymin = lb, ymax = ub, linetype = NA, fill = 'grey40'),
                fill = 'grey10',
                alpha = 0.1,
                show.legend = FALSE
    ) +
    
    # Coef. line
    geom_line(aes(x = Date, y = coef, colour = model)) +
    
    # Unc. coef dash line
    geom_line(aes(x = Date, y = unc_coef, colour = model),
              size = 0.5, linetype = 2
    )+
    
    # Theme options 
    theme_Publication() +
    scale_colour_Publication()+
    theme(legend.position = 'none')+
    annotate("segment",x=tail(plot_df$Date,1),xend=plot_df$Date[1],y=Inf,yend=Inf,color="black",lwd=1)+
    
    # Axes options
    ylab('Correlation') +
    xlab('Year') +
    coord_cartesian(ylim = c(-1, 1), xlim = c(plot_df$Date[1], tail(plot_df$Date,1)))+
    scale_x_date(date_labels = "'%y") +
    
    ggtitle(sprintf("%s - %s", pair[1], pair[2]))
  
  g
  
}

plot_roll_appendix <- function(pair) {
  # Subset data to pair
  plot_df <- bind_rows(roll_corr_stdres, roll_corr_returns) %>%
    filter(factor1 == pair[1], factor2 == pair[2])
  
  # Set model ordering
  plot_df$model <- factor(plot_df$model, levels = MODEL_ORDER)
  
  # Do plot
  g <- ggplot(plot_df, aes(x = Date, y = coef)) +
    # Conf. bounds
    geom_ribbon(aes(x = Date, ymin = lb, ymax = ub, linetype = NA, fill = 'grey40', colour = model),
                fill = 'grey10',
                alpha = 0.1,
                show.legend = FALSE
    ) +
    
    # Coef. lines
    geom_line(aes(x = Date, y = coef, colour = model)) +
    
    # Unc. coef dash line
    # geom_line(aes(x = Date, y = unc_coef, colour = model),
    #           linetype = 2, size = 0.5
    # )+
    
    # Theme options 
    theme_Publication() +
    scale_colour_Publication()+
    theme(legend.position = 'none')+
    annotate("segment",x=tail(plot_df$Date,1),xend=plot_df$Date[1],y=Inf,yend=Inf,color="black",lwd=1)+
    
    # Axes options
    ylab('Correlation') +
    xlab('Year') +
    coord_cartesian(ylim = c(-1, 1), xlim = c(plot_df$Date[1], tail(plot_df$Date,1)))+
    scale_x_date(date_labels = "'%y") +
    
    ggtitle(sprintf("%s - %s", pair[1], pair[2]))
  
  g
  
}

plot_roll_simulated <- function(pair) {
  # Subset data to pair
  plot_df <- bind_rows(roll_corr_stdres, roll_corr_simulated) %>%
    filter(factor1 == pair[1], factor2 == pair[2], model %in% c('stdres','std'))
  
  # Set model ordering
  plot_df$model <- factor(plot_df$model, levels = MODEL_ORDER)
  
  # Do plot
  g <- ggplot(plot_df, aes(x = Date, y = coef)) +
    # Conf. bounds
    geom_ribbon(aes(x = Date, ymin = lb, ymax = ub, linetype = NA, fill = 'grey40', colour = model),
                fill = 'grey10',
                alpha = 0.1,
                show.legend = FALSE
    ) +
    
    # Coef. lines
    geom_line(aes(x = Date, y = coef, colour = model)) +
    
    # Unc. coef dash line
    # geom_line(aes(x = Date, y = unc_coef, colour = model),
    #           linetype = 2, size = 0.5
    # )+
    
    # Theme options 
    theme_Publication() +
    scale_colour_Publication()+
    theme(legend.position = 'none')+
    annotate("segment",x=tail(plot_df$Date,1),xend=plot_df$Date[1],y=Inf,yend=Inf,color="black",lwd=1)+
    
    # Axes options
    ylab('Correlation') +
    xlab('Year') +
    coord_cartesian(ylim = c(-1, 1), xlim = c(plot_df$Date[1], tail(plot_df$Date,1)))+
    scale_x_date(date_labels = "'%y") +
    
    ggtitle(sprintf("%s - %s", pair[1], pair[2]))
  
  g
  
}


# Page structure related --------------------------------------------------

# Function to generate pages
generate_page <- function(order, func, labels = NULL, width = 16, height = 20) {
  # Generate list of the plots
  g_list <- lapply(order, func)
  g <- grid.arrange(g_list[[1]],
                    g_list[[2]],
                    g_list[[3]],
                    g_list[[4]],
                    g_list[[5]],
                    g_list[[6]],
                    nrow = 3,
                    ncol = 2,
                    as.table = FALSE
  )
  # If labels, add the legend
  if(!(is.null(labels))) {
    
    # Get legend
    legend <- g_list[[1]] +
      theme(legend.position = 'bottom') +
      scale_colour_manual(
        labels = labels, 
        values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")
      )
    
    legend <- gtable_filter(ggplotGrob(legend), "guide-box")
    
    # Append and print
    g <- grid.arrange(
      g,
      legend,
      nrow = 2,
      heights = c(height-1, 1)
    )
  }
  
  g
}