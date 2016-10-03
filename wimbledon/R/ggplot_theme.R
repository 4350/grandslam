#' Theme base from
#' https://rpubs.com/Koundy/71792
#' Minion Pro theme imported using
#' extrafont package
#' library(extrafont)
#' import_font('fontpath'), but font needs to be .tty
#' afterwards run fonts() to see its inside :)

theme_Publication <- function(base_size=11, base_family="Minion Pro") {
  library(grid)
  library(ggthemes)
  half_line <- base_size/2
  (theme_foundation(base_size = base_size, base_family = base_family)
  + theme(plot.title = element_text(family = "Minion Pro", 
                                    face = "bold",
                                    size = rel(1.2), 
                                    hjust = 0.5),
          text = element_text(family = "Minion Pro", 
                              face = "plain",
                              colour = "black", 
                              size = base_size,
                              lineheight = 0.9,  
                              hjust = 0.5,
                              vjust = 0.5, 
                              angle = 0, 
                              margin = margin(), 
                              debug = FALSE),
          panel.background = element_rect(colour = NA, fill = 'grey92'),
          plot.background = element_rect(colour = NA),
          panel.border = element_rect(colour = NA),
          axis.title = element_text(face = "bold",size = rel(1)),
          axis.title.x = element_text(margin = margin(t = 0.8 * half_line,
                                                      b = 0.8 * half_line/2)),
          axis.title.y = element_text(angle = 90,
                                      margin = margin(r = 0.8 * half_line,
                                                      l = 0.8 * half_line/2)),
          axis.text = element_text(size = rel(0.8), colour = "grey30"), 
          axis.line = element_blank(),
          axis.ticks = element_line(colour = "grey20"),
          axis.text.x = element_text(margin = margin(t = 0.8*half_line/2), 
                                     vjust = 1), 
          axis.text.y = element_text(margin = margin(r = 0.8*half_line/2),
                                     hjust = 1),
          panel.grid.major = element_line(colour = "white"), 
          panel.grid.minor = element_blank(), 
          legend.key = element_rect(colour = NA),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.justification = "center", 
          legend.key.size= unit(0.2, "cm"),
          legend.margin = unit(0, "cm"),
          legend.title = element_blank(), 
          legend.text = element_text(size = rel(0.8)),
          plot.margin=unit(c(10,5,5,5),"mm"),
          strip.background=element_rect(colour = NA, fill="grey85"),
          strip.text = element_text(face = "bold", colour = "grey10", size = rel(0.8)),
          strip.text.x = element_text(margin = margin(t = half_line,
                                                      b = half_line)), 
          strip.text.y = element_text(angle = -90, 
                                      margin = margin(l = half_line, 
                                                      r = half_line)),
          strip.switch.pad.grid = unit(0.1, "cm"),
          strip.switch.pad.wrap = unit(0.1, "cm")
  ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}
