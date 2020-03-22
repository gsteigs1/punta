plot_bin_forest <- function(data, x.lab = "Response probability", y.lab = ""){
  fig <- ggplot(data = data) +
    geom_linerange(aes(x = ro, ymin = ci_lo, ymax = ci_up)) +
    geom_text(aes(x = ro, y = -0.05, label = node), hjust = "right") +
    theme_bw() +
    theme(#panel.border = element_blank(), 
      panel.grid = element_blank(),
      #axis.line.x = element_line(),
      axis.text.y = element_blank(), 
      axis.ticks.y = element_blank()) +
    scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(-0.75, 1)) +
    xlab(y.lab) +
    ylab(x.lab) +
    coord_flip()  
  
  return(fig)
}

