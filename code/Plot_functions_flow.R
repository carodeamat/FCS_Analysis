# Plot_functions_flow.R
# utilities file with functions to generate plots in flow cytometry analysis.

#Plot clusters in reduced dimension

plot.feature <- function(data,
                          x = "UMAP_1",
                          y = "UMAP_2",
                          color.by = "Clusters_pheno",
                          split.by = NULL,
                          order = NULL,
                          size = 1,
                          alpha = 1,
                          show.cluser.id = TRUE,
                          show.cluser.id.size = 4,
                          colors = "default",
                          main = "UMAP clusters",
                          aspect.ratio = 1,
                          x.range = "default",
                          y.range = "default",
                          plot.theme = theme_bw()) {
  
  
  if (length(color.by) > 1) {
    warning(Sys.time(), " color.by has more than one elements. Only the first one will be used")
    color.by <- color.by[1]
  }
  
  plot.x = plot.y =NULL
  
  plot.data <- data.frame(plot.x = data[,x],
                          plot.y = data[,y],
                          plot.color.by = as.factor(data[,color.by]))
  
  if (!is.null(split.by)){
    plot.data[,split.by] <- data[,split.by]
  }
  
  if (is.null(order)) {
    plot.data$plot.color.by <- factor(plot.data$plot.color.by)
  } else {
    plot.data$plot.color.by <- factor(as.character(plot.data$plot.color.by), levels = order)
  }
  
  
  # plot
  gg <- ggplot(plot.data) + geom_point(aes(x=plot.x, y=plot.y, color = plot.color.by), size = size, alpha = alpha)
  gg <- gg + plot.theme
  gg <- gg + labs(x = x, y = y, title = paste0(main))
  gg <- gg + labs(color = color.by)
  gg <- gg + theme(aspect.ratio = aspect.ratio)
  
  if (show.cluser.id) {
    pos <- aggregate(  plot.data[, seq_len(2)], list( pos = plot.data$plot.color.by ), mean)
    
    for ( i in seq_along(pos$pos)) {
      gg <- gg + annotate(geom="text", x = pos$plot.x[i], y = pos$plot.y[i], 
                          label = pos$pos[i],
                          size = show.cluser.id.size)
    }
  }
  
  if (all(colors!="default")&is.vector(colors)&length(colors)==length(levels(plot.data$plot.color.by))){
    gg <- gg +  scale_colour_manual(values = colors)
  }
  
  if (all(x.range!="default")&is.vector(x.range)&is.numeric(x.range)&length(x.range)==2){
    gg <- gg +  scale_x_continuous(limits = x.range)
  }else{
    spc <- 0.1*(max(plot.data$plot.x) - min(plot.data$plot.x))
    gg <- gg +  scale_x_continuous(limits = c(min(plot.data$plot.x)-spc, max(plot.data$plot.x)+spc))
  }
  
  if (all(y.range!="default")&is.vector(y.range)&is.numeric(y.range)&length(y.range)==2){
    gg <- gg +  scale_y_continuous(limits = y.range)
  }else{
    spc <- 0.1*(max(plot.data$plot.y) - min(plot.data$plot.y))
    gg <- gg +  scale_y_continuous(limits = c(min(plot.data$plot.y)-spc, max(plot.data$plot.y)+spc))
  }
  
  if (!is.null(split.by)){
    gg <- gg +  facet_grid(. ~ .data[[split.by]])
  }
  
  return(gg)
  
}


#Plot expression in reduced dimension

plot.expr <- function(feat_data,
                      expr_matrix,
                      x = "UMAP_1",
                      y = "UMAP_2",
                      color.by = "CD4",
                      split.by = NULL,
                      size = 1,
                      alpha = 1,
                      main = "UMAP expression",
                      colors = viridis(100),
                      limits = NULL,
                      aspect.ratio = 1,
                      x.range = "default",
                      y.range = "default",
                      plot.theme = theme_bw()) {
  
  
  if (length(color.by) > 1) {
    warning(Sys.time(), " color.by has more than one elements. Only the first one will be used")
    color.by <- color.by[1]
  }
  
  
  plot.x = plot.y =NULL
  
  plot.data <- data.frame(plot.x = feat_data[,x],
                          plot.y = feat_data[,y],
                          plot.color.by = expr_matrix[,color.by])
  
  if (!is.null(split.by)){
    plot.data[,split.by] <- feat_data[,split.by]
  }
  
  # plot
  gg <- ggplot(plot.data) + geom_point(aes(x=plot.x, y=plot.y, color = plot.color.by), size = size, alpha = alpha)
  gg <- gg + plot.theme
  gg <- gg + labs(x = x, y = y, title = paste0(main))
  gg <- gg + labs(color = color.by) 
  gg <- gg + theme(aspect.ratio=aspect.ratio)
  
  if (!is.null(limits)){
    gg <- gg + scale_colour_gradientn(colors = colors, limits = limits)
  } else{
    gg <- gg + scale_colour_gradientn(colors = colors)
  }
  
  if (all(x.range!="default")&is.vector(x.range)&is.numeric(x.range)&length(x.range)==2){
    gg <- gg +  scale_x_continuous(limits = x.range)
  }else{
    spc <- 0.1*(max(plot.data$plot.x) - min(plot.data$plot.x))
    gg <- gg +  scale_x_continuous(limits = c(min(plot.data$plot.x)-spc, max(plot.data$plot.x)+spc))
  }
  
  if (all(y.range!="default")&is.vector(y.range)&is.numeric(y.range)&length(y.range)==2){
    gg <- gg +  scale_y_continuous(limits = y.range)
  }else{
    spc <- 0.1*(max(plot.data$plot.y) - min(plot.data$plot.y))
    gg <- gg +  scale_y_continuous(limits = c(min(plot.data$plot.y)-spc, max(plot.data$plot.y)+spc))
  }
  
  if (!is.null(split.by)){
    gg <- gg +  facet_grid(. ~ .feat_data[[split.by]])
  }
  
  return(gg)
  
}

# Overlaid density plots

dens.plot <- function(data,
                      x = "CD4",
                      color.by,
                      linew = 1.5,
                      alpha = 0.6,
                      colors = "default",
                      xlabel = "marker expression",
                      aspect.ratio = 1.5) {
  
  plot.data <- data.frame(plot.x = data[,x],
                          plot.color.by = data[,color.by])
  
  gg <- ggplot(plot.data, aes(x=plot.x, fill=plot.color.by, y=after_stat(scaled))) + 
    geom_density(alpha=alpha, size = linew) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = expand_scale(mult = c(0, 0.05), 
                                             add = c(0, 0))) +
    xlab(xlabel) +
    theme(aspect.ratio=aspect.ratio)
  #theme_classic()
  if (all(colors!="default")&is.vector(colors)&length(colors)==length(unique(plot.data$plot.color.by))){
    gg <- gg + scale_fill_manual(color.by, values = colors)
  }
  
}
