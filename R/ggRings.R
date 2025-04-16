#' Circular heatmap
#' 
#' generates a radial heatmap (circular rings)
#' expression data using ggplot2. It displays multiple layers (rings) representing
#' different variables.
#' 
#' @param data A data frame containing gene expression data. Each column represents
#' a different variable
#' @param ring_cols A vector of column names from the `data` frame to be plotted as rings.
#' @param annotation_col The column name for annotations (e.g., SNP IDs or gene names).
#' @param ring_colors A vector of colors to use for each ring.
#' @param ring_width The width of each ring.
#' @param title The title for the plot.
#' @param annotation_distance The distance from the center for placing annotations.
#' @param add_ticks Logical indicating whether to add ticks on the outer ring.
#' @param tick_length The length of the ticks.
#' @param fill_low The color used for the lowest values in the rings.
#' @import ggplot2
#' @import dplyr
#' @import ggnewscale
#' @import ggforce
#' @import RColorBrewer
#' @examples
#' data(geneOfInterest)
#' plot_circular_rings(
#'   data = geneOfInterest,
#'   ring_cols = c("ATAC", "R2", "H3K27ac"),
#'   annotation_col = "rs",
#'   ring_colors = c("darkred", "darkblue", "seagreen", "yellow"),
#'   annotation_distance = 1,
#'   add_ticks = TRUE,
#'   tick_length = 0.1,
#'   title = NULL,
#'   fill_low = "white"
#' )
#' @return A `ggplot` object representing the radial heatmap.
#' @export

plot_circular_rings <- function(data, 
                                ring_cols, 
                                annotation_col = "rs",     
                                ring_colors = NULL, 
                                ring_width = 0.5, 
                                title = "title",
                                annotation_distance = 1,  
                                add_ticks = TRUE,           
                                tick_length = 0.1,
                                fill_low = "white") {      
  if (length(fill_low) == 1) {
    fill_low <- rep(fill_low, length(ring_cols))
  } else if (length(fill_low) != length(ring_cols)) {
    stop("Length of fill_low must be 1 or equal to the number of ring_cols.")
  }
  n <- nrow(data)
  data <- data %>%
    mutate(
      rs = as.character(rs),  
      start = 2 * pi * (row_number() - 1) / n,
      end   = 2 * pi * row_number() / n
    )
  if (is.null(ring_colors)) {
    ring_colors <- RColorBrewer::brewer.pal(length(ring_cols), "Set2")
  }
  p <- ggplot()
  for (i in seq_along(ring_cols)) {
    r0 <- 1 + (i - 1) * (ring_width + 0.1)  
    r  <- r0 + ring_width                   
    col <- ring_cols[i]
    high_color <- ring_colors[i]
    low_color <- fill_low[i]
    
    data_layer <- data %>% 
      mutate(r0 = r0, r = r)
    
    p <- p +
      geom_arc_bar(data = data_layer,
                   aes(x0 = 0, y0 = 0, r0 = r0, r = r, 
                       start = start, end = end,
                       fill = !!sym(col))) +
      scale_fill_gradient(
        low = low_color, 
        high = high_color, 
        name = col,
        guide = guide_colorbar(ticks.colour = "black", frame.colour = "black")
      )
    if (i < length(ring_cols)) {
      p <- p + new_scale_fill()
    }
  }
  outer_ring_outer <- 1 + (length(ring_cols) - 1) * (ring_width + 0.1) + ring_width
  annotation_radius <- outer_ring_outer + annotation_distance
  angles <- seq(0, 2 * pi, length.out = n + 1)[-1]   
  annotations_df <- data.frame(
    annotation = data[[annotation_col]],
    angle = angles,
    x = annotation_radius * cos(angles),
    y = annotation_radius * sin(angles)
  )
  p <- p +
    geom_text(data = annotations_df,
              aes(x = x, y = y, label = annotation),
              size = 3, fontface = "bold", color = "black",
              angle = ifelse(annotations_df$angle > pi/2 & annotations_df$angle < 3*pi/2, 
                             annotations_df$angle * 180/pi + 180, 
                             annotations_df$angle * 180/pi),
              hjust = ifelse(annotations_df$angle > pi/2 & annotations_df$angle < 3*pi/2, 1, 0))
  
  if (add_ticks) {
    tick_angles <- seq(0, 2 * pi, length.out = n + 1)
    ticks_df <- data.frame(
      angle = tick_angles,
      x_start = outer_ring_outer * cos(tick_angles),
      y_start = outer_ring_outer * sin(tick_angles),
      x_end = (outer_ring_outer + tick_length) * cos(tick_angles),
      y_end = (outer_ring_outer + tick_length) * sin(tick_angles)
    )
    
    p <- p +
      geom_segment(data = ticks_df,
                   aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                   color = "black", size = 0.5)
  }
  
  p <- p +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "right") +
    labs(title = title)
  
  if (!is.null(title)) {
    p <- p + labs(title = title)
  }
  
  return(p)
}

