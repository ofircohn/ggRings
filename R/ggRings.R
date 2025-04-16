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
                                annotation_distance = 0.3,
                                add_ticks = TRUE,           
                                tick_length = 0.1,
                                fill_low = "white",
                                legend_position = "right",
                                legend_margin = 0.15,       
                                text_size = 3) {      
  if (length(fill_low) == 1) {
    fill_low <- rep(fill_low, length(ring_cols))
  } else if (length(fill_low) != length(ring_cols)) {
    stop("Length of fill_low must be 1 or equal to the number of ring_cols.")
  }
  
  n <- nrow(data)
  
  data <- data %>%
    mutate(
      rs = as.character(rs),  
      id = 1:n,
      segment_angle = 2 * pi / n,
      start = 2 * pi * (row_number() - 1) / n,
      end = 2 * pi * row_number() / n,
      mid = (start + end) / 2  
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
  
  if (add_ticks) {
    tick_data <- data %>%
      mutate(
        x_start = outer_ring_outer * cos(mid),
        y_start = outer_ring_outer * sin(mid),
        x_end = (outer_ring_outer + tick_length) * cos(mid),
        y_end = (outer_ring_outer + tick_length) * sin(mid)
      )
    
    p <- p +
      geom_segment(data = tick_data,
                   aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                   color = "black", size = 0.5)
  }
  
  annotations_df <- data %>%
    mutate(
      angle_rad = mid,
      angle_deg = (angle_rad * 180 / pi) %% 360,
      x = (outer_ring_outer + tick_length + annotation_distance) * cos(angle_rad),
      y = (outer_ring_outer + tick_length + annotation_distance) * sin(angle_rad),
      angle_text = ifelse(angle_deg > 90 & angle_deg < 270, angle_deg + 180, angle_deg),
      hjust = ifelse(angle_deg > 90 & angle_deg < 270, 1, 0),
      vjust = 0.5
    )
  
  p <- p +
    geom_text(data = annotations_df,
              aes(x = x, y = y, label = !!sym(annotation_col)),
              angle = annotations_df$angle_text,
              hjust = annotations_df$hjust,
              vjust = annotations_df$vjust,
              size = text_size, 
              fontface = "bold", 
              color = "black")
  
  max_radius <- max(outer_ring_outer + tick_length + annotation_distance, 
                    max(abs(annotations_df$x), abs(annotations_df$y))) * 1.1
  
  if (legend_position == "right") {
    x_limits <- c(-max_radius, max_radius * (1 + legend_margin))
    y_limits <- c(-max_radius, max_radius)
  } else if (legend_position == "left") {
    x_limits <- c(-max_radius * (1 + legend_margin), max_radius)
    y_limits <- c(-max_radius, max_radius)
  } else if (legend_position == "top") {
    x_limits <- c(-max_radius, max_radius)
    y_limits <- c(-max_radius, max_radius * (1 + legend_margin))
  } else if (legend_position == "bottom") {
    x_limits <- c(-max_radius, max_radius)
    y_limits <- c(-max_radius * (1 + legend_margin), max_radius)
  } else {
    x_limits <- c(-max_radius, max_radius)
    y_limits <- c(-max_radius, max_radius)
  }
  
  p <- p +
    coord_fixed(xlim = x_limits, ylim = y_limits) +
    theme_void() +
    theme(
      legend.position = legend_position,
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
  
  if (!is.null(title)) {
    p <- p + labs(title = title)
  }
  
  return(p)
}
