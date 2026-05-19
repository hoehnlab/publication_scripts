library(ggplot2)
library(patchwork)

colors = RColorBrewer::brewer.pal(3,"Set1")
names(colors) = c("GC", "Other", "Ambig.")
eo_line_color_scale <- scale_color_gradient(
  low = colors["Other"],
  high = colors["GC"],
  breaks = c(0, 0.25, 0.5, 0.75, 1),
  limits = c(0, 1)
)
eo_node_fill_scale <- scale_fill_gradient(
  low = colors["Other"],
  high = colors["GC"],
  breaks = c(0, 0.25, 0.5, 0.75, 1),
  limits = c(0, 1)
)

.get_xvals <- function(p) {
  xmax <- max(p$data$x)
  xmin <- min(p$data$x)
  return(c(xmin, xmax))
}

# Extract x-range from a ggplot object
.get_xrange <- function(p) {
  b <- ggplot_build(p)
  xr <- b$layout$panel_params[[1]]$x.range
  diff(xr)
}

font_size <- theme(axis.text = element_text(size=7),
                   axis.title = element_text(size=9),
                   plot.title = element_text(size=9),
                   strip.text = element_text(size=8),
                   legend.text = element_text(size=7),
                   legend.title = element_text(size=8))

scale_x_timeline <- function(
    xlim,
    extra=50,
    y0 = -1,
    minor_length = 0.6,
    major_length = 1.0,
    cross = TRUE,
    expand = c(0, 0)
) {
  
  stopifnot(length(xlim) == 2)
  
  xmin <- xlim[1]
  xmax <- xlim[2]
  
  major_x <- c(
    min(xmin, 0),
    max(xmin, 0),
    xmax
  )

  
  make_df <- function(x, length) {
    if (cross) {
      data.frame(
        x = x,
        xend = x,
        y = y0 - length,
        yend = y0 + length
      )
    } else {
      data.frame(
        x = x,
        xend = x,
        y = y0,
        yend = y0 - length
      )
    }
  }
  
  list(
    geom_segment(
      data = make_df(major_x, major_length),
      aes(x = x, xend = x, y = y, yend = yend),
      lineend = "square",
      linewidth = 0.25,
      inherit.aes = FALSE
    ),
    
    # scale
    scale_x_continuous(
      limits = c(min(xmin, 0)-extra,xmax+extra),
      breaks = major_x,
      labels = function(x) ifelse(x %in% c(0, xmax), x, ""),
      expand = c(0,0)
    ),
    
    geom_segment(data=data.frame(test=1), x = xmin, xend = xmax, y = y0, yend = y0, lineend="square", linewidth=0.25, inherit.aes = FALSE),
    
    # hide default ticks so only yours show
    theme(axis.ticks = element_blank(), axis.title.x = element_blank(), axis.line.x = element_blank(), axis.text.x = element_text(margin = margin(t = 0.8 * 11/4), vjust = 1)),
    
    font_size,
    theme(axis.text.x = element_text(color="black"))
  )
}

# Main helper: align plots so 1 data unit = same physical length
align_plots_by_xscale <- function(..., extra=50, normalize = TRUE) {
  plots <- list(...)
  
  plots <- lapply(plots, function(p) {
    xlim <- .get_xvals(p)
    p + scale_x_timeline(xlim = xlim, extra = extra)
  })
  
  widths <- vapply(plots, .get_xrange, numeric(1))
  
  if (normalize) {
    widths <- widths / min(widths)
  }
  
  wrap_plots(plots, nrow = 1, widths = widths)
}