#' theme_gwaspr_col.
#'
#' ggplot theme.
#' @param x ggplot
#' @param bgFill background fill color
#' @param lineColor color of axis lines
#' @param lineSize border line size
#' @param stripFill strip background color
#' @param caption.adj hjust for the caption
#' @param caption.size size of the caption
#' @return ggplot with gwaspr theme
#' @export

theme_gwaspr_col <- function(x, bgFill = "white", lineColor = "grey90", lineSize = 0.75, stripFill = "white",
                             caption.adj = 1, caption.size = 7, ...) {
    theme(panel.background = element_rect(color = "black", fill = bgFill, linewidth = lineSize),
          panel.grid = element_line(color = lineColor),
          panel.border = element_rect(color = "black", fill = NA, linewidth = lineSize),
          strip.background = element_rect(color = "black", fill = stripFill, linewidth = lineSize),
          legend.key = element_rect(color = NA),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          plot.caption = element_markdown(hjust = caption.adj, size = caption.size),
          ...)
}
