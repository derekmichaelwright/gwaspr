#' theme_gwaspr_pie.
#'
#' ggplot theme.
#' @param x ggplot
#' @param bgFill background fill color
#' @param lineSize border line size
#' @param stripFill strip background color
#' @param caption.adj hjust for the caption
#' @param caption.size size of the caption
#' @return ggplot with agData theme for pie graphs
#' @export

theme_gwaspr_pie <- function(x, bgFill = "white", lineSize = 0.75, stripFill = "white",
                             caption.adj = 1, caption.size = 7, ...) {
  theme(panel.background = element_rect(color = "black", fill = bgFill, linewidth = lineSize),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = lineSize),
        strip.background = element_rect(color = "black", fill = stripFill, linewidth = lineSize),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.key = element_rect(color = NA),
        plot.caption = element_text(hjust = caption.adj, size = caption.size),
        ...)
}
