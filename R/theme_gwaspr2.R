#' theme_gwaspr2.
#'
#' ggplot theme.
#' @param x ggplot
#' @param linesize border line size
#' @param bgFill background fill color
#' @param stripFill strip background color
#' @param lineColor color of axis lines
#' @return ggplot with gwaspr theme
#' @export

theme_gwaspr2 <- function(x, bgFill = "grey95", lineColor = "white", linesize = 0.75, stripFill = "white", ...) {
  theme(panel.background = element_rect(color = "black", fill = bgFill, size = linesize),
        panel.grid = element_line(color = lineColor),
        panel.border = element_rect(color = "black", fill = NA, size = linesize),
        strip.background = element_rect(color = "black", fill = stripFill, size = linesize),
        legend.key = element_rect(color = NA),
        ...)
}
