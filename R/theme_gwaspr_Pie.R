#' theme_gwaspr_pie.
#'
#' ggplot theme.
#' @param x ggplot
#' @return ggplot with agData theme for pie graphs
#' @export

theme_agData_pie <- function(x, bgFill = "white", lineColor = "grey95", linesize = 0.75, stripFill = "white", ...) {
  theme(panel.background = element_rect(color = "black", fill = bgFill, size = linesize),
        panel.grid = element_blank(), # panel.grid = element_line(color = lineColor),
        panel.border = element_rect(color = "black", fill = NA, size = linesize),
        strip.background = element_rect(color = "black", fill = stripFill, size = linesize),
        #panel.background = element_blank(),
        #panel.border = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.key = element_rect(color = NA),
        ...)
}
