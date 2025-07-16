#' gg_plotly
#'
#' Creates an interactive plotly plot. Note: This function will require `plotly` and `htmlwidgets` packages.
#' @param mp ggplot object.
#' @param filename Filename for plot.
#' @param showlegend Logical, should the legend be displayed.
#' @return A plotly summary GWAS plot.
#' @export

gg_plotly <- function(mp, filename = "my_gg_plotly.html", showlegend = T) {
  #
  mpp <- plotly::ggplotly(mp) %>%
    plotly::layout(showlegend = showlegend)
  #
  htmlwidgets::saveWidget(plotly::as_widget(mpp),
                          filename,
                          selfcontained = T)
  mpp
}
