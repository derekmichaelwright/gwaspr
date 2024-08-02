#' gg_GWAS_plotly
#'
#' Creates an interactive plotly plot.
#' @param mp ggplot object.
#' @param filename Filename for plot.
#' @return A plotly summary GWAS plot.
#' @export

gg_GWAS_plotly <- function(mp, filename = "GWAS_plotly.html") {
  #
  mpp <- plotly::ggplotly(mp) %>% plotly::layout(showlegend = F)
  htmlwidgets::saveWidget(plotly::as_widget(mpp),
                          filename,
                          selfcontained = T)
  mpp
}
