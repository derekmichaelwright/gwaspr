#' gg_GWAS_Summary_plotly
#'
#' Creates a summary GWAS plot of significant associations.
#' @param mp Plot from gg_GWAS_Summary.
#' @param width Plot width.
#' @param height Plot height.
#' @param filename Filename for plot.
#' @return A plotly summary GWAS plot.
#' @export

gg_GWAS_Summary_plotly <- function(mp, width = 10, height = 8,
                                   filename = "GWAS_Summary.html") {
  #
  mpp <- plotly::ggplotly(mp)
  htmlwidgets::saveWidget(plotly::as_widget(mp),
                          filename,
                          knitrOptions = list(fig.width = width, fig.height = height),
                          selfcontained = T)
  mpp
}
