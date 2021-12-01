#' gg_GWAS_Summary_plotly
#'
#' Creates a summary GWAS plot of significant associations.
#' @param mp Plot from gg_GWAS_Summary.
#' @param width Plot width.
#' @param height Plot height.
#' @param filename Filename for plot.
#' @return A plotly summary GWAS plot.
#' @export

gg_GWAS_Summary_plotly <- function(mp, filename = "GWAS_Summary.html", width = 10, height = 8 ) {
  #
  mpp <- ggplotly(mp)
  saveWidget(as_widget(mpp),
             filename,
             knitrOptions = list(fig.width = width, fig.height = height),
             selfcontained = T)
  mpp
}
