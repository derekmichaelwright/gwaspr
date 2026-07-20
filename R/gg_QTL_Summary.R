#' gg_QTL_Summary
#'
#' Creates a summary QTL plot of significant associations.
#' @param myG Genetic map. Must contain the following 3 columns c("Marker","Chr","Pos").
#' @param myQ Table of QTL results. Must contain the following 5 columns c("Trait","Chr","Pos","Pos_lo","Pos_hi","lod").
#' @param title Custom title for the plot.
#' @param lodFill Logical. If true, fill color will be a gradient based on lod values.
#' @param fillColor Color for filling points.
#' @param fillColor_low If `lodfill = T`, This will be the low end of the gradeint color scale for filling points.
#' @param xLab Custom x lab.
#' @param facetLab Custom label for facets.
#' @return A QTL summary plot.
#' @export

gg_QTL_Summary <- function(
    myG, myQ,
    title = "Summary of QTL Results",
    lodFill = F,
    fillColor = "darkgreen",
    fillColor_low = "grey50",
    xLab = "Pos",
    facetLab = "lg"
    ) {
  #
  mp <- ggplot(myG, aes(x = Pos)) +
    geom_blank() +
    facet_grid(. ~ paste(facetLab, Chr), scales = "free_x", space = "free_x") +
    theme_gwaspr() +
    labs(title = title, x = xLab, y = NULL) +
    #
    geom_segment(data = myQ, linewidth = 1, alpha = 0.7,
                 aes(x = Pos_lo, xend = Pos_hi, y = Trait, yend = Trait))
  if(lodFill == T) {
    mp <- mp +
      geom_point(data = myQ, aes(x = Pos, y = Trait, fill = lod), pch = 23, size = 2.5) +
      scale_fill_gradient(low = fillColor_low, high = fillColor)
    } else {
      mp <- mp +
        geom_point(data = myQ, aes(x = Pos, y = Trait), fill = fillColor, pch = 23, size = 2.5)
      }
  mp
}
