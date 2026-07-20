#' gg_QTL_Summary_Groups
#'
#' Creates a summary QTL plot of significant associations.
#' @param myG Genetic map. Must contain the following 3 columns c("Marker","Chr","Pos").
#' @param myQ Table of QTL results. Must contain the following 5 columns c("Trait","Chr","Pos","Pos_lo","Pos_hi","lod"). then it must also include grouping columns set by the user for `yGroup`,`facetGroup` & `colorGroup`.
#' @param title Custom title for the plot.
#' @param yGroup Name of column for y axis grouping.
#' @param yLab Custom y-axis label.
#' @param facetGroup Name of column for facet grouping.
#' @param colorGroup Name of column for color grouping.
#' @param colorName Title for the color legend.
#' @param fillColors Color for filling points.
#' @param xLab Custom x-axis label.
#' @param facetLab Custom label for facets.
#' @return A QTL summary plot.
#' @export

gg_QTL_Summary_Groups <- function(
    myQ, myG,
    title = "Summary of QTL Results",
    yGroup,
    yLab = NULL,
    facetGroup,
    colorGroup,
    colorName = NULL,
    fillColors = gwaspr_Colors,
    xLab = "Pos",

    facetLab = "lg"
    ) {
  #
  mp <- ggplot(myG, aes(x = Pos)) +
    geom_blank() +
    facet_grid(get(facetGroup) ~ paste(facetLab, Chr), scales = "free_x", space = "free_x") +
    theme_gwaspr(legend = "bottom") +
    labs(title = title, x = xLab, y = yLab) +
    #
    geom_segment(data = myQ, linewidth = 1, alpha = 0.7,
                 aes(x = Pos_lo, xend = Pos_hi, y = get(yGroup), yend = get(yGroup))) +
    geom_point(data = myQ, aes(x = Pos, y = get(yGroup), fill = get(colorGroup)), pch = 23, size = 2.5) +
    scale_fill_manual(name = colorName, values = fillColors)
  mp
}
