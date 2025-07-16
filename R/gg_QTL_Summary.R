#' gg_QTL_Summary
#'
#' Creates a summary QTL plot of significant associations.
#' Note: this function requires the GWAS results files to be ordered
#' @param xx Table of QTL results.
#' @param myG genotype file
#' @param groups Grouping for the traits. Should be equal length to `traits`.
#' @return A QTL summary plot.
#' @export

gg_QTL_Summary <- function(
    xx = myQTLs,
    myG = myG,
    title = "Summary of QTL Results",
    caption = "derek was here"
    ) {
  #
  ggplot(xx, aes(x = Pos / 100000000)) + 
    geom_blank(data = myG) +
    geom_segment(aes(y = Pop, yend = Pop, xend = Pos_2 / 100000000, color = Pop), size = 2 ) +
    facet_grid(trait ~ Chr, space = "free_x", scales = "free_x") +
    scale_color_manual(values = c("steelblue","darkred","purple4")) + 
    theme_gwaspr(legend.position = "none") +
    labs(title = title, caption = caption,
         x = "100 MBp", y = "Population")
}
