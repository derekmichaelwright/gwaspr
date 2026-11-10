#' gg_GWAS_Hits
#'
#' [Creates a summary GWAS plot of significant associations.](https://derekmichaelwright.github.io/gwaspr/articles/gg_GWAS_Hits.html)
#' @param xx Table of significant GWAS results. See ?table_GWAS_Results().
#' @param xG Genotype data in hapmap format.
#' @param xCV Optional, for filtering if you have a "CV" column.
#' @param traits List of traits to use.
#' @param range Range for binning GWAS hits.
#' @param title Title for horizontal facet.
#' @param sigMin Minimum number of hits to plot.
#' @param models Models to read.
#' @param model.colors Colors for each model.
#' @param vlines Markers to be labelled with a vertical red line.
#' @param vline.colors colors for each vertical line.
#' @param vline.types lty for each vertical line.
#' @param legend.rows Number of rows for the legend.
#' @param facet Logical. facet by GWAS model.
#' @param caption Logical. should a caption be added.
#' @return A GWAS Hits plot.
#' @export

gg_GWAS_Hits <- function(
    folder = "GWAS_Results/",
    traits = list_Traits(folder),
    models =  c("BLINK", "FarmCPU", "MLMM", "MLM", "GLM", "CMLM", "SUPER"),
    model.colors = gwaspr_ModelColors,
    model.shapes = c(21:25,21:22),
    groups = NULL,
    chr = NULL, pos1 = NULL, pos2 = NULL,
    range = 2000000,
    title = "",
    sigMin = 0,
    vlines = NULL,
    vline.colors = rep("red", length(vlines)),
    vline.types = rep(1, length(vlines)),
    vline.legend = T,
    legend.rows = 1,
    chr.unit = "100 Mbp",
    facet = F,
    caption = T,
    skyline = "NYC"
    ) {
  #
  #
  check1 <- run_Summary(folder = folder) %>% dropNAcol()
  check2 <- is_Ordered(folder = folder) %>%
    select(colnames(check1))
  if(sum(is.na(check2)) > 0) {
    warning("Some of your GWAS results files might not be ordered by pvalue")
    warning("use order_GWAS_Results() before making these summary plots")
    warning("use is_ordered() to check if GWAS results are properly ordered for use with this function")
  }
  #
  #myCaption <- paste("Bin size =", range, "\n", "Traits =",
  #                   paste(unique(xx$Trait), collapse = ", "))
  #
  if(chr.unit == "1 kbp")   { x.unit = 1000 }
  if(chr.unit == "10 kbp")  { x.unit = 10000 }
  if(chr.unit == "50 kbp")  { x.unit = 50000 }
  if(chr.unit == "100 kbp") { x.unit = 100000 }
  if(chr.unit == "1 Mbp")   { x.unit = 1000000 }
  if(chr.unit == "10 Mbp")  { x.unit = 10000000 }
  if(chr.unit == "50 Mbp")  { x.unit = 50000000 }
  if(chr.unit == "100 Mbp") { x.unit = 100000000 }
  if(chr.unit == "1 Gbp")   { x.unit = 1000000000 }
  if(!chr.unit %in% c("1 kbp", "10 kbp", "50 kbp", "100 kbp",
                      "1 Mbp", "10 Mbp", "50 Mbp", "100 Mbp", "1 Gbp")) { print("error in chr.unit") }
  #
  #
  xR <- table_GWAS_Results(
    folder = folder,
    traits = traits,
    threshold = 6.7 )
  #
  xRS <- NULL
  #i<-"MLM"
  for(i in unique(models)) {
    xRi <- xR %>% filter(Model == i)
    if(nrow(xRi)>0) {
      xRSi <- table_GWAS_Results_Summary(
        xx = xRi,
        binMarkers = F,
        binSize = 100000 )
      xRS <- bind_rows(xRS, xRSi)
    }
  }
  #
  #mp <- ggplot()# + geom_blank(data = xG)
  #
  mp <- ggplot() +
    geom_point(data = xRS, aes(x = Pos / x.unit, y = Hits)) +
    facet_grid(Models ~ Chr, scales = "free", space = "free_x") +
    scale_shape_manual(values = model.shapes) +
    #scale_fill_manual(values = model.colors) +
    theme(legend.position = "bottom") +
    theme_gwaspr() +
    guides(fill = guide_legend(nrow = legend.rows, override.aes = list(size = 1.5)),
           color = guide_legend(nrow = legend.rows, byrow = T),
           size = "none",
           shape = guide_legend(nrow = legend.rows, byrow = T) ) +
    labs(title = title, x = chr.unit, y = "Significant Associations")
  #
  if(!is.null(vlines)) {
    xGM <- xG %>%
      filter(SNP %in% vlines) %>%
      mutate(SNP = factor(SNP, levels = vlines)) %>%
      arrange(SNP)
    mp <- mp +
      geom_vline(data = xGM, alpha = 0.7,
                 aes(xintercept = Pos / x.unit, color = SNP, lty = SNP)) +
      scale_color_manual(name = NULL, values = vline.colors) +
      scale_linetype_manual(name = NULL, values = vline.types)
  }
  #
  mp
}

