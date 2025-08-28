#' gg_GWAS_Hits
#'
#' Creates a summary GWAS plot of significant associations.
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
#' @return A GWAS Hits plot.
#' @export

gg_GWAS_Hits <- function(
    xx, xG, xCV = NULL,
    traits,
    range = 2000000,
    title = "",
    sigMin = 0,
    models =  c("BLINK", "FarmCPU", "MLMM", "MLM", "GLM", "CMLM", "SUPER"),
    model.colors = c("steelblue", "darkorange3", "darkred", "darkgreen", "darkorchid4", "burlywood4", "darkseagreen4"),
    model.shapes = c(21:25,21:22),
    vlines = NULL,
    vline.colors = rep("red", length(vlines)),
    vline.types = rep(1, length(vlines)),
    legend.rows = 1
    ) {
  #
  xx <- xx %>% filter(Trait %in% traits, Model %in% models) %>% mutate(Hits = NA)
  xG <- xG %>% dplyr::select(SNP=1, Chr=3, Pos=4)
  #
  for(i in 1:nrow(xx)) {
    myChr <- xx$Chr[i]
    myPos <- xx$Pos[i]
    myMod <- xx$Model[i]
    xi <- xx %>%
      filter(Chr == myChr, Model == myMod,
             Pos > myPos-range, Pos < myPos+range)
    xx$Hits[i] <- nrow(xi)
    xi <- xi %>% filter(!(SNP == xx$SNP[i] & Trait == xx$Trait[i]))
    xx$Pos[xx$SNP %in% xi$SNP & xx$Trait %in% xi$Trait] <- NA
  }
  xx <- xx %>% filter(!is.na(Pos), Hits > sigMin) %>%
    mutate(Model = factor(Model, levels = models))
  #
  if(!is.null(xCV)) { xx <- xx %>% filter(CV %in% xCV) }
  ymax <- max(xx$Hits)
  #
  mp <- ggplot(xx, aes(x = Pos / 100000000) ) +
    geom_blank(data = xG)
  #
  if(!is.null(vlines)) {
    xGM <- xG %>%
      filter(SNP %in% vlines) %>%
      mutate(SNP = factor(SNP, levels = vlines)) %>%
      arrange(SNP)
    mp <- mp +
      geom_vline(data = xGM, alpha = 0.7,
                 aes(xintercept = Pos/1e+08, color = SNP, lty = SNP)) +
      scale_color_manual(name = NULL, values = vline.colors) +
      scale_linetype_manual(name = NULL, values = vline.types)
  }
  mp <- mp +
    geom_point(aes(y = Hits, size = Hits, key1 = SNP,
                   fill = Model, shape = Model), alpha = 0.7) +
    facet_grid(. ~ paste("Chr", Chr), space = "free_x", scales = "free_x") +
    scale_shape_manual(values = model.shapes)+
    scale_fill_manual(values = model.colors)+
    scale_size_continuous(range = c(0.5,3)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 1),
                       limits = c(0,ymax)) +
    theme_gwaspr() +
    theme(legend.position = "bottom") +
    #guides(color = F, size = F) +
    guides(fill = guide_legend(nrow = legend.rows, override.aes = list(size = 1.5)),
           color = guide_legend(nrow = legend.rows, byrow = T),
           size = "none",
           shape = guide_legend(nrow = legend.rows, byrow = T) ) +
    labs(title = title, x = "100 Mbp", y = "Significant Associations")
  mp
}

