#' gg_Manhattan_Zoom
#'
#' Creates a manhattan plot zoomed in to a particular region.
#' @param folder Folder containing GWAS results.
#' @param trait The trait to read.
#' @param chrom Chromosome to plot.
#' @param pos1 Start position on chromosome.
#' @param pos2 End position on chromosome.
#' @param title A title for the plot.
#' @param threshold Significant Threshold.
#' @param sug.threshold Suggested threshold.
#' @param markers Markers to be labelled.
#' @param labels Labels to be used for markers.
#' @param vlines Markers which will be used as a location for a vertical lines.
#' @param vline.colors colors for each vertical line.
#' @param sig.col Color for significant assoctiations.
#' @param models Models to read.
#' @param facet Logical, whether or not to produce a facetted or multi-model plot. Default is `facet = F`.
#' @param highlight.sig Logical, whether or not to highlight significant associations with a black circle. Used if `facet = F`.
#' @param legend.rows Number of rows for the legend.
#' @param plotHBPvalues Logical, should H.B.P.Values be uses.
#' @return A manhattan plot.
#' @export

gg_Manhattan_Zoom <- function(
    folder = "GWAS_Results/",
    trait = list_Traits(folder)[1],
    chrom, pos1, pos2,
    title = trait,
    threshold = NULL,
    sug.threshold = NULL,
    markers = NULL,
    labels = markers,
    vlines = markers,
    vline.colors = "red",
    sig.col = "red",
    models = c("MLM", "FarmCPU", "BLINK", "MLMM", "GLM", "CMLM", "SUPER"),
    model.colors = c("darkgreen", "darkorange3", "steelblue", "darkred", "darkorchid4", "burlywood4", "darkseagreen4"),
    facet = F,
    highlight.sig = F,
    highlight.marker.color = "red",
    legend.rows = 1,
    plotHBPvalues = F
    ) {
  #
  # Read in files
  #
  fnames <- list_Result_Files(folder)
  fnames <- fnames[grepl(paste(paste0(trait,".csv"),collapse="|"),fnames)]
  fnames <- fnames[grepl(paste(models,collapse="|"),fnames)]
  #
  xx <- NULL
  #
  for(i in fnames) {
    mod <- substr(i, gregexpr("GWAS_Results", i)[[1]][1]+13, gregexpr(".csv", i)[[1]][1]-1)
    mod <- substr(mod, 1, gregexpr("\\.", mod)[[1]][1]-1)
    xi <- read.csv(paste0(folder, i))
    if(sum(colnames(xi)=="nobs")>0) { xi <- select(xi, -nobs) }
    xi <- xi %>% mutate(Model = mod, negLog10_P = -log10(P.value))
    xx <- bind_rows(xx, xi)
  }
  #
  # Prep data
  #
  xx <- xx %>% mutate(Model = factor(Model, levels = models)) %>%
    filter(Chr == chrom, Pos > pos1, Pos < pos2, !is.na(Model))
  #
  if (is.null(threshold)) {
    threshold <- -log10(0.05/nrow(xi))
  }
  #
  xx <- xx %>%
    mutate(Sig.level = ifelse(negLog10_P >= threshold, "Sig","Not Sig"))
  if(!is.null(sug.threshold)) {
    xx <- xx %>%
      mutate(Sig.level = ifelse(negLog10_P < threshold & negLog10_P >= sug.threshold, "Sug", Sig.level))
  }
  #
  if(plotHBPvalues == T) {
    xx <- xx %>% mutate(pvals = -log10(H.B.P.Value))
  } else {
    xx <- xx %>% mutate(pvals = negLog10_P)
  }
  #
  x2 <- xx %>% filter(negLog10_P >= threshold)
  x3 <- xx %>% filter(SNP %in% markers)
  #
  # Start Plots
  #
  mp <- ggplot(xx, aes(x = Pos / 1000000, y = pvals))
  #
  # Add vlines
  #
  if(!is.null(vlines)) {
    mp <- mp +
      geom_vline(data = xx %>% filter(SNP %in% markers),
                 aes(xintercept = Pos / 1000000, color = SNP),
                 alpha = 0.5) +
      scale_color_manual(name = NULL, values = vline.colors)
  }
  #
  # Add threshold lines
  #
  mp <- mp +
    geom_hline(yintercept = threshold, color = "red", alpha = 0.8, linewidth = 0.5) +
    geom_hline(yintercept = sug.threshold, color = "blue", alpha = 0.8, linewidth = 0.5)
  #
  # Plot the Rest
  #
  mp <- mp +
    geom_hline(yintercept = threshold, color = "red", alpha = 0.6) +
    geom_point(aes(fill = Model, size = Sig.level), pch = 21, color = alpha("white",0)) +
    scale_size_manual(name = NULL, values = c(0.5,1.25,0.75), guide = "none") +
    guides(color = guide_legend(nrow = legend.rows, byrow = T, override.aes = list(alpha = 1))) +
    theme_gwaspr(legend.position = "bottom",
                 axis.title.y = element_markdown()) +
    labs(title = trait, y = "-log<sub>10</sub>(*p*)", x = "Mbp")
  if(highlight.sig == T) {
    mp <- mp + geom_point(data = x2, fill = alpha("white", 0), pch = 21, size = 1.25)
  }
  #
  # Add Marker labels
  #
  if(!is.null(markers)) {
    xx <- xx %>% mutate(Label = ifelse(SNP %in% markers, SNP, NA),
                        Label = plyr::mapvalues(Label, markers, labels))
    mp <- mp + geom_text_repel(data = xx %>% filter(SNP %in% markers),
                               aes(label = Label), size = 2)
  }
  #
  mp <- mp +
    geom_point(data = x3, aes(fill = Model, size = Sig.level), color = highlight.marker.color, pch = 21)
  #
  # Facet plot
  #
  if(facet == T) {
    mp <- mp + facet_grid(Model ~ Chr, scales = "free") +
      scale_fill_manual(name = NULL, values = alpha(model.colors,0.8), guide = "none")
  } else {
    mp <- mp +
      scale_fill_manual(name = NULL, values = alpha(model.colors,0.8)) +
      guides(fill = guide_legend(nrow = legend.rows, override.aes = list(size = 2)))
  }
  mp
}
