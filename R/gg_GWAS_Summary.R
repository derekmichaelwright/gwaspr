#' gg_GWAS_Summary
#'
#' [Creates a summary GWAS plot of significant associations.](https://derekmichaelwright.github.io/gwaspr/articles/gg_GWAS_Summary.html)
#' Note: this function requires the GWAS results files to be ordered
#' @param folder Folder containing GWAS results.
#' @param traits The traits to read.
#' @param groups Grouping for the traits. Should be equal length to `traits`.
#' @param threshold Significant threshold.
#' @param sug.threshold Suggestive threshold.
#' @param chr Chromosomes to plot.
#' @param pos1 starting position to plot.
#' @param pos2 ending position to plot.
#' @param models Models to read.
#' @param model.colors Colors for each model.
#' @param shapes The shape values to use for the different models. e.g., 21:25.
#' @param vlines Markers to be labelled with a vertical red line.
#' @param vline.colors colors for each vertical line.
#' @param vline.types lty for each vertical line.
#' @param vline.legend Logical, display of vline color legend.
#' @param genes Table of genes with "Chr", Start" and "End" columns.
#' @param threshold.legend Logical, display of vline color legend.
#' @param hlines Locations for horizontal lines. e.g., hlines = c(1.5,2.5).
#' @param title A title for the plot.
#' @param caption A caption for the plot.
#' @param rowread Number of rows to read for each GWAS results file.
#' @param legend.rows Number of rows for the legend.
#' @param chr.unit Unit for the x-axis. Can be one of c("kbp","100 kbp","Mbp","100 Mbp","Gbp").
#' @param plotHBPvalues Logical, should H.B.P.Values be uses.
#' @param skyline Which skyline type to use. Can be "NYC" or "Kansas". If NULL, it will use the highest P.value.
#' @return A GWAS summary plot.
#' @export

gg_GWAS_Summary <- function(
    folder = "GWAS_Results/",
    traits = list_Traits(folder),
    groups = NULL,
    threshold = round(-log10(0.00000005),1),
    sug.threshold = round(-log10(0.000005),1),
    chr = NULL, pos1 = NULL, pos2 = NULL,
    models =  c("MLM", "MLMM", "FarmCPU", "BLINK", "CMLM", "SUPER"),
    model.colors = gwaspr_ModelColors,
    shapes = 21:25,
    hlines = NULL,
    vlines = NULL,
    vline.colors = rep("red",length(vlines)),
    vline.types = rep(1, length(vlines)),
    vline.legend = T,
    genes = NULL,
    threshold.legend = T,
    title = "Summary of Significant GWAS Results",
    caption = NULL,
    rowread = 2000,
    legend.position = "bottom",
    legend.rows = 1,
    chr.unit = "100 Mbp",
    plotHBPvalues = F,
    skyline = "NYC"
    ) {
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
  #if(is.null(sug.threshold) &
  #   caption == paste0("Significant Threshold = ", threshold, " = Large\nSuggestive = ", sug.threshold, " = Small") ) {
  #   caption <- paste0("Significant Threshold = ", threshold)
  #}
  #
  fnames <- list_Result_Files(folder)
  fnames <- fnames[grepl(paste(models,collapse="|"),fnames)]
  fnames <- fnames[grepl(paste(c(paste0(traits, ".csv"), paste0(traits, "\\(")), collapse="|"), fnames)]
  #
  if(!is.null(skyline)) {
    if(skyline == "NYC") { fnames <- fnames[!grepl("\\(Kansas\\)", fnames)] }
    if(skyline == "Kansas") {
      fnames <- fnames[!(grepl("\\(NYC\\)", fnames) & grepl("FarmCPU", fnames))]
      fnames <- fnames[!(grepl("\\(NYC\\)", fnames) & grepl("BLINK",   fnames))]
    }
  }
  #
  xP <- NULL
  #
  for(i in fnames) {
    xPi <- table_GWAS_Results(folder = folder,
                               fnames = i,
                               threshold = threshold,
                               sug.threshold = sug.threshold,
                               skyline = skyline)
    #
    if(nrow(xPi)>0) { xP <- bind_rows(xP, xPi) }
  }
  #
  if(plotHBPvalues == T) {
    xP <- xP %>% mutate(Pvalue = negLog10_HBP)
  } else {
    xP <- xP %>% mutate(Pvalue = negLog10_P)
  }
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
  xP <- xP %>%
    filter(!is.na(SNP)) %>%
    arrange(Chr, Pos, P.value, Trait) %>%
    mutate(Model = factor(Model, levels = models),
           Trait = factor(Trait, levels = traits),
           Group = NULL) %>%
    filter(!is.na(Trait)) %>%
    arrange(desc(Model))
  #
  xG <- read.csv(paste0(folder, fnames[grepl(traits[1],fnames)][1])) %>%
    mutate(Trait = xP$Trait[1],
           Trait = factor(Trait, levels = traits),
           Group = NULL, Model = NULL)
  #
  if(!is.null(chr)) {
    xG <- xG %>% filter(Chr %in% chr)
    xP <- xP %>% filter(Chr %in% chr)
    if(!is.null(genes)) { genes <- genes %>% filter(Chr %in% chr) }

    if(!is.null(pos1)) {
      xG <- xG %>% filter(Pos > pos1, Pos < pos2)
      xP <- xP %>% filter(Pos > pos1, Pos < pos2)
      if(!is.null(genes)) { genes <- genes %>% filter(Start > pos1, End < pos2) }
    }
  }
  #
  if(!is.null(groups)) {
    xP <- xP %>%
      mutate(Group = plyr::mapvalues(Trait, traits, groups),
             Group = factor(Group, levels = unique(groups)))
    xG <- xG %>%
      mutate(Group = plyr::mapvalues(Trait, traits, groups),
             Group = factor(Group, levels = unique(groups)))
  }
  xP <- xP %>%
    mutate(Threshold = ifelse(Threshold == "Significant",
                              paste("Significant =",threshold), Threshold),
           Threshold = ifelse(Threshold == "Suggestive",
                              paste("Suggestive =",sug.threshold), Threshold),
           Threshold = factor(Threshold))
  #
  if(is.null(pos2)) { myBreaks <- 0:(round(max(xP$Pos)/x.unit))
  } else { myBreaks <- (round(pos1/x.unit)):(round(pos2/x.unit)) }
  #
  mp <- ggplot(xP) + geom_blank(data = xG)
  #
  if(!is.null(vlines)) {
    xGM <- xG %>% filter(SNP %in% vlines) %>%
      mutate(SNP = factor(SNP, levels = vlines)) %>%
      select(SNP, Chr, Pos)
    #
    mp <- mp +
      geom_vline(data = xGM, alpha = 0.5,
                 aes(xintercept = Pos / x.unit, color = SNP)) +
      scale_color_manual(name = "Marker", values = vline.colors) +
      scale_linetype_manual(name = "Marker", values = vline.types)
  }
  #
  if(!is.null(hlines)) {
    mp <- mp + geom_hline(yintercept = hlines, alpha = 0.7)
  }
  #
  mp <- mp +
    geom_point(color = "black", alpha = 0.5,
               aes(x = Pos / x.unit, y = Trait, shape = Model,
                   fill = Model, size = Threshold,
                   key1 = SNP, key2 = negLog10_P, key3 = negLog10_HBP))
  #
  if(!is.null(groups)) {
    mp <- mp + facet_grid(Group ~ Chr, scales = "free", space = "free")
  } else { mp <- mp + facet_grid(. ~ Chr, scales = "free", space = "free") }
  #
  mp <- mp +
    scale_fill_manual(values = model.colors, breaks = models) +
    scale_shape_manual(values = shapes, breaks = models) +
    scale_size_manual(name = NULL, values = c(2.25,0.75)) +
    scale_y_discrete(limits = rev) +
    scale_x_continuous(breaks = myBreaks, minor_breaks = myBreaks) +
    theme_gwaspr(legend.position = legend.position) +
    guides(shape = guide_legend(nrow = legend.rows, override.aes = list(size = 4)),
           color = guide_legend(nrow = legend.rows),
           fill = guide_legend(nrow = legend.rows),
           size = guide_legend(nrow = legend.rows) ) +
    labs(title = title, y = NULL, x = chr.unit, caption = caption)
  #
  if(vline.legend == F) {
    mp <- mp + guides(color = vline.legend)
  }
  #
  if(threshold.legend == F) {
    mp <- mp + guides(color = threshold.legend)
  }
  #
  if(!is.null(genes)) {
    mp <- mp +
      geom_rect(data = genes, ymin = -0.5, ymax = 0.5, color = "black",
                aes(xmin = Start/x.unit, xmax = End/x.unit, key = Name) )
  }
  mp
}

#folder = "GWAS_Results/"; traits = list_Traits(folder); groups = NULL
#threshold = round(-log10(0.00000005),1); sug.threshold = round(-log10(0.000005),1)
#chr = NULL; pos1 = NULL; pos2 = NULL
#models =  c("MLM", "FarmCPU", "BLINK", "MLMM", "GLM", "CMLM", "SUPER")
#model.colors = c("darkgreen", "darkorange3", "steelblue", "darkred", "darkorchid4", "burlywood4", "darkseagreen4")
#shapes = 21:25; hlines = NULL; vlines = NULL; vline.colors = rep("red",length(vlines));
#vline.types = rep(1, length(vlines)); vline.legend = T; title = "Summary of Significant GWAS Results";
#caption = paste0("Significant Threshold = ", threshold, " = Large\nSuggestive Threshold = ", sug.threshold," = Small")
#rowread = 2000; legend.position = "bottom"; legend.rows = 1; plotHBPvalues = F; skyline = F
#threshold.legend = T
