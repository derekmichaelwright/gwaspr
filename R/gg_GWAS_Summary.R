#' gg_GWAS_Summary
#'
#' Creates a summary GWAS plot of significant associations.
#' Note: this function requires the GWAS results files to be ordered
#' @param folder Folder containing GWAS results.
#' @param traits The traits to read.
#' @param groups Grouping for the traits. Should be equal length to `traits`.
#' @param threshold Significant threshold.
#' @param sug.threshold Suggestive threshold.
#' @param chroms Chromosomes to plot.
#' @param pos1 starting position to plot.
#' @param pos2 ending position to plot.
#' @param models Models to read.
#' @param model.colors Colors for each model.
#' @param shapes The shape values to use for the different models. e.g., 21:25.
#' @param vlines Markers to be labelled with a vertical red line.
#' @param vline.colors colors for each vertical line.
#' @param vline.types lty for each vertical line.
#' @param vline.legend Logical, display of vline color legend.
#' @param hlines Locations for horizontal lines. e.g., hlines = c(1.5,2.5).
#' @param title A title for the plot.
#' @param caption A caption for the plot.
#' @param rowread Number of rows to read for each GWAS results file.
#' @param legend.rows Number of rows for the legend.
#' @param plotHBPvalues Logical, should H.B.P.Values be uses.
#' @param skyline Which skyline type to use. Can be "NYC" or "Kansas". If left NULL, it will use the highest P.value.
#' @return A GWAS summary plot.
#' @export

gg_GWAS_Summary <- function(
    folder = "GWAS_Results/",
    traits = list_Traits(folder),
    groups = NULL,
    threshold = round(-log10(0.00000005),1),
    sug.threshold = round(-log10(0.000005),1),
    chroms = NULL, pos1 = NULL, pos2 = NULL,
    models =  c("MLM", "MLMM", "FarmCPU", "BLINK",  "GLM", "CMLM", "SUPER"),
    model.colors = gwaspr_Colors,
    shapes = 21:25,
    hlines = NULL,
    vlines = NULL,
    vline.colors = rep("red",length(vlines)),
    vline.types = rep(1, length(vlines)),
    vline.legend = T,
    title = "Summary of Significant GWAS Results",
    caption = paste0("Significant Threshold = ", threshold, " = Large\nSuggestive Threshold = ", sug.threshold," = Small"),
    rowread = 2000,
    legend.position = "bottom",
    legend.rows = 1,
    plotHBPvalues = F,
    skyline = "Kansas"
    ) {
  #
  check1 <- is_ran(folder = folder) %>% dropNAcol()
  check2 <- is_ordered(folder = folder) %>%
    dplyr::select(colnames(check1))
  if(sum(is.na(check2)) > 0) {
    warning("Some of your GWAS results files might not be ordered by pvalue")
    warning("use order_GWAS_Results() before making these summary plots")
    warning("use is_ordered() to check if GWAS results are properly ordered for use with this function")
  }
  #
  if(is.null(sug.threshold) &
     caption == paste0("Significant Threshold = ", threshold, " = Large\nSuggestive = ", sug.threshold, " = Small") ) {
     caption <- paste0("Significant Threshold = ", threshold)
  }
  #
  fnames <- list_Result_Files(folder)
  fnames <- fnames[grepl(paste(models,collapse="|"),fnames)]
  fnames <- fnames[grepl(paste(c(paste0(traits, ".csv"), paste0(traits, "\\(")), collapse="|"), fnames)]
  #
  if(!is.null(skyline)) {
    if(skyline == "NYC") { fnames <- fnames[!grepl("\\(Kansas\\)", fnames)] }
    if(skyline == "Kansas") {
      fnames <- fnames[!grepl("\\(NYC\\)&FarmCPU", fnames)]
      fnames <- fnames[!grepl("\\(NYC\\)&BLINK", fnames)]
    }
  }
  #
  myP <- NULL
  #
  for(i in fnames) {
    myPi <- table_GWAS_Results(folder = folder,
                               fnames = i,
                               threshold = threshold,
                               sug.threshold = sug.threshold,
                               skyline = skyline)
    #
    if(nrow(myPi)>0) { myP <- bind_rows(myP, myPi) }
  }
  #
  if(plotHBPvalues == T) {
    myP <- myP %>% mutate(Pvalue = negLog10_HBP)
  } else {
    myP <- myP %>% mutate(Pvalue = negLog10_P)
  }
  #
  myP <- myP %>%
    filter(!is.na(SNP)) %>%
    arrange(Chr, Pos, P.value, Trait) %>%
    mutate(Model = factor(Model, levels = models),
           Trait = factor(Trait, levels = traits),
           Group = NULL) %>%
    filter(!is.na(Trait)) %>%
    arrange(desc(Model))
  #
  myG <- read.csv(paste0(folder, fnames[1])) %>%
    mutate(Trait = myP$Trait[1],
           Trait = factor(Trait, levels = traits),
           Group = NULL, Model = NULL)
  #
  if(!is.null(chroms)) {
    myG <- myG %>% filter(Chr %in% chroms, Pos > pos1, Pos < pos2)
    x1 <- x1 %>% filter(Chr %in% chroms, Pos > pos1, Pos < pos2)
    x2 <- x2 %>% filter(Chr %in% chroms, Pos > pos1, Pos < pos2)
  }
  #
  if(!is.null(groups)) {
    myP <- myP %>%
      mutate(Group = plyr::mapvalues(Trait, traits, groups),
             Group = factor(Group, levels = rev(unique(myGroups))))
    myG <- myG %>%
      mutate(Group = plyr::mapvalues(Trait, traits, groups),
             Group = factor(Group, levels = rev(unique(myGroups))))
  }
  #
  #
  x1 <- myP %>% filter(Pvalue >= threshold)
  x2 <- myP %>% filter(Pvalue < threshold)
  #
  if(is.null(sug.threshold)) { x2 <- x2[0,] }
  #
  mp <- ggplot(x1, aes(x = Pos / 100000000)) + geom_blank(data = myG)
  #
  if(!is.null(vlines)) {
    myGM <- myG %>% filter(SNP %in% vlines) %>%
      mutate(SNP = factor(SNP, levels = vlines)) %>%
      dplyr::select(SNP, Chr, Pos) %>%
      arrange(SNP)
    #
    mp <- mp +
      geom_vline(data = myGM, alpha = 0.5,
                 aes(xintercept = Pos / 100000000, color = SNP)) +
      scale_color_manual(name = NULL, values = vline.colors) +
      scale_linetype_manual(name = NULL, values = vline.types)
  }
  #
  if(!is.null(hlines)) {
    mp <- mp + geom_hline(yintercept = hlines, alpha = 0.7)
  }
  #
  mp <- mp +
    geom_point(data = x2,
               size = 0.75, color = "black", alpha = 0.5,
               aes(y = Trait, shape = Model, fill = Model,
                   key1 = SNP, key2 = negLog10_P, key3 = negLog10_HBP)) +
    geom_point(size = 2.25, color = "black", alpha = 0.5,
               aes(y = Trait, shape = Model, fill = Model,
                   key1 = SNP, key2 = negLog10_P, key3 = negLog10_HBP))
  if(!is.null(groups)) {
    mp <- mp + facet_grid(Group ~ Chr, scales = "free", space = "free")
  } else { mp <- mp + facet_grid(. ~ Chr, scales = "free", space = "free") }
  mp <- mp +
    scale_fill_manual(name = NULL, values = model.colors, breaks = models) +
    scale_shape_manual(name = NULL, values = shapes, breaks = models) +
    scale_y_discrete(limits = rev) +
    scale_x_continuous(breaks = 0:20) +
    theme_gwaspr(legend.position = legend.position) +
    guides(shape = guide_legend(nrow = legend.rows, override.aes = list(size = 4)),
           color = guide_legend(nrow = legend.rows),
           fill = guide_legend(nrow = legend.rows)) +
    labs(title = title, y = NULL, x = "100 Mbp", caption = caption)
  #
  if(vline.legend == F) {
    mp <- mp + guides(color = vline.legend)
  }
  #
  mp
}

#folder = "GWAS_Results/"; traits = list_Traits(folder); groups = NULL
#threshold = round(-log10(0.00000005),1); sug.threshold = round(-log10(0.000005),1)
#chroms = NULL; pos1 = NULL; pos2 = NULL
#models =  c("MLM", "FarmCPU", "BLINK", "MLMM", "GLM", "CMLM", "SUPER")
#model.colors = c("darkgreen", "darkorange3", "steelblue", "darkred", "darkorchid4", "burlywood4", "darkseagreen4")
#shapes = 21:25; hlines = NULL; vlines = NULL; vline.colors = rep("red",length(vlines));
#vline.types = rep(1, length(vlines)); vline.legend = T; title = "Summary of Significant GWAS Results";
#caption = paste0("Significant Threshold = ", threshold, " = Large\nSuggestive Threshold = ", sug.threshold," = Small")
#rowread = 2000; legend.position = "bottom"; legend.rows = 1; plotHBPvalues = F; skyline = F
