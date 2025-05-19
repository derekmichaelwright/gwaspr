#' gg_Manhattan
#'
#' Creates a manhattan plot.
#' @param folder Folder containing GWAS results.
#' @param trait The trait to read.
#' @param title A title for the plot.
#' @param threshold Significant Threshold.
#' @param sug.threshold Suggested threshold.
#' @param chrom Chromosomes to plot. Use if you want to plot a single chromosome.
#' @param pos1 starting position to plot.
#' @param pos2 ending position to plot.
#' @param vlines Markers which will be used as a location for a vertical lines.
#' @param vline.colors colors for each vertical line.
#' @param vline.types lty for each vertical line.
#' @param vline.legend Logical, whether or not to add a legend for the vlines.
#' @param markers Markers to be labelled.
#' @param labels Labels to be used for markers.
#' @param facet Logical, whether or not to produce a facetted or multi-model plot. Default is `facet = F`
#' @param addQQ Logical, whether or not to add a QQ plot
#' @param pmax A max value for the y-axis.
#' @param models Models to read.
#' @param model.colors Colors for each model. Used if `facet = F`.
#' @param highlight.sig Logical, whether or not to highlight significant associations with a black circle. Used if `facet = F`.
#' @param sig.col Color for significant assoctiations. Used if `facet = T`.
#' @param chrom.colors Colors for each chromosome. Used if `facet = T`.
#' @param chrom.unit Unit for the x-axis. Can be one of c("kbp","100 kbp","Mbp","100 Mbp","Gbp").
#' @param legend.rows Number of rows for the legend.
#' @return A manhattan plot.
#' @export

gg_Manhattan <- function (
    folder = "GWAS_Results/",
    trait = list_Traits(folder)[1],
    title = trait,
    threshold = NULL,
    sug.threshold = NULL,
    chrom = NULL,
    vlines = markers,
    vline.colors = rep("red", length(vlines)),
    vline.types = rep(1, length(vlines)),
    vline.legend = T,
    markers = NULL,
    labels = markers,
    facet = F,
    addQQ = T,
    pmax = NULL,
    models = c("MLM", "FarmCPU", "BLINK", "MLMM", "GLM", "CMLM", "SUPER"),
    model.colors = c("darkgreen", "darkorange3", "steelblue", "darkred", "darkorchid4", "burlywood4", "darkseagreen4"),
    highlight.sig = F,
    sig.col = "darkred",
    chrom.colors = rep(c("darkgreen", "darkgoldenrod3"), 30),
    chrom.unit = "100 Mbp",
    legend.rows = 1
    ) {
  #
  # Read in files
  #
  fnames <- list_Result_Files(folder)
  fnames <- fnames[grepl(paste(trait, collapse="|"), fnames)]
  fnames <- fnames[grepl(paste(models,collapse="|"), fnames)]
  #
  xx <- NULL
  for (i in fnames) {
    mod <- substr(i, gregexpr("GWAS_Results", i)[[1]][1] +
                    13, gregexpr(".csv", i)[[1]][1] - 1)
    mod <- substr(mod, 1, gregexpr("\\.", mod)[[1]][1] -
                    1)
    xi <- read.csv(paste0(folder, i))
    if (sum(colnames(xi) == "nobs") > 0) {
      xi <- select(xi, -nobs)
    }
    xi <- xi %>% mutate(Model = mod, `-log10(p)` = -log10(P.value),
                        `-log10(p)_exp` = -log10((rank(P.value, ties.method = "first") -
                                                    0.5)/nrow(.)))
    xx <- bind_rows(xx, xi)
  }
  #
  # Prep data
  #
  xx <- xx %>% filter(Model %in% models) %>%
    mutate(Model = factor(Model, levels = models)) %>%
    arrange(desc(Model))
  #
  if (is.null(threshold)) {
    threshold <- -log10(0.05/nrow(xi))
  }
  #
  if (!is.null(pmax)) {
    xx <- xx %>%
      mutate(`-log10(p)` = ifelse(`-log10(p)` > pmax, pmax, `-log10(p)`))
  }
  #
  xx <- xx %>%
    mutate(Sig.level = ifelse(`-log10(p)` >= threshold, "Sig","Not Sig"))
  if(!is.null(sug.threshold)) {
    xx <- xx %>%
      mutate(Sig.level = ifelse(`-log10(p)` < threshold & `-log10(p)` >= sug.threshold, "Sug", Sig.level))
  }
  x2 <- xx %>% filter(`-log10(p)` > threshold)
  #
  if(chrom.unit == "kbp")     { x.unit = 1000 }
  if(chrom.unit == "100 kbp") { x.unit = 100000 }
  if(chrom.unit == "Mbp")     { x.unit = 1000000 }
  if(chrom.unit == "100 Mbp") { x.unit = 100000000 }
  if(chrom.unit == "Gbp")     { x.unit = 1000000000 }
  if(!chrom.unit %in% c("kbp", "100 kbp", "Mbp", "100 Mbp", "Gbp")) { print("error in chrom.unit") }
  #
  if(!is.null(chrom)) {
    xx <- xx %>% filter(Chr %in% chrom)
    x2 <- x2 %>% filter(Chr %in% chrom)
  }
  #
  myBreaks <- 0:(round(max(xx$Pos)/x.unit))
  #
  # Start Plots
  #
  mp1 <- ggplot(xx, aes(x = Pos/x.unit, y = `-log10(p)`)) +
    theme_gwaspr(axis.title.y = element_markdown()) +
    labs(title = title, y = "-log<sub>10</sub>(*p*)", x = chrom.unit)
  #
  mp2 <- ggplot(xx, aes(y = `-log10(p)`, x = `-log10(p)_exp`)) +
    theme_gwaspr() +
    labs(title = "", y = NULL, x = "Expected")
  #
  # Add vlines
  #
  if (!is.null(vlines)) {
    vv <- xx %>% filter(SNP %in% vlines) %>% filter(!duplicated(SNP)) %>%
      mutate(SNP = factor(SNP, levels = vlines)) %>% select(-Model)
    mp1 <- mp1 +
      geom_vline(data = vv, aes(xintercept = Pos/x.unit, color = SNP, lty = SNP), alpha = 0.7)
  }
  #
  # Add threshold lines
  #
  mp1 <- mp1 +
    geom_hline(yintercept = threshold, color = "red", alpha = 0.8, linewidth = 0.5) +
    geom_hline(yintercept = sug.threshold, color = "blue", alpha = 0.8, linewidth = 0.5)
  mp2 <- mp2 +
    geom_hline(yintercept = threshold, color = "red", alpha = 0.8, linewidth = 0.5) +
    geom_hline(yintercept = sug.threshold, color = "blue", alpha = 0.8, linewidth = 0.5)
  #
  # Add Marker labels
  #
  if (!is.null(markers)) {
    xm <- xx %>%
      mutate(SNP = ifelse(SNP %in% markers, SNP, NA),
             Label = plyr::mapvalues(SNP, markers, labels)) %>%
      filter(!is.na(Label), `-log10(p)` > min(threshold, sug.threshold))
    mp1 <- mp1 +
      geom_text_repel(data = xm, aes(label = Label), size = 2)
  }
  # vline legends
  if (vline.legend == T) {
    mp1 <- mp1 +
      scale_color_manual(name = NULL, values = vline.colors) +
      scale_linetype_manual(name = NULL, values = vline.types)
  } else {
    mp1 <- mp1 +
      scale_color_manual(name = NULL, values = vline.colors, guide = "none") +
      scale_linetype_manual(name = NULL, values = vline.types, guide = "none")
  }
  #
  # Plot facetted by model
  #
  if (facet == T) {
    mp1 <- mp1 +
      geom_point(aes(fill = factor(Chr), size = Sig.level), pch = 21, color = alpha("white", 0)) +
      geom_point(data = x2, pch = 21, size = 1.5, color = "black", fill = sig.col, alpha = 0.8) +
      facet_grid(Model ~ Chr, scales = "free", space = "free_x") +
      scale_fill_manual(name = NULL, values = alpha(chrom.colors, 0.8), guide = "none") +
      scale_size_manual(name = NULL, values = c(0.4,1.25,0.75), guide = "none") +
      scale_x_continuous(breaks = myBreaks, minor_breaks = myBreaks) +
      guides(color = guide_legend(nrow = legend.rows, byrow = T, override.aes = list(alpha = 1))) +
      theme(legend.position = "bottom")
    #
    if(addQQ == T) {
      mp2 <- mp2 +
        geom_point(pch = 1, color = chrom.colors[1], alpha = 0.8) +
        geom_point(data = x2, pch = 21, color = "black", fill = "darkred", alpha = 0.8) +
        geom_abline() +
        facet_grid(Model ~ "QQ", scales = "free_y")
      #
      mp <- ggarrange(mp1, mp2, ncol = 2, widths = c(4,1), align = "h",
                      legend = "bottom", common.legend = T)
    } else { mp <- mp1 }
  } else {
    #
    # Plot models together
    #
    mp1 <- mp1 +
      geom_point(aes(fill = Model, size = Sig.level), pch = 21, color = alpha("white", 0)) +
      facet_grid(. ~ Chr, scales = "free", space = "free_x") +
      scale_fill_manual(name = NULL, values = alpha(model.colors,0.8)) +
      scale_size_manual(name = NULL, values = c(0.3,1.25,0.75), guide = "none") +
      scale_x_continuous(breaks = myBreaks, minor_breaks = myBreaks) +
      guides(fill = guide_legend(nrow = legend.rows, override.aes = list(size = 2)),
             color = guide_legend(nrow = legend.rows, byrow = T, override.aes = list(alpha = 1))) +
      theme(legend.position = "bottom")
    if(highlight.sig == T) {
      mp1 <- mp1 + geom_point(data = x2, fill = alpha("white", 0), pch = 21, size = 1.25)
    }
    #
    if(addQQ == T) {
      mp2 <- mp2 +
        geom_point(pch = 1, aes(color = Model)) +
        geom_point(data = x2, aes(color = Model)) +
        geom_abline() +
        facet_grid(. ~ "QQ", scales = "free_y") +
        scale_color_manual(name = NULL, values = alpha(model.colors,0.8)) +
        guides(fill = guide_legend(nrow = legend.rows, override.aes = list(size = 2)),
               color = guide_legend(nrow = legend.rows, byrow = T))
      if(highlight.sig == T) {
        mp2 <- mp2 + geom_point(data = x2, fill = alpha("white", 0), pch = 21, size = 1.25)
      }
      #
      mp <- ggarrange(mp1, mp2, ncol = 2, widths = c(4,1), align = "h",
                      legend = "bottom", common.legend = T)
    } else { mp <- mp1 }
  }
  #
  # Output Plot
  #
  mp
}
