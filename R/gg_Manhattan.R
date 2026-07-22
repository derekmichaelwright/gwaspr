#' gg_Manhattan
#'
#' [Create manhattan plots from GAPIT GWAS results.](https://derekmichaelwright.github.io/gwaspr/articles/04_gg_Manhattan.html)
#' @param folder Folder containing GWAS results.
#' @param trait The trait to read.
#' @param title A title for the plot.
#' @param threshold Significant Threshold.
#' @param sug.threshold Suggested threshold.
#' @param chr Chromosomes to plot. Use if you want to plot a single chromosome.
#' @param markers Markers to be labelled.
#' @param labels Labels to be used for markers.
#' @param vlines Markers which will be used as a location for a vertical lines.
#' @param vline.colors colors for each vertical line.
#' @param vline.types lty for each vertical line.
#' @param legend Logical, whether or not to add a legend.
#' @param legend.rows Number of rows for the legend.
#' @param legend.box Alignment of the legend. Default is "horizontal", but it can be changed to "vertical".
#' @param point.sizes Sizes for the points. c("Not Sig", "Sig", "Sug").
#' @param pmax A max value for the y-axis. Markers with higher values will be lowered to pmax..
#' @param pmin A min Value for plotting. Markers with lower values will be removed.
#' @param models Models to read.
#' @param model.colors Colors for each model. Used if `facet = F`.
#' @param facet Logical, whether or not to produce a facetted or multi-model plot. Default is `facet = F`.
#' @param highlight.sig Logical, whether or not to highlight significant associations with a black circle. Used if `facet = F`.
#' @param sig.color Color for significant assoctiations. Used if `facet = T`.
#' @param chr.colors Colors for each chromosome. Used if `facet = T`.
#' @param chr.unit Unit for the x-axis. Can be one of c("kbp","100 kbp","Mbp","100 Mbp","Gbp").
#' @param plotHBPvalues Logical, if TRUE, H.B.P.values be uses.
#' @param skyline Which skyline type to use. Can be "NYC" or "Kansas". If NULL, it will use the highest P.value.
#' @param addQQ Logical, whether or not to add a QQ plot.
#' @return A manhattan plot.
#' @export

gg_Manhattan <- function(
    folder = "GWAS_Results/",
    trait = list_Traits(folder)[1],
    title = trait,
    threshold = NULL,
    sug.threshold = NULL,
    chr = NULL,
    markers = NULL,
    labels = markers,
    vlines = markers,
    vline.colors = rep("red", length(vlines)),
    vline.types = rep(1, length(vlines)),
    legend = T,
    legend.rows = 1,
    legend.box = "horizontal",
    point.sizes = c(0.3,1,0.75),
    pmax = NULL,
    pmin = 0,
    models =  c("MLM", "MLMM", "FarmCPU", "BLINK", "GLM", "CMLM", "SUPER"),
    model.colors = c("darkgreen", "darkred", "darkorange3", "steelblue", "darkorchid4", "blue2", "magenta3"),
    facet = F,
    highlight.sig = F,
    sig.color = "black",
    chr.colors = rep(c("darkgreen", "darkgoldenrod3"), 30),
    chr.unit = "100 Mbp",
    plotHBPvalues = F,
    skyline = NULL,
    addQQ = T
    ) {
  #
  # Read in files
  #
  fnames <- list_Result_Files(folder)
  fnames <- fnames[grepl(paste0(models, ".", trait, collapse="|"), fnames)]
  fnames <- fnames[grepl(paste0(trait, c(".csv","\\("), collapse="|"), fnames)]
  #
  if(!is.null(skyline)) {
    if(skyline == "NYC") { fnames <- fnames[!grepl("\\(Kansas\\)", fnames)] }
    if(skyline == "Kansas") {
      fnames <- fnames[!grepl("\\(NYC\\)&FarmCPU", fnames)]
      fnames <- fnames[!grepl("\\(NYC\\)&BLINK", fnames)]
    }
  }
  #
  xx <- NULL
  for(i in fnames) {
    mod <- substr(i, gregexpr("GWAS_Results", i)[[1]][1] + 13, gregexpr(".csv", i)[[1]][1] - 1)
    mod <- substr(mod, 1, gregexpr("\\.", mod)[[1]][1] - 1)
    sky <- substr(i, gregexpr("\\(", i)[[1]][1] + 1, gregexpr("\\)", i)[[1]][1] - 1)
    #
    xi <- read.csv(paste0(folder, i))
    if(sum(colnames(xi) == "nobs") > 0) { xi <- dplyr::select(xi, -nobs) }
    if(sum(colnames(xi) == "effect") > 0) { xi <- rename(xi, Effect=effect) }
    xi <- xi %>%
      mutate(Model = mod,
             negLog10_P = -log10(P.value),
             negLog10_Exp = -log10((rank(P.value, ties.method = "first") - 0.5)/nrow(.)),
             Type = sky )
    xx <- bind_rows(xx, xi)
  }
  #
  xx <- xx %>% filter(P.value > 0) %>%
    arrange(desc(P.value)) %>%
    filter(!duplicated(paste(SNP, Model)))
  #
  # Prep data
  #
  xx <- xx %>%
    mutate(Model = factor(Model, levels = models)) %>%
    arrange(desc(Model))
  #
  models <- unique(xx$Model)
  if(paste(models,collapse=",") != "MLM,MLMM,FarmCPU,BLINK,GLM,CMLM,SUPER" &
     paste(model.colors,collapse=",") == "darkgreen,darkred,darkorange3,steelblue,darkorchid4,burlywood4,darkseagreen4") {
    model.colors <- model.colors[c("MLM", "MLMM", "FarmCPU", "BLINK", "GLM", "CMLM", "SUPER") %in% models]
  }
  #
  if(is.null(threshold)) { threshold <- -log10(0.05/nrow(xi)) }
  #
  if(!is.null(pmax)) { xx <- xx %>% mutate(negLog10_P = ifelse(negLog10_P > pmax, pmax, negLog10_P)) }
  #
  xx <- xx %>%
    mutate(Sig.level = ifelse(negLog10_P >= threshold, "Sig","Not Sig"))
  if(!is.null(sug.threshold)) {
    xx <- xx %>%
      mutate(Sig.level = ifelse(negLog10_P < threshold & negLog10_P >= sug.threshold, "Sug", Sig.level))
  }
  #
  if(plotHBPvalues == T) {
    xx <- xx %>% mutate(Pvalue = -log10(H.B.P.Value))
  } else {
    xx <- xx %>% mutate(Pvalue = negLog10_P)
  }
  #
  x2 <- xx %>% filter(negLog10_P > threshold)
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
  if(!chr.unit %in% c("1 kbp", "10 kbp", "50 kbp", "100 kbp", "1 Mbp", "10 Mbp", "50 Mbp", "100 Mbp", "1 Gbp")) { print("error in chr.unit") }
  #
  if(!is.null(chr)) {
    xx <- xx %>% filter(Chr %in% chr)
    x2 <- x2 %>% filter(Chr %in% chr)
    #addQQ = F
  }
  #
  if(!is.null(vlines)) {
    vv <- xx %>%
      filter(SNP %in% vlines) %>%
      filter(!duplicated(SNP)) %>%
      mutate(SNP = factor(SNP, levels = vlines)) %>%
      dplyr::select(SNP, Chr, Pos)
  }
  #
  myBreaks <- 0:(round(max(xx$Pos)/x.unit))
  ylabel <- ifelse(plotHBPvalues == T, "-log<sub>10</sub>(*HBp*)", "-log<sub>10</sub>(*p*)")
  xx <- xx %>% filter(Pvalue >= pmin)
  #
  if(is.null(pmax)) { pmax <- max(xx$negLog10_P) }
  #
  # Start Plots
  #
  mp1 <- ggplot(xx, aes(x = Pos/x.unit, y = Pvalue)) +
    theme_gwaspr(axis.title.y = element_markdown()) +
    labs(title = title, y = ylabel, x = chr.unit)
  #
  mp2 <- ggplot(xx, aes(y = negLog10_P, x = negLog10_Exp)) +
    theme_gwaspr() +
    labs(y = NULL, x = "Expected")
  #
  # Add vlines
  #
  if(!is.null(vlines)) {
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
  if(!is.null(markers)) {
    xm <- xx %>%
      mutate(SNP = ifelse(SNP %in% markers, SNP, NA),
             Label = plyr::mapvalues(SNP, markers, labels)) %>%
      filter(!is.na(Label), negLog10_P > min(threshold, sug.threshold))
    mp1 <- mp1 +
      geom_text_repel(data = xm, aes(label = Label), size = 2)
  }
  #
  # Plot facetted by model
  #
  if(facet == T) {
    mp1 <- mp1 +
      geom_point(aes(fill = factor(Chr), size = Sig.level, key1 = SNP), pch = 21, color = alpha("white", 0)) +
      facet_grid(Model ~ paste("Chr", Chr), scales = "free", space = "free_x") +
      scale_fill_manual(values = chr.colors, guide = "none") +
      scale_size_manual(values = point.sizes, guide = "none") +
      scale_color_manual(name = "Marker", values = vline.colors) +
      scale_linetype_manual(name = "Marker", values = vline.types) +
      scale_x_continuous(breaks = myBreaks, minor_breaks = myBreaks) +
      guides(color = guide_legend(nrow = legend.rows, byrow = T)) +
      theme(legend.position = "bottom", legend.box=legend.box, legend.margin=margin())
    #
    if(highlight.sig == T) { mp1 <- mp1 + geom_point(data = x2, aes(key1 = SNP), pch = 21, size = point.sizes[2], fill = sig.color) }
    # legends
    if(legend == F) { mp1 <- mp1 + theme(legend.position = "none") }
    #
    if(addQQ == T) {
      mp2 <- mp2 +
        geom_point(pch = 16, size = point.sizes[2], color = chr.colors[1]) +
        geom_abline() +
        facet_grid(Model ~ "QQ", scales = "free")
      if(highlight.sig == T) { mp2 <- mp2 + geom_point(data = x2, aes(key1 = SNP), pch = 21, size = point.sizes[2], fill = sig.color) }
      #
      mp <- ggarrange(mp1, mp2, ncol = 2, widths = c(4,1), align = "h",
                      legend = "bottom", common.legend = T)
    } else { mp <- mp1 }
  } else {
    #
    # Plot models together
    #
    mp1 <- mp1 +
      geom_point(aes(fill = Model, size = Sig.level, key1 = SNP), pch = 21, color = alpha("white", 0)) +
      facet_grid(. ~ paste("Chr", Chr), scales = "free_x", space = "free_x") +
      scale_fill_manual(values = model.colors) +
      scale_size_manual(values = point.sizes, guide = "none") +
      scale_color_manual(name = "Marker", values = vline.colors) +
      scale_linetype_manual(name = "Marker", values = vline.types) +
      scale_x_continuous(breaks = myBreaks, minor_breaks = myBreaks) +
      scale_y_continuous(limits = c(pmin, (pmax+pmax*0.03))) +#, expand = c(0,0)
      guides(fill = guide_legend(nrow = legend.rows, override.aes = list(size = 2)),
             color = guide_legend(nrow = legend.rows, byrow = T, override.aes = list(alpha = 1))) +
      theme(legend.position = "bottom", legend.box=legend.box, legend.margin=margin())
    if(highlight.sig == T) { mp1 <- mp1 + geom_point(data = x2, aes(key1 = SNP), pch = 21, size = point.sizes[2], color = sig.color) }
    #
    if(legend == F) { mp1 <- mp1 + theme(legend.position = "none") }
    #
    if(addQQ == T) {
      mp2 <- mp2 +
        geom_point(aes(fill = Model), pch = 21, size = point.sizes[2], color = alpha("White",0)) +
        geom_abline() +
        facet_grid(. ~ "QQ") +
        scale_fill_manual(name = NULL, values = alpha(model.colors,0.8)) +
        scale_y_continuous(limits = c(pmin, (pmax+pmax*0.03))) +#, expand = c(0,0)
        guides(fill = guide_legend(nrow = legend.rows, override.aes = list(size = 2)),
               color = guide_legend(nrow = legend.rows, byrow = T, override.aes = list(alpha = 1)))
      if(highlight.sig == T) { mp2 <- mp2 + geom_point(data = x2, aes(key1 = SNP), pch = 21, size = point.sizes[2], color = sig.color) }
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

#folder = "GWAS_Results/"; trait = "DTF_Sask_2017"; title = trait; threshold = NULL; sug.threshold = NULL
#chr = NULL; markers = "Lcu.1GRN.Chr1p365986872"; labels = markers
#vlines = markers; vline.colors = rep("red", length(vlines)); vline.types = rep(1, length(vlines)); legend = T
#facet = F; addQQ = T; pmax = NULL; pmin = 0
#models = c("MLM", "FarmCPU", "BLINK", "MLMM", "GLM", "CMLM", "SUPER")#c("MLM")#
#models = "BLINK"
#model.colors = c("darkgreen", "darkorange3", "steelblue", "darkred", "darkorchid4", "burlywood4", "darkseagreen4")
#highlight.sig = F; sig.color = "darkred"; chr.colors = rep(c("darkgreen", "darkgoldenrod3"), 30)
#chr.unit = "100 Mbp"; legend.rows = 1; plotHBPvalues = F; skyline = NULL
#point.sizes = c(0.3,1.25,0.75)
