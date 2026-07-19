#' gg_Manhattan_Zoom
#'
#' [Create manhattan plots from GAPIT GWAS results zoomed into a specific region on a chromosome.](https://derekmichaelwright.github.io/gwaspr/articles/05_gg_Manhattan_Zoom.html)
#' @param folder Folder containing GWAS results.
#' @param trait The trait to read.
#' @param title A title for the plot.
#' @param chr Chromosome to plot.
#' @param pos1 Start position on chromosome.
#' @param pos2 End position on chromosome.
#' @param threshold Significant Threshold.
#' @param sug.threshold Suggested threshold.
#' @param markers Markers to be labelled.
#' @param labels Labels to be used for markers.
#' @param vlines Markers which will be used as a location for a vertical lines.
#' @param vline.colors colors for each vertical line.
#' @param vline.types lty for each vertical line.
#' @param vline.legend Logical, whether or not to add a legend for the vlines.
#' @param facet Logical, whether or not to produce a facetted or multi-model plot. Default is `facet = F`.
#' @param pmax A max value for the y-axis.
#' @param models Models to read.
#' @param model.colors Colors for each model. Used if `facet = F`.
#' @param sig.color Color for significant assoctiations.
#' @param legend.rows Number of rows for the legend.
#' @param legend.box Alignment of the legend. Default is "horizontal", but it can be changed to "vertical".
#' @param point.sizes Sizes for the points. c("Not Sig", "Sig", "Sug").
#' @param plotHBPvalues Logical, should H.B.P.Values be uses.
#' @param skyline Which skyline type to use. Can be "NYC" or "Kansas". If left NULL, it will use the highest P.value.
#' @return A manhattan plot.
#' @export

gg_Manhattan_Zoom <- function(
    folder = "GWAS_Results/",
    trait = list_Traits(folder)[1],
    title = trait,
    chr = 1, pos1 = NULL, pos2 = NULL,
    addGenome = T,
    threshold = NULL,
    sug.threshold = NULL,
    markers = NULL,
    labels = markers,
    vlines = markers,
    vline.colors = rep("red", length(vlines)),
    vline.types = rep(1, length(vlines)),
    vline.legend = T,
    facet = F,
    pmax = NULL,
    models =  c("MLM", "MLMM", "FarmCPU", "BLINK",  "GLM", "CMLM", "SUPER"),
    model.colors = c("darkgreen", "darkred", "darkorange3", "steelblue", "darkorchid4", "blue2", "magenta3"),
    sig.color = "black",
    legend.rows = 1,
    legend.box = "horizontal",
    point.sizes = c(0.3,1,0.75),
    plotHBPvalues = F,
    skyline = "Kansas"
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
  #
  for(i in fnames) {
    mod <- substr(i, gregexpr("GWAS_Results", i)[[1]][1]+13, gregexpr(".csv", i)[[1]][1]-1)
    mod <- substr(mod, 1, gregexpr("\\.", mod)[[1]][1]-1)
    sky <- substr(i, gregexpr("\\(", i)[[1]][1] + 1, gregexpr("\\)", i)[[1]][1] - 1)
    xi <- read.csv(paste0(folder, i))
    if(sum(colnames(xi)=="nobs")>0) { xi <- dplyr::select(xi, -nobs) }
    xi <- xi %>% mutate(Model = mod, Type = sky, negLog10_P = -log10(P.value))
    xx <- bind_rows(xx, xi)
  }
  #
  xx <- xx %>% filter(P.value > 0) %>%
    arrange(desc(P.value)) %>%
    filter(!duplicated(paste(SNP, Model)))
  #
  # Prep data
  #
  #
  if(is.null(pos1)) { pos1 <- 0 }
  if(is.null(pos2)) { pos2 <- max(xx %>% filter(Chr == chr) %>% pull(Pos)) }
  myChroms <- xx %>% group_by(Chr) %>%
    summarise(Chr_max = max(Pos)) %>%
    ungroup() %>%
    mutate(Reg_min = ifelse(Chr == chr, 0, NA),
           Reg_max = ifelse(Chr == chr, Chr_max, NA),
           pos1 = ifelse(Chr == chr, pos1, NA),
           pos2 = ifelse(Chr == chr, pos2, NA))
  #
  xx <- xx %>%
    mutate(Model = factor(Model, levels = models)) %>%
    filter(Chr == chr, Pos > pos1, Pos < pos2, !is.na(Model)) %>%
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
  x2 <- xx %>% filter(negLog10_P >= threshold)
  x3 <- xx %>% filter(SNP %in% markers) %>%
    mutate(Label = plyr::mapvalues(SNP, markers, labels))
  x4 <- xx %>% filter(SNP %in% vlines) %>%
    filter(!duplicated(paste(SNP, Model)))
  #
  # Start Plots
  #
  mp1 <- ggplot(myChroms) +
    #chromosomes
    geom_rect(aes(xmax = Chr_max), xmin = 0, ymin = 0, ymax = 1,
              fill = "white", color = "black") +
    #
    #chr
    geom_rect(aes(xmin = Reg_min, xmax = Reg_max), ymin = 0, ymax = 1,
              fill = "black", alpha = 0.2) +
    #region
    geom_rect(aes(xmin = pos1, xmax = pos2), ymin = 0, ymax = 1,
              fill = "black", alpha = 0.2, color = "black") +
    facet_grid(. ~ paste("Chr", Chr), space = "free_x", scales = "free_x") +
    theme_void() + labs(subtitle = "Genome")
  mp2 <- ggplot(myChroms %>% filter(!is.na(pos1))) +
    #chr
    geom_rect(aes(xmin = Reg_min, xmax = Reg_max), ymin = 0, ymax = 1,
              fill = "black", alpha = 0.2, color = "black") +
    #region
    geom_rect(aes(xmin = pos1, xmax = pos2), ymin = 0, ymax = 1,
              fill = "black", alpha = 0.2, color = "black") +
    theme_void() + labs(subtitle = paste("Chromsome",chr))
  mp <- ggplot(xx, aes(x = Pos / 1000000, y = Pvalue))
  #
  # Add vlines
  #
  if(!is.null(vlines)) {
    mp <- mp +
      geom_vline(data = x4, aes(xintercept = Pos / 1000000, color = SNP), alpha = 0.5) +
      scale_color_manual(name = "Marker", values = vline.colors)
  }
  #
  # Add threshold lines
  #
  if(is.null(pmax)) { pmax <- max(xx$negLog10_P) }
  mp <- mp +
    geom_hline(yintercept = threshold, color = "red", alpha = 0.8, linewidth = 0.5) +
    geom_hline(yintercept = sug.threshold, color = "blue", alpha = 0.8, linewidth = 0.5)
  #
  # Plot the Rest
  #
  mp <- mp +
    geom_point(aes(fill = Model, size = Sig.level, key1 = SNP), pch = 21, color = alpha("white",0)) +
    scale_size_manual(values = point.sizes, guide = "none") +
    guides(color = guide_legend(nrow = legend.rows, byrow = T, override.aes = list(alpha = 1))) +
    theme_gwaspr(legend.position = "bottom", legend.box=legend.box, legend.margin=margin(),
                 axis.title.y = element_markdown()) +
    labs(title = trait, y = "-log<sub>10</sub>(*p*)", x = "Mbp")
  #
  # Add Marker labels
  #
  if(!is.null(markers)) {
    mp <- mp + geom_text_repel(data = x3, aes(label = Label), size = 2)
  }
  #
  # Facet plot
  #
  if(facet == T) {
    mp <- mp +
      facet_grid(Model ~ paste("Chromosome", Chr), scales = "free") +
      scale_fill_manual(values = alpha(model.colors,0.8), guide = "none")
    #if(nrow(x3)>0) { mp <- mp + geom_point(data = x3, aes(fill = Model, size = Sig.level), color = sig.color, pch = 21) }
    if(addGenome == T) { mp <- ggarrange(mp, mp2, mp1, ncol = 1, nrow = 3, heights = c(3.5*length(models),1,1.25)) }
  } else {
    mp <- mp +
      facet_grid(. ~ paste("Chromosome", Chr)) +
      scale_fill_manual(values = alpha(model.colors,0.8)) +
      guides(fill = guide_legend(nrow = legend.rows, override.aes = list(size = 2)))
    #if(nrow(x3)>0) { mp <- mp + geom_point(data = x3, aes(fill = Model, size = Sig.level), color = sig.color, pch = 21) }
    if(addGenome == T) { mp <- ggarrange(mp, mp2, mp1, ncol = 1, nrow = 3, heights = c(6.75,1,1.25)) }
  }
  mp
}

#folder = "vignettes/GWAS_Results/"; trait = "Cotyledon_RedvsYellow"; title = trait
#chr = 1; pos1 = 350000000; pos2 = 400000000;
#threshold = NULL; sug.threshold = NULL
#markers = c("Lcu.1GRN.Chr1p365986872", "Lcu.1GRN.Chr1p361840399"); labels = c("365Mbp","361Mbp"); vlines = markers
#vline.colors = rep("red", length(vlines)); vline.types = rep(1, length(vlines))
#vline.legend = T; models = c("MLM", "FarmCPU", "BLINK")
#model.colors = c("darkgreen", "darkorange3", "steelblue", "darkred", "darkorchid4", "burlywood4", "darkseagreen4")
#facet = T; sig.color = "red"; legend.rows = 1
#plotHBPvalues = F; pmax = NULL;legend.box = "horizontal";point.sizes = c(0.3,1,0.75);skyline = "Kansas"
#addGenome = T
