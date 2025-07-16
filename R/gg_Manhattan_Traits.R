#' gg_Manhattan_Traits
#'
#' Creates a manhattan plot.
#' @param folder Folder containing GWAS results.
#' @param traits The traits to read.
#' @param title A title for the plot.
#' @param threshold Significant Threshold.
#' @param sug.threshold Suggested threshold.
#' @param chrom Chromosomes to plot. Use if you want to plot a single chromosome.
#' @param markers Markers to be labelled.
#' @param labels Labels to be used for markers.
#' @param vlines Markers which will be used as a location for a vertical lines.
#' @param vline.colors colors for each vertical line.
#' @param vline.types lty for each vertical line.
#' @param vline.legend Logical, whether or not to add a legend for the vlines.
#' @param addQQ Logical, whether or not to add a QQ plot.
#' @param pmax A max value for the y-axis.
#' @param models Model to read.
#' @param trait.colors Colors for each trait.
#' @param chrom.unit Unit for the x-axis. Can be one of c("kbp","100 kbp","Mbp","100 Mbp","Gbp").
#' @param legend.rows Number of rows for the legend.
#' @param skyline Which skyline type to use. Can be "NYC" or "Kansas". If left NULL, it will use the highest P.value.
#' @return A manhattan plot.
#' @export

gg_Manhattan_Traits <- function (
    folder = "GWAS_Results/",
    traits = list_Traits(folder)[1:2],
    title = NULL,
    threshold = NULL,
    sug.threshold = NULL,
    chrom = NULL,
    markers = NULL,
    labels = markers,
    vlines = markers,
    vline.colors = rep("red", length(vlines)),
    vline.types = rep(1, length(vlines)),
    vline.legend = F,
    addQQ = T,
    pmax = NULL,
    models = c("MLM","MLMM","FarmCPU","BLINK"),
    trait.colors = c("darkgreen", "darkorange3", "steelblue", "darkred", "darkorchid4", "burlywood4", "darkseagreen4"),
    chrom.unit = "100 Mbp",
    legend.rows = 1,
    plotHBPvalues = F,
    skyline = NULL
    ) {
  #
  # Read in files
  #
  fnames <- list_Result_Files(folder)
  fnames <- fnames[grepl(paste(traits, collapse="|"), fnames)]
  fnames <- fnames[grepl(paste(models, collapse="|"), fnames)]
  #if(removeKansas == T) { fnames <- fnames[!grepl("Kansas", fnames)] }
  #
  xx <- NULL
  for(i in fnames) {
    trt <- substr(i, gregexpr("GWAS_Results", i)[[1]][1]+13,
                    gregexpr(".csv", i)[[1]][1]-1 )
    trt <- gsub("\\(Kansas\\)|\\(NYC\\)", "", trt)
    #mod <- substr(trt, gregexpr("GWAS_Results", i)[[1]][1] + 13, gregexpr(".csv", i)[[1]][1] - 1)
    mod <- substr(trt, 1, gregexpr("\\.", trt)[[1]][1] - 1)
    trt <- substr(trt, gregexpr("\\.", trt)[[1]][1]+1, nchar(trt))
    sky <- substr(i, gregexpr("\\(", i)[[1]][1] + 1, gregexpr("\\)", i)[[1]][1] - 1)
    #
    xi <- read.csv(paste0(folder, i))
    if(sum(colnames(xi) == "nobs") > 0) { xi <- select(xi, -nobs) }
    if(sum(colnames(xi) == "effect") > 0) { xi <- rename(xi, Effect=effect) }
    xi <- xi %>%
      mutate(Model = mod, Trait = trt,
             negLog10_P = -log10(P.value),
             negLog10_Exp = -log10((rank(P.value, ties.method = "first") - 0.5)/nrow(.)),
             Type = sky )
    xx <- bind_rows(xx, xi)
  }
  #
  #if(removeKansas == F) { xx <- xx %>% arrange(desc(P.value)) %>% filter(!duplicated(paste(SNP, Model, P.value))) }
  if(!is.null(skyline)) {
    if(skyline == "NYC")    { xx <- xx %>% filter(!paste(Model, Type) %in% c("FarmCPU Kansas", "BLINK Kansas")) }
    if(skyline == "Kansas") { xx <- xx %>% filter(!paste(Model, Type) %in% c("FarmCPU NYC", "BLINK NYC")) }
  }
  #
  xx <- xx %>% arrange(desc(P.value)) %>% filter(!duplicated(paste(SNP, Model, P.value)))
  #
  # Prep data
  #
  xx <- xx %>% filter(Model %in% models) %>%
    mutate(Model = factor(Model, levels = models),
           Trait = factor(Trait, levels = traits)) %>%
    arrange(desc(Model), desc(Trait))
  #
  if(is.null(threshold)) {
    threshold <- -log10(0.05/nrow(xi))
  }
  #
  if(!is.null(pmax)) {
    xx <- xx %>%
      mutate(negLog10_P = ifelse(negLog10_P > pmax, pmax, negLog10_P))
  }
  #
  xx <- xx %>%
    mutate(Sig.level = ifelse(negLog10_P >= threshold, "Sig","NA"))
  if(!is.null(sug.threshold)) {
    xx <- xx %>%
      mutate(Sig.level = ifelse(negLog10_P < threshold & negLog10_P >= sug.threshold, "Sug", Sig.level))
  }
  xx <- xx %>% mutate(Sig.level = factor(Sig.level, levels = c("Sig","Sug","NA")))
  #
  if(plotHBPvalues == T) {
    xx <- xx %>% mutate(Pvalue = -log10(H.B.P.Value))
  } else {
    xx <- xx %>% mutate(Pvalue = negLog10_P)
  }
  #
  x2 <- xx %>% filter(negLog10_P > threshold)
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
  ylabel <- ifelse(plotHBPvalues == T, "-log<sub>10</sub>(*HBp*)", "-log<sub>10</sub>(*p*)")
  #
  # Start Plots
  #
  mp1 <- ggplot(xx, aes(x = Pos/x.unit, y = Pvalue)) +
    theme_gwaspr(axis.title.y = element_markdown()) +
    labs(title = title, y = ylabel, x = chrom.unit)
  #
  mp2 <- ggplot(xx, aes(y = negLog10_P, x = negLog10_Exp)) +
    theme_gwaspr() +
    labs(title = "", y = NULL, x = "Expected")
  #
  # Add vlines
  #
  if(!is.null(vlines)) {
    vv <- xx %>%
      filter(SNP %in% vlines) %>%
      filter(!duplicated(SNP)) %>%
      mutate(SNP = factor(SNP, levels = vlines)) %>%
      select(SNP, Chr, Pos)
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
  # vline legends
  if(vline.legend == T) {
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
  mp1 <- mp1 +
    geom_point(aes(fill = Trait, size = Sig.level), pch = 21, color = alpha("white", 0)) +
    facet_grid(Model ~ Chr, scales = "free", space = "free_x") +
    scale_fill_manual(name = NULL, values = alpha(trait.colors, 0.7)) +
    scale_size_manual(name = NULL, values = c(1.25,0.75,0.4), guide = "none") +
    scale_x_continuous(breaks = myBreaks, minor_breaks = myBreaks) +
    guides(color = guide_legend(nrow = legend.rows, byrow = T, override.aes = list(alpha = 1))) +
    theme(legend.position = "bottom")
    #
  if(addQQ == T) {
    mp2 <- mp2 +
      geom_point(aes(fill = Trait), pch = 1, size = 1.25) +
      geom_abline() +
      facet_grid(Model ~ "QQ", scales = "free_y") +
      scale_color_manual(name = NULL, values = alpha(trait.colors, 0.7)) +
      theme(legend.position = "none")
    #
    mp <- ggarrange(mp1, mp2, ncol = 2, widths = c(4,1), align = "h",
                    legend = "bottom", common.legend = T)
  } else { mp <- mp1 }
  #
  # Output Plot
  #
  mp
}

#folder = "GWAS_Results/"; traits = list_Traits(folder)[3:4];
#title = NULL; threshold = 6.8; sug.threshold = 5;
#chrom = NULL;
#markers = "Lcu.1GRN.Chr1p352153929"; labels = "352"; vlines = markers;
#vline.colors = rep("red", length(vlines)); vline.types = rep(1, length(vlines)); vline.legend = F;
#addQQ = T; pmax = NULL; models = "MLM";
#highlight.sig = F; sig.col = "darkred";
#trait.colors = c("darkgreen", "darkorange3", "steelblue", "darkred", "darkorchid4", "burlywood4", "darkseagreen4")
#chrom.unit = "100 Mbp"; legend.rows = 1; plotHBPvalues = F; skyline = NULL
