#' gg_GWAS_Summary
#'
#' Creates a summary GWAS plot of significant associations.
#' @param folder Folder containing GWAS results.
#' @param traits The traits to read.
#' @param threshold Significant threshold.
#' @param sug.threshold Suggestive threshold.
#' @param models Models to read.
#' @param colors Colors for each model.
#' @param shapes The shape values to use for the different models. e.g., 21:25
#' @param vlines Markers to be labelled with a vertical red line.
#' @param vline.colors colors for each vertical line.
#' @param vline.legend Logical, display of vline color legend.
#' @param hlines Locations for horizontal lines. e.g., hlines = c(1.5,2.5).
#' @param title A title for the plot.
#' @param caption A caption for the plot.
#' @param rowread Number of rows to read for each GWAS results file.
#' @param lrows Number of rows for the legend.
#' @return A GWAS summary plot.
#' @export

gg_GWAS_Summary <- function(folder = NULL, traits = list_Traits(),
                            threshold = -log10(0.00000005),
                            sug.threshold = -log10(0.000005),
                            models =  c("MLM", "MLMM", "FarmCPU", "BLINK", "GLM"),
                            colors = c("darkgreen", "darkred", "darkorange3", "steelblue", "darkgoldenrod2"),
                            shapes = 21:25,
                            hlines = NULL,
                            vlines = NULL,
                            vline.colors = rep("red",length(vlines)),
                            vline.legend = T,
                            title = NULL,
                            caption = paste0("Sig Threshold = ", threshold, " = Large\nSuggestive = ", sug.threshold," = Small"),
                            rowread = 2000,
                            legend.position = "bottom",
                            lrows = 1 ) {
  #
  files <- list_Result_Files(folder)
  files <- files[grepl(paste(traits,collapse="|"), files)]
  files <- files[grepl(paste(models,collapse="|"), files)]
  #
  myP <- NULL
  #
  #i<-files[20]
  for(i in files) {
    myPi <- table_GWAS_Results(folder = folder, files = i,
              threshold = threshold, sug.threshold = sug.threshold)
    if(nrow(myPi)>0) { myP <- bind_rows(myP, myPi) }
  }
  #
  myP <- myP %>% filter(!is.na(SNP)) %>%
    arrange(Chr, Pos, P.value, Trait) %>%
    mutate(Model = factor(Model, levels = models),
           Trait = factor(Trait, levels = rev(traits))) %>%
    filter(!is.na(Trait)) %>%
    arrange(desc(Model))
  #
  x1 <- myP %>% filter(`-log10(p)` > threshold)
  x2 <- myP %>% filter(`-log10(p)` < threshold)
  #
  myG <- read.csv(paste0(folder, files[1])) %>%
    mutate(Trait = myP$Trait[1],
           Trait = factor(Trait, levels = rev(traits)))
  #colnames(myG)[c(2:3,5)] <- c("Chromosome","Position","MAF")
  #
  mp <- ggplot(x1, aes(x = Pos / 100000000, y = Trait)) +
    geom_blank(data = myG)
  if(!is.null(vlines)) {
    myGM <- myG %>% filter(SNP %in% vlines) %>%
      mutate(SNP = factor(SNP, levels = vlines)) %>%
      arrange(SNP)
    #for(i in 1:nrow(myGM)) {
    #  mp <- mp + geom_vline(data = myGM[i,], alpha = 0.5, color = vline.colors[i],
    #                        aes(xintercept = Position / 100000000))
    #}
    mp <- mp +
      geom_vline(data = myGM, alpha = 0.5,
                 aes(xintercept = Pos / 100000000, color = SNP)) +
      scale_color_manual(name = NULL, values = vline.colors)
  }
  if(!is.null(hlines)) {
    mp <- mp + geom_hline(yintercept = hlines, alpha = 0.7)
  }
  mp <- mp +
    geom_point(data = x2,
               size = 0.75, color = "black", alpha = 0.5,
               aes(shape = Model, fill = Model, key1 = SNP, key2 = `-log10(p)`)) +
    geom_point(size = 2.25, color = "black", alpha = 0.5,
               aes(shape = Model, fill = Model, key1 = SNP, key2 = `-log10(p)`)) +
    #geom_point(color = "black", alpha = 0.5,
    #           aes(shape = Model, fill = Model, size = `-log10(p)`)) +
    #scale_size_continuous(range = c(0.75, 3), guide = "none") +
    facet_grid(. ~ Chr, drop = F, scales = "free_x", space = "free_x") +
    scale_fill_manual(name = NULL, values = colors, breaks = models) +
    scale_shape_manual(name = NULL, values = shapes, breaks = models) +
    scale_x_continuous(breaks = 0:20) +
    scale_y_discrete(drop = F) +
    theme_gwaspr(legend.position = legend.position#,
                 #legend.box="vertical",
                 #legend.spacing.y = unit(0.4, "cm"),
                 #legend.margin = margin(-0.2,0,0,0, unit="cm")
                 ) +
    guides(shape = guide_legend(nrow = lrows, override.aes = list(size = 4)),
           color = guide_legend(nrow = lrows),
           fill = guide_legend(nrow = lrows)) +
    labs(title = title, y = NULL, x = "100 Mbp", caption = caption)
  if(vline.legend == F) {
    mp <- mp + guides(color = vline.legend)
  }
  #
  mp
}
