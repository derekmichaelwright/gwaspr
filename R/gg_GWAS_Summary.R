#' gg_GWAS_Summary
#'
#' Creates a summary GWAS plot of significant associations.
#' @param folder Folder containing GWAS results.
#' @param traits The traits to read.
#' @param threshold Significant threshold.
#' @param threshold2 Suggestive threshold.
#' @param models Models to read.
#' @param colors Colors for each model.
#' @param markers Markers to be labelled with a vertical red line.
#' @param hlines Locations for horizontal lines. e.g., hlines = c(1.5,2.5).
#' @param title A title for the plot.
#' @param caption A caption for the plot.
#' @param rowread Number of rows to read for each GWAS results file.
#' @return A GWAS summary plot.
#' @export

gg_GWAS_Summary <- function(folder = NULL, traits = list_Traits(),
                            threshold = -log10(0.00000005),
                            threshold2 = -log10(0.000005),
                            models =  c("MLM","MLMM","FarmCPU","Blink","GLM"),
                            colors = c("darkgreen", "darkred", "darkorange3", "steelblue", "darkgoldenrod2"),
                            shapes = 21:25, hlines = NULL,
                            markers = NULL, markers2 = NULL, markers3 = NULL,
                            title = NULL,
                            caption = paste0("Sig Threshold = ", threshold, " = Large\nSuggestive = ", threshold2," = Small"),
                            rowread = 2000
                            ) {
  #
  files <- list_Result_Files(folder)
  files <- files[grepl(paste(traits,collapse="|"), files)]
  files <- files[grepl(paste(models,collapse="|"), files)]
  #
  myP <- NULL
  #
  #i<-files[1]
  for(i in files) {
    myPi <- table_GWAS_Results(folder = folder, file = i,
              threshold = threshold, threshold2 = threshold2)
    if(nrow(myPi)>0) { myP <- bind_rows(myP, myPi) }
  }
  #
  myP <- myP %>% filter(!is.na(SNP)) %>%
    arrange(Chromosome, Position, P.value, Trait) %>%
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
  #
  mp <- ggplot(x1, aes(x = Position / 100000000, y = Trait, key1 = Chromosome, key2 = Position, key3 = P.value)) +
    geom_blank(data = myG)
  if(!is.null(markers)) {
    myGM <- myG %>% filter(SNP %in% markers)
    mp <- mp + geom_vline(data = myGM, alpha = 0.5, color = "red",
                          aes(xintercept = Position / 100000000))
  }
  if(!is.null(markers2)) {
    myGM <- myG %>% filter(SNP %in% markers2)
    mp <- mp + geom_vline(data = myGM, alpha = 0.5, color = "blue",
                          aes(xintercept = Position / 100000000))
  }
  if(!is.null(markers3)) {
    myGM <- myG %>% filter(SNP %in% markers3)
    mp <- mp + geom_vline(data = myGM, alpha = 0.5, color = "green",
                          aes(xintercept = Position / 100000000))
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
    facet_grid(. ~ Chromosome, drop = F, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = colors, breaks = models) +
    scale_shape_manual(values = shapes, breaks = models) +
    scale_x_continuous(breaks = 0:20) +
    scale_y_discrete(drop = F) +
    theme_gwaspr(legend.position = "bottom") +
    guides(shape = guide_legend(override.aes = list(size = 4))) +
    labs(title = title, y = NULL, x = "100 Mbp", caption = caption)
  #
  mp
}
