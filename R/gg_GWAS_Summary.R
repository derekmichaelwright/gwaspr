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
                            models =  c("GLM","MLM","CMLM","MLMM","SUPER","FarmCPU","Blink"),
                            colors = c("darkgreen", "darkred", "darkorange3", "steelblue", "darkorchid4", "darkgoldenrod2"),
                            markers = NULL, hlines = NULL,
                            title = NULL,
                            caption = paste0("Threshold = ", threshold, "\nSuggestive = ", threshold2),
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
           Trait = factor(Trait, levels = rev(traits)))
  #
  x1 <- myP %>% filter(`-log10(p)` > threshold)
  x2 <- myP %>% filter(`-log10(p)` < threshold)
  #
  myG <- read.csv(paste0(folder, files[1])) %>%
    mutate(Model = myP$Model[1],
           Model = factor(Model, levels = models),
           Trait = myP$Trait[1],
           Trait = factor(Trait, levels = rev(traits)))
  #
  mp <- ggplot(x1, aes(x = Position / 100000000, y = Trait)) +
    geom_blank(data = myG)
  if(!is.null(markers)) {
    myGM <- myG %>% filter(SNP %in% markers)
    mp <- mp + geom_vline(data = myGM, alpha = 0.5, color = "red",
                          aes(xintercept = Position / 100000000))
  }
  if(!is.null(hlines)) {
    mp <- mp + geom_hline(yintercept = hlines, alpha = 0.7)
  }
  mp <- mp +
    geom_point(data = x2,
               size = 1, color = "black", alpha = 0.5,
               aes(shape = Model, fill = Model)) +
    geom_point(size = 2, color = "black", alpha = 0.5,
               aes(shape = Model, fill = Model)) +
    facet_grid(. ~ Chromosome, drop = F, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = colors) +
    scale_shape_manual(values = c(21:26)) +
    scale_y_discrete(drop = F) +
    theme_gwaspr(legend.position = "bottom") +
    guides(shape = guide_legend(override.aes = list(size = 4))) +
    labs(title = title, y = NULL, x = "100 Mbp", caption = caption)
  #
  mp
}

#setwd("C:/gitfolder/gwas_tutorial")
#library(gwaspr)
#folder <- "Results/"
#files <- list_Result_Files(folder)
#files <- files[!grepl("GLM", files)]
#traits <- unique(gsub("GAPIT.|GLM.|MLM.|CMLM.|MLMM.|SUPER.|FarmCPU.|Blink.|.GWAS.Results.csv", "", files))[c(1,3)]
#models <- c("MLM","MLMM","FarmCPU","Blink")
#threshold <- 6.7
#threshold2 <- 5
#markers <- "Lcu.2RBY.Chr6p12212845"
#title <- NULL
#hlines <- c(1.5,2.5)
#rowread <- 2000
#colors <- c("darkgreen", "darkred", "darkorange3", "steelblue", "darkorchid4", "darkgoldenrod2")
#caption = paste0("Threshold = ", threshold, "\nSuggestive = ", threshold2)


#setwd("C:/gitfolder/gwaspr")
#library(gwaspr)
#folder = "GWAS_Results/"
#traits = c("DTF_Nepal_2017", "Testa_Pattern")
#threshold = 7.3
#threshold2 = 6.7
