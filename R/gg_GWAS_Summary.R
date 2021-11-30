#' gg_GWAS_Summary
#'
#' Creates a summary GWAS plot of significant associations.
#' @param folder Folder containing GWAS results.
#' @param trait The trait to read.
#' @param subtitle A subtitle for the plot.
#' @param markers Markers to be labelled.
#' @param labels Labels to be used for markers.
#' @param lines Logical value of whether or not to include vertical lines with markers.
#' @param models Models to read.
#' @param colors1 Colors for each chromosome
#' @param colors2 Colors for each model
#' @return A manhattan plot.
#' @export

#setwd("C:/gitfolder/gwas_tutorial")
#library(gwaspr)
#library(tidyverse)
#folder <- "Results/"
#files <- list_Results(folder)
#files <- files[!grepl("GLM", files)]
#threshold <- 6.7
#threshold2 <- 5
#g.range<- 3000000
#rowread <- 2000

gg_GWAS_Summary <- function(folder, files,
                            threshold, threshold2 = NULL,
                            markers = NULL, hlines = F,
                            g.range = 3000000,
                            rowread = 2000,
                            plotly = F,
                            plotly_filename = "GWAS_Summary.html",
                            caption = paste0("Threshold = ", threshold,
                                             "\nSuggestive = ", threshold2)) {
  #
  myP <- NULL
  #
  #i<-files[1]
  for(i in files) {
    myP <- bind_rows(myP, list_GWAS_Results(folder = folder, file = i,
                                            threshold = threshold, threshold2 = threshold2)) #%>% dropNAcol())
  }
  #
  myModels <- c("GLM","MLM","CMLM","MLMM","SUPER","FarmCPU","Blink")
  myColors <- c("darkgreen", "darkred", "darkorange3", "steelblue", "darkorchid4", "darkgoldenrod2")
  myP <- myP %>% filter(!is.na(SNP)) %>%
    arrange(Chromosome, Position, P.value, Trait) %>%
    mutate(Model = factor(Model, levels = myModels))
  #
  x1 <- myP %>% filter(`-log10(p)` > threshold)
  x2 <- myP %>% filter(`-log10(p)` < threshold)
  #
  myG <- read.csv(paste0(folder, files[1])) %>% mutate(Trait = myP$Trait[1])
  mp <- ggplot(x1, aes(x = Position / 100000000, y = Trait)) +
    geom_blank(data = myG) +
    #geom_vline(data = myGM, alpha = 0.5, color = "red",
    #           aes(xintercept = Position / 100000000)) +
    geom_point(data = x2,
               size = 1, color = "black", alpha = 0.5,
               aes(shape = Model, fill = Model)) +
    geom_point(size = 2, color = "black", alpha = 0.5,
               aes(shape = Model, fill = Model)) +
    #geom_hline(yintercept = hlines, alpha = 0.7) +
    facet_grid(. ~ Chromosome, drop = F, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = myColors) +
    scale_shape_manual(values = c(21:26)) +
    #scale_y_discrete(drop = F) +
    theme_gwaspr(legend.position = "bottom") +
    guides(shape = guide_legend(override.aes = list(size = 4))) +
    labs(y = NULL, x = "100 Mbp", caption = caption)
  #
  mp
  if(plotly == T) {
    mpp <- plotly::ggplotly(mp)
    htmlwidgets::saveWidget(plotly::as_widget(mp),
                            plotly_filename,
                            knitrOptions = list(fig.width = 10, fig.height = 10),
                            selfcontained = T)
  }
  mp
}
