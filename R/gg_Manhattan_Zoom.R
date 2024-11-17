#' gg_Manhattan_Zoom
#'
#' Creates a manhattan plot zoomed in to a particular region.
#' @param folder Folder containing GWAS results.
#' @param trait The trait to read.
#' @param chr Chromosome to plot.
#' @param start Start position on chromosome.
#' @param end End position on chromosome.
#' @param title A title for the plot.
#' @param markers Markers to be labelled.
#' @param labels Labels to be used for markers.
#' @param vlines Markers which will be used as a location for a vertical lines.
#' @param vline.colors colors for each vertical line.
#' @param models Models to read.
#' @return A manhattan plot.
#' @export

gg_Manhattan_Zoom <- function(folder, trait, chr, start, end,
                              title = trait,
                              markers = NULL, labels = markers,
                              vlines = markers, vline.colors = "red",
                              models = c("GLM","MLM","CMLM","MLMM","SUPER","FarmCPU","BLINK") ) {
  fnames <- list.files(folder)[grepl("GWAS_Results", list.files(folder))]
  fnames <- fnames[grepl(paste0(trait,".csv"), fnames)]
  xx <- NULL
  # i <- fnames[3]
  for(i in fnames) {
    mod <- substr(i, gregexpr("GWAS_Results", i)[[1]][1]+13, gregexpr(".csv", i)[[1]][1]-1)
    mod <- substr(mod, 1, gregexpr("\\.", mod)[[1]][1]-1)
    xi <- read.csv(paste0(folder, i))
    if(sum(colnames(xi)=="nobs")>0) { xi <- select(xi, -nobs) }
    xi <- xi %>%
      mutate(Model = mod,
             `-log10(p)`     = -log10(P.value),
             `-log10(p)_exp` = -log10((rank(P.value, ties.method="first")-.5)/nrow(.)))
    xx <- bind_rows(xx, xi)
  }
  xx <- xx %>% mutate(Model = factor(Model, levels = models)) %>%
    filter(Chr == chr, Pos > start, Pos < end, !is.na(Model))
  threshold <- -log10(0.05 / nrow(xi))
  x1 <- xx %>% filter(`-log10(p)` < threshold)
  x2 <- xx %>% filter(`-log10(p)` > threshold)
  # Man plot
  mp <- ggplot(xx, aes(x = Pos / 1000000, y = `-log10(p)`))
  if(!is.null(vlines)) {
    mp <- mp +
      geom_vline(data = xx %>% filter(SNP %in% markers),
                 aes(xintercept = Pos / 1000000, color = SNP),
                 alpha = 0.5) +
      scale_color_manual(values = vline.colors)
  }
  mp <- mp +
    geom_hline(yintercept = threshold, color = "red", alpha = 0.6) +
    geom_point(alpha = 0.8, color = "darkgreen", pch = 16) +
    geom_point(data = x2, pch = 16, size = 1.5, color = "darkred", alpha = 0.8) +
    facet_grid(Model ~ Chr, scales = "free") +

    theme_gwaspr(legend.position = "none",
                 axis.title.y = element_markdown()) +
    labs(title = trait, y = "-log<sub>10</sub>(*p*)", x = "Mbp")
  if(!is.null(markers)) {
    xx <- xx %>% mutate(Label = ifelse(SNP %in% markers, SNP, NA),
                        Label = plyr::mapvalues(Label, markers, labels))
    mp <- mp + geom_text_repel(data = xx %>% filter(SNP %in% markers),
                               aes(label = Label), size = 2)
  }
  mp
}
