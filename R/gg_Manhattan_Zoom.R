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
#' @param lines Logical value of whether or not to include vertical lines with markers.
#' @param models Models to read.
#' @param colors Colors for each chromosome
#' @return A manhattan plot.
#' @export

gg_Manhattan_Zoom <- function(folder, trait, chr, start, end,
                              title = trait,
                              markers = NULL, labels = markers, lines = F,
                              models = c("GLM","MLM","CMLM","MLMM","SUPER","FarmCPU","Blink"),
                              colors = c("darkgreen","darkgoldenrod3","darkgreen","darkgoldenrod3",
                                         "darkgreen", "darkgoldenrod3","darkgreen") ) {
  fnames <- grep(paste0(trait,".GWAS.Results"), list.files(folder))
  fnames <- list.files(folder)[fnames]
  xx <- NULL
  # i <- fnames[3]
  for(i in fnames) {
    mod <- substr(i, gregexpr("\\.", i)[[1]][1]+1, gregexpr("\\.", i)[[1]][2]-1)
    xi <- read.csv(paste0(folder, i))
    if(sum(colnames(xi)=="nobs")>0) { xi <- select(xi, -nobs) }
    xi <- xi %>%
      mutate(Model = mod,
             `-log10(p)`     = -log10(P.value),
             `-log10(p)_exp` = -log10((rank(P.value, ties.method="first")-.5)/nrow(.)))
    xx <- bind_rows(xx, xi)
  }
  xx <- xx %>% mutate(Model = factor(Model, levels = models)) %>%
    filter(Chromosome == chr, Position > start, Position < end, !is.na(Model))
  threshold <- -log10(0.05 / nrow(xi))
  x1 <- xx %>% filter(`-log10(p)` < threshold)
  x2 <- xx %>% filter(`-log10(p)` > threshold)
  # Man plot
  mp <- ggplot(xx, aes(x = Position / 1000000, y = `-log10(p)`))
  if(!is.null(markers) & lines == T) {
    mp <- mp +
      geom_vline(data = xx %>% filter(SNP %in% markers),
                 aes(xintercept = Position / 1000000), alpha = 0.5)
  }
  if(lines == T) { mp <- mp + geom_line() }
  mp <- mp +
    geom_hline(yintercept = threshold, alpha = 0.6) +
    geom_point(aes(color = factor(Chromosome)), pch = 1, alpha = 0.8) +
    geom_point(data = x2, pch = 21, size = 1.5, color = "black", fill = "darkred", alpha = 0.8) +
    facet_grid(Model ~ Chromosome, scales = "free") +
    scale_color_manual(values = colors) +
    theme_gwaspr(legend.position = "none",
                 #axis.text.x = element_text(angle = 90, hjust = 0.5),
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
