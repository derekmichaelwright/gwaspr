#' gg_Manh_Zoom
#'
#' Removes columns without any values.
#' @param folder Folder containing GWAS results.
#' @param trait The trait to read.
#' @param subtitle A subtitle for the plot.
#' @param markers Markers to be labelled.
#' @param labels Labels to be used for markers.
#' @param lines Logical value of whether or not to include vertical lines with markers.
#' @param models Models to read.
#' @param colors Colors for each chromosome
#' @return A manhattan plot.
#' @export

gg_Manh_Zoom <- function(folder = "Results_Add/", 
                        trait = "Testa_Pattern",
                        chr = 1, start = 100000, end = 1000000,
                        markers = NULL, labels = markers, lines = F,
                        models = c("GLM","MLM","CMLM","MLMM","SUPER","FarmCPU","Blink"), 
                        colors = c("darkgreen","darkgoldenrod3","darkgreen","darkgoldenrod3",
                                   "darkgreen", "darkgoldenrod3","darkgreen")) {
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
             negLog10     = -log10(P.value),
             negLog10_exp = -log10((rank(P.value, ties.method="first")-.5)/nrow(.)))
    xx <- bind_rows(xx, xi)
  }
  xx <- xx %>% mutate(Model = factor(Model, levels = models)) %>%
    filter(Chromosome == chr, Position > start, Position < end, !is.na(Model))
  threshold <- -log10(0.05 / nrow(xi))
  x1 <- xx %>% filter(negLog10 < threshold)
  x2 <- xx %>% filter(negLog10 > threshold)
  # Man plot
  mp <- ggplot(xx, aes(x = Position / 1000000, y = -log10(P.value)))
  if(!is.null(markers) & lines == T) {
    mp <- mp + 
      geom_vline(data = xx %>% filter(SNP %in% markers), 
                 aes(xintercept = Position / 1000000), alpha = 0.5)
  }
  mp <- mp +
    geom_hline(yintercept = threshold, alpha = 0.6) +
    geom_point(aes(color = factor(Chromosome)), pch = 1, alpha = 0.8) +
    geom_line() +
    geom_point(data = x2, pch = 21, size = 1.5, color = "black", fill = "darkred", alpha = 0.8) +
    facet_grid(Model ~ Chromosome, scales = "free") +
    #scale_x_continuous(breaks = seq(100, 700, by = 100)) +
    scale_color_manual(values = colors) +
    theme_gwaspr(legend.position = "none",
                 axis.text.x = element_text(angle = 90, hjust = 0.5)) +
    labs(title = trait, y = "-log10(p)", x = "Mbp")
  if(!is.null(markers)) {
    xx <- xx %>% mutate(Label = ifelse(SNP %in% markers, SNP, NA),
                        Label = plyr::mapvalues(Label, markers, labels))
    mp <- mp + geom_text_repel(data = xx %>% filter(SNP %in% markers),
                               aes(label = Label), size = 2)
  }
  mp
}
