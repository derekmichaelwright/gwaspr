#' gg_Marker
#'
#' Creates a marker plot with myG and myY objects.
#' @param myG GWAS genotype object. Note: needs to be in hapmap format.
#' @param myY GWAS phenotype object.
#' @param trait Trait to plot.
#' @param markers Markers to plot.
#' @param points Logical, whether or not to plot points
#' @param colors Color palette.
#' @return Marker plot.
#' @export

gg_Marker_Box <- function (myG, myY, trait,
                           markers,
                           points = T,
                           colors = c("darkgreen", "darkgoldenrod3", "darkred", "steelblue4",
                                      "darkslategray", "maroon4", "purple4", "darkblue"))
  {
  #markers = c("Lcu.2RBY.Chr2p42543877", "Lcu.2RBY.Chr5p1069654")
  title <- paste(markers, collapse = "-")
  #
  xx <- myG %>% filter(rs %in% markers) %>%
    select(-2,-3,-4,-5,-6,-7,-8,-9,-10,-11) %>%
    column_to_rownames("rs") %>%
    t() %>% as.data.frame() %>% mutate(Alleles = NA)
  for(i in 1:length(markers)) { xx <- xx[xx[,i] %in% c("A","T","G","C","AA","TT","GG","CC"),]}
  #i<-1
  for(i in 1:nrow(xx)) { xx$Alleles[i] <- paste(xx[i,1:length(markers)], collapse = "-") }
  xx <- xx %>% rownames_to_column("Name") %>%
    left_join(myY, by = "Name") %>%
    filter(!is.na(get(trait)))
  # Plot
  mp <- ggplot(xx, aes(x = Alleles, y = get(trait))) +
    geom_violin(aes(fill = Alleles), alpha = 0.5) +
    geom_boxplot(width = 0.1, outlier.shape = NA) +
    scale_fill_manual(name = NULL, values = colors) +
    theme_gwaspr(legend.position = "none") +
    labs(y = trait, x = title)
    if (points == T) { mp <- mp + geom_quasirandom(alpha = 0.5, pch = 16) }
  mp
}
