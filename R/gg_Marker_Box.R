#' gg_Marker
#'
#' Creates a marker plot with myG and myY objects.
#' @param xG GWAS genotype object. Note: needs to be in hapmap format.
#' @param xY GWAS phenotype object.
#' @param traits Traits to plot.
#' @param markers Markers to plot.
#' @param box.width width for the boxplot.
#' @param violin Logical, whether or not to plot violin.
#' @param points Logical, whether or not to plot points.
#' @param marker.colors Color palette.
#' @return Marker plot.
#' @export

gg_Marker_Box <- function (
    xG,
    xY,
    traits,
    myncol = NULL,
    markers,
    box.width = 0.1,
    points = T,
    violin = T,
    marker.colors = c("darkgreen", "darkgoldenrod3", "darkred", "steelblue4",
                      "darkslategray", "maroon4", "purple4", "darkblue") ) {
  #
  title <- paste(markers, collapse = "\n")
  #
  xY <- xY %>% select(1, traits) %>%
    gather(Trait, Value, traits)
  #
  xx <- xG %>% rename(SNP=1) %>%
    filter(SNP %in% markers) %>%
    select(-2,-3,-4,-5,-6,-7,-8,-9,-10,-11) %>%
    column_to_rownames("SNP") %>%
    t() %>% as.data.frame() %>% mutate(Alleles = NA)
  #
  for(i in 1:length(markers)) { xx <- xx[xx[,i] %in% c("A","T","G","C","AA","TT","GG","CC"),] }
  #
  for(i in 1:nrow(xx)) { xx$Alleles[i] <- paste(xx[i,1:length(markers)], collapse = "-") }
  #
  xx <- xx %>% rownames_to_column("Name") %>%
    left_join(xY, by = "Name") %>%
    filter(!is.na(Value))
  #
  yy <- xx %>% filter(Trait == traits[1]) %>% 
    group_by(Alleles) %>%
    summarise(Value = mean(Value, na.rm = T))
  xx <- xx %>% 
    mutate(Alleles = factor(Alleles, levels = rev(yy$Alleles)),
           Trait = factor(Trait, levels = traits))
  # Plot
  mp <- ggplot(xx, aes(x = Alleles, y = Value))
  if(violin == T) { mp <- mp + geom_violin(aes(fill = Alleles), alpha = 0.5) }
  mp <- mp +
    geom_boxplot(width = box.width, outlier.shape = NA) +
    facet_wrap(Trait ~ ., scales = "free_y", ncol = myncol) +
    scale_fill_manual(name = NULL, values = marker.colors) +
    theme_gwaspr(legend.position = "none",
                 axis.text.x = element_text(angle = 45, hjust = 1) ) +
    labs(title = title, y = NULL, x = NULL)
  if (points == T) { mp <- mp + geom_quasirandom(alpha = 0.5, pch = 16) }
  mp
}

#xG = myG; xY = myY; 
#traits = c("Canopy.Height_Ba16", "Canopy.Width_Ba16", "Canopy.Height_Ba17", "Canopy.Width_Ba17")
#markers = myMarkers
#box.width = 0.1; points = T; violin = T; myncol = 4
#marker.colors = c("darkgreen", "darkgoldenrod3", "darkred", "steelblue4", "darkslategray", "maroon4", "purple4", "darkblue")
