#' gg_Marker_Bar
#'
#' Creates a marker plot with myG and myY objects.
#' @param xG GWAS genotype object. Note:  needs to be in hapmap format.
#' @param xY GWAS phenotype object.
#' @param traits Traits to plot.
#' @param markers Markers to plot.
#' @param marker.colors Color palette.
#' @param plot.histogram Logical, if true will plot histogram bars.
#' @param plot.density Logical, if true will plot density bands.
#' @param plot.counts Logical, if true will make a plot of counts, if false will make a density plot.
#' @param myncol Number of columns for facetting when plotting multiple traits.
#' @return Marker plot.
#' @export

gg_Marker_Bar <- function (
    xG,
    xY,
    traits,
    markers,
    marker.colors = gwaspr_Colors,
    plot.histogram = T,
    plot.density = T,
    plot.counts = F,
    myncol = NULL,
    line.color = F
    ) {
  #
  title <- paste(markers, collapse = "\n")
  #
  xY <- xY %>% dplyr::select(1, traits) %>%
    gather(Trait, Value, traits)
  #
  xx <- xG %>% rename(SNP=1) %>%
    filter(SNP %in% markers) %>%
    dplyr::select(-2,-3,-4,-5,-6,-7,-8,-9,-10,-11) %>%
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
    summarise(Value = mean(Value, na.rm = T)) %>%
    arrange(Value)
  xx <- xx %>%
    mutate(Alleles = factor(Alleles, levels = rev(yy$Alleles)),
           Trait = factor(Trait, levels = traits))
  # Plot
  if(plot.counts == F) { mp <- ggplot(xx, aes(x = Value, y=after_stat(density), fill = Alleles)) }
  if(plot.counts == T) { mp <- ggplot(xx, aes(x = Value, y=after_stat(count), fill = Alleles)) }
  if(plot.density == T) { mp <- mp + geom_density(alpha = 0.5, color = line.color) }
  if(plot.histogram == T) { mp <- mp + geom_histogram(position = "dodge", alpha = 0.5, color = line.color) }
  mp <- mp +
    facet_wrap(Trait ~ ., scales = "free", ncol = myncol) +
    scale_fill_manual(name = NULL, values = marker.colors) +
    theme_gwaspr_col(legend.position = "bottom"
                     ) +
    labs(title = title, x = NULL)
  mp
}

#xG = myG; xY = myY;  markers = myMarkers[1]
#traits = c("Disease.Score_Ba16","Lodging.Score_Ba16", "Stem.Blight_Ba17", "Stem.Blight_Ba16")
#marker.colors = c("darkgreen", "darkgoldenrod3", "darkred", "steelblue4", "darkslategray", "maroon4", "purple4", "darkblue")
#myncol = NULL; plot.histogram = T; plot.density = T
