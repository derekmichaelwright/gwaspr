#' gg_Marker_Scatter
#'
#' [Creates a marker scatter plot with myG and myY objects.](https://derekmichaelwright.github.io/gwaspr/articles/gg_Marker.html)
#' @param xG GWAS genotype object. Note: needs to be in hapmap format.
#' @param xY GWAS phenotype object.
#' @param trait.x Trait for the x-axis.
#' @param trait.y Trait for the y-axis.
#' @param markers Markers to plot.
#' @param marker.colors Colors to fill in the violin and boxplots.
#' @param remove.hets Logical, Whether to remove hets or not. advised if plotting multiple markers.
#' @param point.size Size for the points.
#' @param title Title for the plot.
#' @param legend.rows Number of rows for the legend.
#' @param subtitle Subtitle for the plot. Defaults to the list of markers.
#' @param yLab Label for the y-axis.
#' @return Marker plot.
#' @export

gg_Marker_Scatter <- function (
    xG,
    xY,
    trait.x,
    trait.y,
    markers,
    marker.colors = gwaspr_Colors,
    remove.hets = T,
    point.size = 1,
    title = NULL,
    legend.rows = 1,
    subtitle = paste(markers, collapse = "\n")
    ) {
  #
  myLab <- paste(markers, collapse = "\n")
  #
  xT <- xY %>% select(Name=1, trait.x, trait.y)# %>%
    #gather(Trait, Value, trait.x, trait.y)
  #
  xx <- xG %>% rename(SNP=1) %>%
    filter(SNP %in% markers) %>%
    dplyr::select(-2,-3,-4,-5,-6,-7,-8,-9,-10,-11) %>%
    column_to_rownames("SNP") %>%
    t() %>% as.data.frame() %>%
    select(markers) %>%
    mutate(Alleles = NA)
  #
  if(remove.hets == T) { for(i in 1:length(markers)) { xx <- xx[xx[,i] %in% c("A","T","G","C","AA","TT","GG","CC"),] } }
  if(remove.hets == F) { for(i in 1:length(markers)) { xx <- xx[!xx[,i] %in% c("N","NN"),] } }
  #
  for(i in 1:nrow(xx)) { xx$Alleles[i] <- paste(xx[i,1:length(markers)], collapse = "-") }
  #
  xx <- xx %>% rownames_to_column("Name") %>%
    left_join(xT, by = "Name") %>%
    filter(!is.na(get(trait.x))) %>%
    filter(!is.na(get(trait.y)))
  #
  #yy <- xx %>% filter(Trait == traits[1]) %>%
  #  group_by(Alleles) %>%
  #  summarise(Value = mean(Value, na.rm = T)) %>%
  #  arrange(Value)
  xx <- xx %>%
    mutate(Alleles = factor(Alleles))#, #, levels = rev(yy$Alleles)),
           #Trait = factor(Trait, levels = traits))
  # Plot
  mp <- ggplot(xx, aes(x = get(trait.x), y = get(trait.y))) +
    geom_point(aes(color = Alleles), pch = 16, alpha = 0.7) +
    scale_color_manual(values = marker.colors) +
    theme_gwaspr(legend.position = "bottom",
                 axis.text.x = element_text(angle = 45, hjust = 1) ) +
    #guides(color = guide_legend(nrow = legend.rows, override.aes = list(size = 2))) +
    labs(title = title, subtitle = subtitle, x = trait.x, y = trait.y)
  mp
}

#xG = myG_LDP; xY = myY_LDP_BELT_blups;
#trait.x = "BELT.height.mm._LDP"; trait.y = "BELT.diameter.mm._LDP"
#markers = c("Lcu.1GRN.Chr1p371692181","Lcu.1GRN.Chr7p593620255")
#marker.colors = gwaspr_Colors
#remove.hets = T; point.size = 1; title = NULL; legend.rows = 1
#subtitle = paste(markers, collapse = "\n")
