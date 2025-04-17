#' gg_Marker
#'
#' Creates a marker plot with myG and myY objects.
#' @param xG GWAS genotype object. Note:  needs to be in hapmap format.
#' @param xY GWAS phenotype object.
#' @param trait Trait to plot.
#' @param markers Markers to plot.
#' @param marker.colors Color palette.
#' @return Marker plot.
#' @export

gg_Marker_Bar <- function (
    xG,
    xY,
    trait,
    markers,
    marker.colors = c("darkgreen", "darkgoldenrod3", "darkred", "steelblue4",
                      "darkslategray", "maroon4", "purple4", "darkblue") ) {
 #
 title <- paste(markers, collapse = "\n")
 #
 xx <- xG %>% rename(SNP=1) %>%
   filter(SNP %in% markers) %>%
   select(-2,-3,-4,-5,-6,-7,-8,-9,-10,-11) %>%
   column_to_rownames("SNP") %>%
   t() %>% as.data.frame() %>% mutate(Alleles = NA)
 for(i in 1:length(markers)) { xx <- xx[xx[,i] %in% c("A","T","G","C","AA","TT","GG","CC"),]}
 #
 for(i in 1:nrow(xx)) { xx$Alleles[i] <- paste(xx[i,1:length(markers)], collapse = "-") }
 xx <- xx %>% rownames_to_column("Name") %>%
   left_join(xY, by = "Name") %>%
   filter(!is.na(get(trait)))
 xx[,trait] <- factor(xx[,trait])
 # Plot
 mp <- ggplot(xx, aes(x = Alleles)) +
   geom_bar(aes(fill = get(trait)), position = "dodge", color = "black", alpha = 0.5) +
   scale_fill_manual(name = NULL, values = marker.colors) +
   theme_gwaspr(legend.position = "bottom") +
   labs(title = title, y = trait, x = NULL)
 mp
}
