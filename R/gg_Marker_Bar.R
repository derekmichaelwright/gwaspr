#' gg_Marker
#'
#' Creates a marker plot with myG and myY objects.
#' @param myG GWAS genotype object. Note:  needs to be in hapmap format.
#' @param myY GWAS phenotype object.
#' @param trait Trait to plot.
#' @param markers Markers to plot.
#' @param colors Color palette.
#' @return Marker plot.
#' @export

gg_Marker_Bar <- function (myG, myY, trait,
                           markers,
                           colors = c("darkgreen", "darkgoldenrod3", "darkred", "steelblue4",
                                      "darkslategray", "maroon4", "purple4", "darkblue"))
  {

 markers <- c("Lcu.2RBY.Chr6p12212845", "Lcu.2RBY.Chr6p12192948", "Lcu.2RBY.Chr6p14410759")
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
 xx[,trait] <- factor(xx[,trait])
 # Plot
 mp <- ggplot(xx, aes(x = Alleles)) +
   geom_bar(aes(fill = get(trait)), position = "dodge", color = "black", alpha = 0.5) +
   scale_fill_manual(name = NULL, values = colors) +
   theme_gwaspr(legend.position = "bottom") +
   labs(x = title)
 mp
}
